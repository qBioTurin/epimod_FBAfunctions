# ---------------------------------------------------------------------
# Validate and repair a PNPRO file against a set of FBA models — DEBUG
# ---------------------------------------------------------------------
#' @export
validate_pnpro <- function(pnpro2validate,
                           hypernode_root,
                           biounit_models,
                           boundary_metabolites,
                           out_dir,
                           hypernode_name) {

  message("\n─── 0) read & parse PNPRO ─────────────────────────────────")
  message("file: ", pnpro2validate)

  ## 0) Read & parse the PNPRO -------------------------------------------------
  xml <- read_xml(pnpro2validate)

  tnodes   <- xml_find_all(xml, "//transition")
  t_names  <- xml_attr(tnodes, "name")
  t_delays <- xml_attr(tnodes, "delay")

  is_call  <- str_detect(t_delays, "^Call\\[")
  is_fba   <- str_detect(t_delays, "^FBA\\[")
  special  <- t_names[is_call | is_fba]

  arc_df <- xml_find_all(xml, "//arc") %>%
    map_df(~{
      a <- xml_attrs(.x)
      tibble(
        head         = a["head"],
        tail         = a["tail"],
        kind         = a["kind"],
        multiplicity = as.integer(coalesce(a["mult"], "1"))
      )
    }) %>%
    filter((kind == "INPUT"  & head %in% special) |
           (kind == "OUTPUT" & tail %in% special)) %>%
    transmute(
      transition = if_else(kind == "INPUT", head, tail),
      direction  = kind,
      place      = if_else(kind == "INPUT", tail, head),
      multiplicity,
      command    = t_delays[match(transition, t_names)]
    )

  message("✔ arcs parsed: ", nrow(arc_df))

  ## 1) Re-derive boundary reactions ------------------------------------------
  message("\n─── 1) scan FBA metadata ──────────────────────────────────")

  abbrs <- map_chr(biounit_models, ~ .x$abbreviation[2])   # original line

  models_df <- tibble(model = biounit_models) %>%
    mutate(
      abbr     = map_chr(model, ~ .x$abbreviation[2]),
      meta_dir = file.path(hypernode_root, "biounits",
                           map_chr(model, ~ .x$FBAmodel))
    ) %>%
    mutate(
      metabolites = map(meta_dir, ~ read_csv(
        file.path(.x, "metabolites_metadata.csv"), show_col_types = FALSE)),
      reactions   = map(meta_dir, ~ read_csv(
        file.path(.x, "reactions_metadata.csv"),   show_col_types = FALSE))
    )

  message("✔ biounits loaded: ", nrow(models_df))

  projectable_df <- models_df %>%
    unnest(metabolites) %>%
    filter(id %in% boundary_metabolites) %>%
    select(abbr, met_id = id)

  boundary_df <- models_df %>%
    unnest(reactions) %>%
    filter(type == "boundary") %>%
    select(abbr, reaction = abbreviation, equation)

  shared_rxns_df <- projectable_df %>%
    left_join(boundary_df, by = "abbr", relationship = "many-to-many") %>%
    filter(str_detect(equation, paste0("\\b", met_id, "\\b"))) %>%
    distinct(abbr, met_id, reaction)

  message("✔ projectable reactions: ", nrow(shared_rxns_df))

  ## 2) Parse (original) FBA / Call commands ----------------------------------
  message("\n─── 2) parse original PN commands ─────────────────────────")

  fba_cmds <- tibble(
    transition = t_names,
    delay      = t_delays
  ) %>%
    filter(str_detect(delay, "^FBA\\[")) %>%
    mutate(
      parts    = str_match(
        delay,
        'FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9\\.]+) *, *"([^"]+)" *, *"([^"]+)"'),
      model    = parts[, 2],
      reaction = parts[, 3],
      abbr     = str_remove(parts[, 5], "^n_")
    ) %>%
    select(transition, delay, reaction, abbr)

  call_cmds <- tibble(
    transition = t_names,
    delay      = t_delays
  ) %>%
    filter(str_detect(delay, "^Call\\[")) %>%
    mutate(
      parts      = str_match(
        delay,
        'Call\\[\\s*"([^"]+)"\\s*,\\s*(.+)\\s*,\\s*([0-9]+)\\s*\\]'),
      fun_name   = parts[, 2],
      param_expr = parts[, 3],
      org_index  = as.integer(parts[, 4])
    ) %>%
    select(transition, command = delay, fun_name, param_expr, org_index)

  message("✔ original FBA cmds : ", nrow(fba_cmds))
  message("✔ original Call cmds: ", nrow(call_cmds))

  ## 2a) Build repaired FBA commands ------------------------------------------
  biomass_rxns <- tibble(
    abbr     = abbrs,
    met_id   = "biomass_e",
    reaction = "EX_biomass_e"
  )

  repaired_fba_cmds <- bind_rows(shared_rxns_df, biomass_rxns) %>%
    mutate(
      model_file    = paste0(abbr, "_model.txt"),
      count_place   = paste0("n_",         abbr),
      biomass_place = paste0("biomass_e_", abbr),
      scaling       = 1L,
      is_biomass    = reaction == "EX_biomass_e"
    ) %>%
    crossing(direction = c("INPUT", "OUTPUT")) %>%
    mutate(
      transition = paste0(
        reaction,
        if_else(direction == "INPUT", "_in_", "_out_"), abbr),
      place   = if_else(reaction == "EX_biomass_e", biomass_place, met_id),
      command = sprintf(
        'FBA["%s","%s",%d,"%s","%s"%s]',
        model_file, reaction, scaling, count_place, biomass_place,
        if_else(is_biomass, ', "true"', ''))
    ) %>%
    select(transition, command, reaction, abbr, direction, place)

  message("✔ repaired FBA cmds : ", nrow(repaired_fba_cmds))

  ## 2b) Build repaired Call commands -----------------------------------------
  func_cols  <- c(Starvation = 0L, Duplication = 1L, Death = 2L)
  prefix_map <- c(Starvation = "Starv", Duplication = "Dup", Death = "Death")

  repaired_call_cmds <- imap_dfr(abbrs, function(abbr, idx) {
    org_index <- idx - 1L
    tibble(
      fun_name  = names(func_cols),
      col_index = unname(func_cols),
      org_index = org_index
    ) %>%
      mutate(
        transition = paste0(prefix_map[fun_name], "_", abbr),
        command    = sprintf(
          'Call["%s", FromTable["population_parameters.csv", %d, %d], %d]',
          fun_name, org_index, col_index, org_index))
  })

  message("✔ repaired Call cmds: ", nrow(repaired_call_cmds))

  ## 3) Build call_arcs --------------------------------------------------------
  call_arcs <- purrr::pmap_dfr(
    list(
      transition = repaired_call_cmds$transition,
      command    = repaired_call_cmds$command,
      fun_name   = repaired_call_cmds$fun_name,
      org_index  = repaired_call_cmds$org_index
    ),
    function(transition, command, fun_name, org_index) {

      abbr          <- abbrs[org_index + 1]
      count_place   <- paste0("n_",         abbr)
      biomass_place <- paste0("biomass_e_", abbr)

      switch(fun_name,
        "Starvation" =
          tibble(transition = transition, direction = "INPUT",
                 place = biomass_place, multiplicity = 1L, command = command),

        "Duplication" =
          bind_rows(
            tibble(transition = transition, direction = "INPUT",
                   place = count_place,   multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "INPUT",
                   place = biomass_place, multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = count_place,   multiplicity = 2L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = biomass_place, multiplicity = 1L, command = command)),

        "Death" =
          bind_rows(
            tibble(transition = transition, direction = "INPUT",
                   place = count_place,   multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "INPUT",
                   place = biomass_place, multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = biomass_place, multiplicity = 1L, command = command))
      )
    })

  ## 4) Build FBA arcs ---------------------------------------------------------
  fba_arcs <- purrr::pmap_dfr(
    list(
      transition = repaired_fba_cmds$transition,
      reaction   = repaired_fba_cmds$reaction,
      abbr       = repaired_fba_cmds$abbr,
      met_place  = repaired_fba_cmds$place,
      command    = repaired_fba_cmds$command
    ),
    function(transition, reaction, abbr, met_place, command) {

      count_place   <- paste0("n_",         abbr)
      biomass_place <- paste0("biomass_e_", abbr)

      if (reaction != "EX_biomass_e") {

        if (endsWith(transition, paste0("_in_", abbr))) {
          # IMPORT
          bind_rows(
            tibble(transition = transition, direction = "INPUT",
                   place = met_place,     multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "INPUT",
                   place = count_place,   multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "INPUT",
                   place = biomass_place, multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = count_place,   multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = biomass_place, multiplicity = 1L, command = command))

        } else {  # EXPORT
          bind_rows(
            tibble(transition = transition, direction = "INPUT",
                   place = count_place,   multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "INPUT",
                   place = biomass_place, multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = met_place,     multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = count_place,   multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = biomass_place, multiplicity = 1L, command = command))
        }

      } else {  # biomass boundary --------------------------------------------
        if (endsWith(transition, paste0("_in_", abbr))) {
          tibble(transition = transition, direction = "INPUT",
                 place = biomass_place, multiplicity = 1L, command = command)
        } else {           # OUT
          bind_rows(
            tibble(transition = transition, direction = "INPUT",
                   place = biomass_place, multiplicity = 1L, command = command),
            tibble(transition = transition, direction = "OUTPUT",
                   place = biomass_place, multiplicity = 2L, command = command))
        }
      }
    })

  ## 5) Combine, write, return -------------------------------------------------
  arc_df_repaired <- bind_rows(call_arcs, fba_arcs) %>%
    distinct(transition, direction, place, multiplicity, command)

  message("✔ repaired arcs: ", nrow(arc_df_repaired))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(arc_df,          file.path(out_dir, "raw_arcs.csv"))
  write_csv(arc_df_repaired, file.path(out_dir, "repaired_arcs.csv"))
  message("✔ CSVs written in ", out_dir)

  invisible(list(
    raw_arcs      = arc_df,
    repaired_arcs = arc_df_repaired))
}

