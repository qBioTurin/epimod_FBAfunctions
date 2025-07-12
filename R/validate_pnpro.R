#' Validate and repair a PNPRO file against a set of FBA models
#'
#' @export
validate_pnpro <- function(pnpro2validate,
                           hypernode_root,
                           biounit_models,
                           boundary_metabolites,
                           out_dir,
                           hypernode_name) {
  # ------------------------------------------------------------
  # Helper vectors
  # ------------------------------------------------------------
  abbrs <- purrr::map_chr(biounit_models, ~ .x$abbreviation[2])

  # ------------------------------------------------------------
  # 0) Read & parse the PNPRO
  # ------------------------------------------------------------
  xml <- xml2::read_xml(pnpro2validate)

  tnodes   <- xml2::xml_find_all(xml, "//transition")
  t_names  <- xml2::xml_attr(tnodes, "name")
  t_delays <- xml2::xml_attr(tnodes, "delay")

  is_call <- stringr::str_detect(t_delays, "^Call\\[")
  is_fba  <- stringr::str_detect(t_delays, "^FBA\\[")
  special <- t_names[is_call | is_fba]

  arc_df <- xml2::xml_find_all(xml, "//arc") %>%
    purrr::map_df(~{
      a <- xml2::xml_attrs(.x)
      tibble::tibble(
        head         = a["head"],
        tail         = a["tail"],
        kind         = a["kind"],
        multiplicity = as.integer(dplyr::coalesce(a["mult"], "1"))
      )
    }) %>%
    dplyr::filter((kind == "INPUT"  & head %in% special) |
                  (kind == "OUTPUT" & tail %in% special)) %>%
    dplyr::transmute(
      transition = dplyr::if_else(kind == "INPUT", head, tail),
      direction  = kind,
      place      = dplyr::if_else(kind == "INPUT", tail, head),
      multiplicity,
      command    = t_delays[match(transition, t_names)]
    )

  # ------------------------------------------------------------
  # 1) Determine the boundary reactions each species should expose
  # ------------------------------------------------------------
  models_df <- tibble::tibble(model = biounit_models) %>%
    dplyr::mutate(
      abbr     = purrr::map_chr(model, ~ .x$abbreviation[2]),
      meta_dir = file.path(hypernode_root, "biounits", purrr::map_chr(model, ~ .x$FBAmodel))
    ) %>%
    dplyr::mutate(
      metabolites = purrr::map(meta_dir, ~ readr::read_csv(file.path(.x, "metabolites_metadata.csv"), show_col_types = FALSE)),
      reactions   = purrr::map(meta_dir, ~ readr::read_csv(file.path(.x, "reactions_metadata.csv"),   show_col_types = FALSE))
    )

  projectable_df <- models_df %>%
    tidyr::unnest(metabolites) %>%
    dplyr::filter(id %in% boundary_metabolites) %>%
    dplyr::select(abbr, met_id = id)

  boundary_df <- models_df %>%
    tidyr::unnest(reactions) %>%
    dplyr::filter(type == "boundary") %>%
    dplyr::select(abbr, reaction = abbreviation, equation)

  shared_rxns_df <- projectable_df %>%
    dplyr::left_join(boundary_df, by = "abbr", relationship = "many-to-many") %>%
    dplyr::filter(stringr::str_detect(equation, paste0("\\b", met_id, "\\b"))) %>%
    dplyr::distinct(abbr, met_id, reaction)

  # ------------------------------------------------------------
  # 2) Parse Call[…] and FBA[…] commands from PN
  # ------------------------------------------------------------
  fba_cmds <- tibble::tibble(
    transition = t_names,
    delay      = t_delays
  ) %>%
    dplyr::filter(stringr::str_detect(delay, "^FBA\\[")) %>%
    dplyr::mutate(
      parts    = stringr::str_match(delay,
                  'FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9\\.]+) *, *"([^"]+)" *, *"([^"]+)"'),
      model    = parts[, 2],
      reaction = parts[, 3],
      abbr     = stringr::str_remove(parts[, 5], "^n_")
    ) %>%
    dplyr::select(transition, delay, reaction, abbr)

  call_cmds <- tibble::tibble(
    transition = t_names,
    delay      = t_delays
  ) %>%
    dplyr::filter(stringr::str_detect(delay, "^Call\\[")) %>%
    dplyr::mutate(
      parts      = stringr::str_match(delay, 'Call\\[\\s*"([^"]+)"\\s*,\\s*(.+)\\s*,\\s*([0-9]+)\\s*\\]'),
      fun_name   = parts[, 2],
      param_expr = parts[, 3],
      org_index  = as.integer(parts[, 4])
    ) %>%
    dplyr::select(transition, command = delay, fun_name, param_expr, org_index)

  # Add biomass reaction for each species
  biomass_rxns <- tibble::tibble(
    abbr     = abbrs,
    met_id   = "biomass_e",
    reaction = "EX_biomass_e"
  )

  shared_rxns_df <- dplyr::bind_rows(shared_rxns_df, biomass_rxns)

  repaired_fba_cmds <- shared_rxns_df %>%
    dplyr::mutate(
      model_file    = paste0(abbr, "_model.txt"),
      count_place   = paste0("n_",         abbr),
      biomass_place = paste0("biomass_e_", abbr),
      scaling       = 1L,
      is_biomass    = reaction == "EX_biomass_e"
    ) %>%
    tidyr::crossing(direction = c("INPUT", "OUTPUT")) %>%
    dplyr::mutate(
      transition = paste0(
        reaction,
        dplyr::if_else(direction == "INPUT", "_in_", "_out_"),
        abbr
      ),
      place = dplyr::if_else(reaction == "EX_biomass_e", biomass_place, met_id),
      command = sprintf(
        'FBA["%s","%s",%d,"%s","%s"%s]',
        model_file, reaction, scaling, count_place, biomass_place,
        dplyr::if_else(is_biomass, ', "true"', '')
      )
    ) %>%
    dplyr::select(transition, command, reaction, abbr, direction, place)

  # ------------------------------------------------------------
  # 2b) Repair Call commands (population dynamics)
  # ------------------------------------------------------------
  func_cols  <- c(Starvation = 0L, Duplication = 1L, Death = 2L)
  prefix_map <- c(Starvation = "Starv", Duplication = "Dup", Death = "Death")

  repaired_call_cmds <- purrr::imap_dfr(abbrs, function(abbr, idx) {
    org_index <- idx - 1L

    tibble::tibble(
      fun_name  = names(func_cols),
      col_index = unname(func_cols),
      org_index = org_index
    ) %>%
      dplyr::mutate(
        transition = paste0(prefix_map[fun_name], "_", abbr),
        command    = sprintf(
          'Call["%s", FromTable["population_parameters.csv", %d, %d], %d]',
          fun_name, org_index, col_index, org_index
        )
      ) %>%
      dplyr::select(transition, command, fun_name, org_index, col_index)
  })

  # ------------------------------------------------------------
  # 3) Build arcs for population‐process transitions
  # ------------------------------------------------------------
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
        Starvation = tibble::tibble(
          transition, direction = "INPUT", place = biomass_place, multiplicity = 1L, command
        ),

        Duplication = dplyr::bind_rows(
          tibble::tibble(transition, direction = "INPUT",  place = count_place,   multiplicity = 1L, command),
          tibble::tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
          tibble::tibble(transition, direction = "OUTPUT", place = count_place,   multiplicity = 2L, command),
          tibble::tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 1L, command)
        ),

        Death = dplyr::bind_rows(
          tibble::tibble(transition, direction = "INPUT",  place = count_place,   multiplicity = 1L, command),
          tibble::tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
          tibble::tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 1L, command)
        )
      )
    }
  )

  # ------------------------------------------------------------
  # 4) Build arcs for FBA reactions
  # ------------------------------------------------------------
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
          dplyr::bind_rows(
            tibble::tibble(transition, direction = "INPUT",  place = met_place,     multiplicity = 1L, command),
            tibble::tibble(transition, direction = "INPUT",  place = count_place,   multiplicity = 1L, command),
            tibble::tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
            tibble::tibble(transition, direction = "OUTPUT", place = count_place,   multiplicity = 1L, command),
            tibble::tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 1L, command)
          )
        } else {
          # EXPORT
          dplyr::bind_rows(
            tibble::tibble(transition, direction = "INPUT",  place = count_place,   multiplicity = 1L, command),
            tibble::tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
            tibble::tibble(transition, direction = "OUTPUT", place = met_place,     multiplicity = 1L, command),
            tibble::tibble(transition, direction = "OUTPUT", place = count_place,   multiplicity = 1L, command),
            tibble::tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 1L, command)
          )
        }
      } else {
        # Biomass boundary
        if (endsWith(transition, paste0("_in_", abbr))) {
          tibble::tibble(transition, direction = "INPUT", place = biomass_place, multiplicity = 1L, command)
        } else {
          dplyr::bind_rows(
            tibble::tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
            tibble::tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 2L, command)
          )
        }
      }
    }
  )

  # ------------------------------------------------------------
  # 5) Combine repaired arcs
  # ------------------------------------------------------------
  arc_df_repaired <- dplyr::bind_rows(call_arcs, fba_arcs) %>%
    dplyr::distinct(transition, direction, place, multiplicity, command)

  # ------------------------------------------------------------
  # 6) Write CSV outputs
  # ------------------------------------------------------------
  readr::write_csv(arc_df,          file.path(out_dir, "raw_arcs.csv"))
  readr::write_csv(arc_df_repaired, file.path(out_dir, "repaired_arcs.csv"))

  # ------------------------------------------------------------
  # 7) Return invisibly
  # ------------------------------------------------------------
  invisible(list(raw_arcs = arc_df, repaired_arcs = arc_df_repaired))
}

