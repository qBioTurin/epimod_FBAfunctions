#' Project boundary reactions across a set of FBA models  — debug build
#'
#' (Same logic as your original; only extra `message()` lines.)
#' @export
project_boundary_reactions <- function(biounit_models,
                                       boundary_metabolites,
                                       out_dir,
                                       hypernode_name,
                                       base_dir = getwd()) {

  # ──────────────────────────────────────────────────────────
  # Build dataframe describing each model
  # ──────────────────────────────────────────────────────────
  models_df <- tibble(model = biounit_models) %>%
    mutate(
      abbr     = map_chr(model, ~ .x$abbreviation[2]),
      FBAmodel = map_chr(model, ~ .x$FBAmodel),
      model_file = file.path("hypernodes", hypernode_name,
                             "biounits", FBAmodel, paste0(FBAmodel, ".mat")),
      meta_dir   = file.path("hypernodes", hypernode_name,
                             "biounits", FBAmodel)
    )

  message("🔍 expecting these .mat files:\n",
          paste(models_df$model_file, collapse = "\n  "))

  models_df <- models_df %>%
    filter(file.exists(model_file))

  message("✔ models found: ", nrow(models_df))

  models_df <- models_df %>%
    mutate(
      metabolites = map(meta_dir, ~ {
        f <- file.path(.x, "metabolites_metadata.csv")
        message("   ↪ reading metabolites: ", f)
        read_csv(f, show_col_types = FALSE)
      }),
      reactions   = map(meta_dir, ~ {
        f <- file.path(.x, "reactions_metadata.csv")
        message("   ↪ reading reactions  : ", f)
        rx <- read_csv(f, show_col_types = FALSE)
        message("     columns: ", paste(names(rx), collapse = ", "))
        rx
      })
    )

  # ──────────────────────────────────────────────────────────
  # Extract boundary reactions
  # ──────────────────────────────────────────────────────────
  boundary_df <- models_df %>%
    dplyr::select(abbr, reactions) %>%
    tidyr::unnest(reactions) %>%
    dplyr::filter(type == "boundary") %>%
    dplyr::select(abbr, reaction = abbreviation,
                  lowbnd, uppbnd, equation)

  message("✔ boundary_df rows: ", nrow(boundary_df))

  # Match projectable metabolites to boundary reactions
  shared_rxns_df <- models_df %>%
    unnest(metabolites) %>%
    filter(id %in% boundary_metabolites) %>%
    distinct(abbr, id) %>%
    left_join(boundary_df, by = "abbr",
              relationship = "many-to-many") %>%
    filter(str_detect(equation, str_c("\\b", id, "\\b"))) %>%
    distinct(abbr, reaction, lowbnd, uppbnd)

  message("✔ shared_rxns_df rows: ", nrow(shared_rxns_df))

  # ──────────────────────────────────────────────────────────
  # Organise projections
  # ──────────────────────────────────────────────────────────
  shared_list <- shared_rxns_df %>%
    group_by(abbr) %>%
    summarize(rxns = list(reaction), .groups = "drop")

  all_models <- shared_list$abbr
  joint_rxns <- unique(shared_rxns_df$reaction)

  reaction_orgs <- shared_rxns_df %>%
    distinct(abbr, reaction) %>%
    group_by(reaction) %>%
    summarize(orgs = list(abbr), .groups = "drop")

  common_rxns <- reaction_orgs %>%
    filter(map_lgl(orgs, ~ setequal(.x, all_models))) %>%
    pull(reaction)

  exclusive_list <- shared_list %>%
    mutate(excl = map(rxns, ~ setdiff(.x, common_rxns)))

  message("✔ common reactions : ", length(common_rxns))
  message("✔ joint   reactions: ", length(joint_rxns))

  # ──────────────────────────────────────────────────────────
  # Write outputs
  # ──────────────────────────────────────────────────────────
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  write_lines(joint_rxns,
              file.path(out_dir, "extracted_boundary_reactions.txt"))

  if (length(common_rxns) > 0)
    write_lines(common_rxns,
                file.path(out_dir, "common_reactions.txt"))

  walk2(
    exclusive_list$excl, exclusive_list$abbr,
    ~ if (length(.x) > 0)
        write_lines(.x,
          file.path(out_dir,
                    sprintf("exclusive_reactions_%s.txt", .y)))
  )

  bounds_tbl <- shared_rxns_df %>%
    dplyr::select(reaction, organism = abbr,
                  lower_bound = lowbnd, upper_bound = uppbnd) %>%
    dplyr::distinct()

  write_csv(bounds_tbl,
            file.path(out_dir, "reaction_bounds.csv"))

  bounds_tbl %>%
    group_by(organism) %>%
    group_walk(~ write_csv(
      .x,
      file.path(out_dir, sprintf("per_bounds_%s.csv", .y$organism))))

  # =============================
  # Generate irreversible network version
  # =============================
  cat("\nGenerating irreversible network references...\n")

  models_df %>%
    dplyr::select(abbr, reactions) %>%
    unnest(reactions) %>%
    group_by(abbr) %>%
    group_walk(~ {
      irr_net <- .x %>%
        mutate(
          is_reversible = lowbnd < 0 & uppbnd > 0,
          forward_reaction = if_else(is_reversible,
                                     paste0(abbreviation, "_f"),
                                     abbreviation),
          reverse_reaction = if_else(is_reversible,
                                     paste0(abbreviation, "_r"),
                                     NA_character_),
          lb_f = if_else(is_reversible, 0, pmax(0, lowbnd)),
          ub_f = if_else(is_reversible, uppbnd, uppbnd),
          lb_r = if_else(is_reversible, 0, NA_real_),
          ub_r = if_else(is_reversible, -lowbnd, NA_real_)
        ) %>%
        dplyr::select(original_reaction = abbreviation,
                      forward_reaction, lb_f, ub_f,
                      reverse_reaction, lb_r, ub_r)

      irr_long <- irr_net %>%
        pivot_longer(cols = c(forward_reaction, reverse_reaction),
                     names_to = "direction",
                     values_to = "reaction") %>%
        filter(!is.na(reaction)) %>%
        mutate(
          lb = if_else(direction == "forward_reaction", lb_f, lb_r),
          ub = if_else(direction == "forward_reaction", ub_f, ub_r)
        ) %>%
        dplyr::select(reaction, lb, ub)

      irr_out <- file.path(out_dir,
                           sprintf("irreversible_reactions_%s.csv",
                                   unique(.y$abbr)))
      write_csv(irr_long, irr_out)
      cat(sprintf("  ✓ irreversible saved for %s: %s\n",
                  unique(.y$abbr), irr_out))
    })

  invisible(list(
    shared_reactions    = set_names(shared_list$rxns, shared_list$abbr),
    common_reactions    = common_rxns,
    exclusive_reactions = set_names(exclusive_list$excl,
                                    exclusive_list$abbr),
    bounds              = bounds_tbl
  ))
}

