#' Project boundary reactions across a set of FBA models
#'
#' @export
project_boundary_reactions <- function(biounit_models,
                                       boundary_metabolites,
                                       out_dir,
                                       hypernode_name) {
  # Build dataframe describing each model
  models_df <- tibble::tibble(model = biounit_models) %>%
    dplyr::mutate(
      abbr     = purrr::map_chr(model, ~ .x$abbreviation[2]),
      FBAmodel = purrr::map_chr(model, ~ .x$FBAmodel),
      model_file = file.path("hypernodes", hypernode_name, "biounits", FBAmodel, paste0(FBAmodel, ".mat")),
      meta_dir   = file.path("hypernodes", hypernode_name, "biounits", FBAmodel)
    ) %>%
    dplyr::filter(file.exists(model_file)) %>%
    dplyr::mutate(
      metabolites = purrr::map(meta_dir, ~ readr::read_csv(file.path(.x, "metabolites_metadata.csv"), show_col_types = FALSE)),
      reactions   = purrr::map(meta_dir, ~ readr::read_csv(file.path(.x, "reactions_metadata.csv"), show_col_types = FALSE))
    )

  boundary_df <- models_df %>%
    dplyr::select(abbr, reactions) %>%
    tidyr::unnest(reactions) %>%
    dplyr::filter(type == "boundary") %>%
    dplyr::select(abbr, reaction = abbreviation, lowbnd, uppbnd, equation)

  # Match projectable metabolites to boundary reactions
  shared_rxns_df <- models_df %>%
    tidyr::unnest(metabolites) %>%
    dplyr::filter(id %in% boundary_metabolites) %>%
    dplyr::distinct(abbr, id) %>%
    dplyr::left_join(boundary_df, by = "abbr", relationship = "many-to-many") %>%
    dplyr::filter(stringr::str_detect(equation, stringr::str_c("\\b", id, "\\b"))) %>%
    dplyr::distinct(abbr, reaction, lowbnd, uppbnd)

  # Organize projections
  shared_list <- shared_rxns_df %>%
    dplyr::group_by(abbr) %>%
    dplyr::summarize(rxns = list(reaction), .groups = "drop")

  all_models <- shared_list$abbr
  joint_rxns <- unique(shared_rxns_df$reaction)

  reaction_orgs <- shared_rxns_df %>%
    dplyr::distinct(abbr, reaction) %>%
    dplyr::group_by(reaction) %>%
    dplyr::summarize(orgs = list(abbr), .groups = "drop")

  common_rxns <- reaction_orgs %>%
    dplyr::filter(purrr::map_lgl(orgs, ~ setequal(.x, all_models))) %>%
    dplyr::pull(reaction)

  exclusive_list <- shared_list %>%
    dplyr::mutate(excl = purrr::map(rxns, ~ setdiff(.x, common_rxns)))

  # Output
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  readr::write_lines(joint_rxns, file.path(out_dir, "extracted_boundary_reactions.txt"))

  if (length(common_rxns) > 0) {
    readr::write_lines(common_rxns, file.path(out_dir, "common_reactions.txt"))
  }

  purrr::walk2(
    exclusive_list$excl, exclusive_list$abbr,
    ~ if (length(.x) > 0) {
        readr::write_lines(.x, file.path(out_dir, sprintf("exclusive_reactions_%s.txt", .y)))
      }
  )

  bounds_tbl <- shared_rxns_df %>%
    dplyr::select(reaction, organism = abbr, lower_bound = lowbnd, upper_bound = uppbnd) %>%
    dplyr::distinct()

  readr::write_csv(bounds_tbl, file.path(out_dir, "reaction_bounds.csv"))

  bounds_tbl %>%
    dplyr::group_by(organism) %>%
    dplyr::group_walk(~ readr::write_csv(.x, file.path(out_dir, sprintf("per_bounds_%s.csv", .y$organism))))

  # =============================
  # Generate irreversible network version
  # =============================
  cat("\nGenerating irreversible network references...\n")

  models_df %>%
    dplyr::select(abbr, reactions) %>%
    tidyr::unnest(reactions) %>%
    dplyr::group_by(abbr) %>%
    dplyr::group_walk(~ {
      irr_net <- .x %>%
        dplyr::mutate(
          is_reversible = lowbnd < 0 & uppbnd > 0,
          forward_reaction = dplyr::if_else(is_reversible, paste0(abbreviation, "_f"), abbreviation),
          reverse_reaction = dplyr::if_else(is_reversible, paste0(abbreviation, "_r"), NA_character_),
          lb_f = dplyr::if_else(is_reversible, 0, pmax(0, lowbnd)),
          ub_f = dplyr::if_else(is_reversible, uppbnd, uppbnd),
          lb_r = dplyr::if_else(is_reversible, 0, NA_real_),
          ub_r = dplyr::if_else(is_reversible, -lowbnd, NA_real_)
        ) %>%
        dplyr::select(
          original_reaction = abbreviation,
          forward_reaction, lb_f, ub_f,
          reverse_reaction, lb_r, ub_r
        )

      irr_long <- irr_net %>%
        tidyr::pivot_longer(
          cols = c(forward_reaction, reverse_reaction),
          names_to = "direction",
          values_to = "reaction"
        ) %>%
        dplyr::filter(!is.na(reaction)) %>%
        dplyr::mutate(
          lb = dplyr::if_else(direction == "forward_reaction", lb_f, lb_r),
          ub = dplyr::if_else(direction == "forward_reaction", ub_f, ub_r)
        ) %>%
        dplyr::select(reaction, lb, ub)

      irr_output_path <- file.path(out_dir, sprintf("irreversible_reactions_%s.csv", unique(.y$abbr)))
      readr::write_csv(irr_long, irr_output_path)

      cat(sprintf("  \u2713 Irreversible network saved for %s: %s\n", unique(.y$abbr), irr_output_path))
    })

  invisible(list(
    shared_reactions    = purrr::set_names(shared_list$rxns, shared_list$abbr),
    common_reactions    = common_rxns,
    exclusive_reactions = purrr::set_names(exclusive_list$excl, exclusive_list$abbr),
    bounds              = bounds_tbl
  ))
}

