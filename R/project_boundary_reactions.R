#' Project boundary reactions across a set of FBA models  — debug build
#'
#' (Same logic as your original; only extra `message()` lines.)
#' @export
#'
project_boundary_reactions <- function(biounit_models,
                                       boundary_metabolites,
                                       out_dir,
                                       hypernode_name,
                                       base_dir = getwd()) {

  message("\n─── building model dataframe ──────────────────────────────")

  # helper: drop any directory and the ".mat" extension
  clean_name <- function(x) tools::file_path_sans_ext(fs::path_file(x))
  models_df <- tibble::tibble(model = biounit_models) %>%
    dplyr::mutate(
      abbr = purrr::map_chr(model, ~ {
        raw <- .x$abbreviation[2]
        val <- sub("^/", "", clean_name(raw))      # remove leading slash
        message("   ↪ abbr extracted: ", val)
        val
      }),

      FBAname = purrr::map_chr(model, ~ clean_name(.x$FBAmodel)),

      model_file = fs::path(base_dir, "hypernodes", hypernode_name,
                           "biounits", FBAname, paste0(FBAname, ".mat")),
      meta_dir   = fs::path(base_dir, "hypernodes", hypernode_name,
                            "biounits", FBAname)
    )

  message("\n🔍 expecting these .mat files:")
  purrr::walk(models_df$model_file, ~ message("   ", .x))

  models_df <- dplyr::filter(models_df, file.exists(model_file))
  message("✔ models FOUND: ", nrow(models_df))

  models_df <- dplyr::mutate(
    models_df,
    metabolites = purrr::map(meta_dir, ~ {
      f <- file.path(.x, "metabolites_metadata.csv")
      message("   ↪ reading metabolites : ", f)
      readr::read_csv(f, show_col_types = FALSE)
    }),
    reactions = purrr::map(meta_dir, ~ {
      f <- file.path(.x, "reactions_metadata.csv")
      message("   ↪ reading reactions    : ", f)
      rx <- readr::read_csv(f, show_col_types = FALSE)
      message("     columns: ", paste(names(rx), collapse = ", "))
      message("     first 5 abbreviations: ",
              paste(head(rx$abbreviation, 5), collapse = ", "))
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
    tidyr::unnest(metabolites) %>%
    dplyr::filter(id %in% boundary_metabolites) %>%
    dplyr::distinct(abbr, id) %>%
    dplyr::left_join(boundary_df, by = "abbr",
                     relationship = "many-to-many") %>%
    dplyr::filter(stringr::str_detect(equation,
               stringr::str_c("\\b", id, "\\b"))) %>%
    dplyr::distinct(abbr, reaction, lowbnd, uppbnd)

  message("✔ shared_rxns_df rows: ", nrow(shared_rxns_df))

  # ──────────────────────────────────────────────────────────
  # Organise projections
  # ──────────────────────────────────────────────────────────
  shared_list <- shared_rxns_df %>%
    dplyr::group_by(abbr) %>%
    dplyr::summarise(rxns = list(reaction), .groups = "drop")

  all_models <- shared_list$abbr
  joint_rxns <- unique(shared_rxns_df$reaction)

  reaction_orgs <- shared_rxns_df %>%
    dplyr::distinct(abbr, reaction) %>%
    dplyr::group_by(reaction) %>%
    dplyr::summarise(orgs = list(abbr), .groups = "drop")

  common_rxns <- reaction_orgs %>%
    dplyr::filter(purrr::map_lgl(orgs, ~ setequal(.x, all_models))) %>%
    dplyr::pull(reaction)

  exclusive_list <- shared_list %>%
    dplyr::mutate(excl = purrr::map(rxns, ~ setdiff(.x, common_rxns)))

  message("✔ common reactions : ", length(common_rxns))
  message("✔ joint   reactions: ", length(joint_rxns))

  # ──────────────────────────────────────────────────────────
  # Write outputs
  # ──────────────────────────────────────────────────────────
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  readr::write_lines(joint_rxns,
                     file.path(out_dir, "extracted_boundary_reactions.txt"))

  if (length(common_rxns) > 0)
    readr::write_lines(common_rxns,
                       file.path(out_dir, "common_reactions.txt"))

  purrr::walk2(
    exclusive_list$excl, exclusive_list$abbr,
    ~ if (length(.x) > 0)
        readr::write_lines(.x,
          file.path(out_dir,
                    sprintf("exclusive_reactions_%s.txt", .y)))
  )

  # Build the bounds table
  bounds_tbl <- shared_rxns_df %>%
    dplyr::select(reaction, organism = abbr,
                  lower_bound = lowbnd, upper_bound = uppbnd) %>%
    dplyr::distinct()

  # always write a (possibly empty) reaction_bounds.csv
  readr::write_csv(bounds_tbl,
                   file.path(out_dir, "reaction_bounds.csv"))

  # and for each organism write per_bounds_<abbr>.csv (header only if empty)
  all_abbrs <- unique(models_df$abbr)
  purrr::walk(all_abbrs, function(ab) {
    out_file <- file.path(out_dir, sprintf("per_bounds_%s.csv", ab))
    df_sub <- bounds_tbl %>% dplyr::filter(organism == ab)
    readr::write_csv(df_sub, out_file)
  })

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
          forward_reaction = dplyr::if_else(is_reversible,
                                            paste0(abbreviation, "_f"),
                                            abbreviation),
          reverse_reaction = dplyr::if_else(is_reversible,
                                            paste0(abbreviation, "_r"),
                                            NA_character_),
          lb_f = dplyr::if_else(is_reversible, 0, pmax(0, lowbnd)),
          ub_f = dplyr::if_else(is_reversible, uppbnd, uppbnd),
          lb_r = dplyr::if_else(is_reversible, 0, NA_real_),
          ub_r = dplyr::if_else(is_reversible, -lowbnd, NA_real_)
        ) %>%
        dplyr::select(original_reaction = abbreviation,
                      forward_reaction, lb_f, ub_f,
                      reverse_reaction, lb_r, ub_r)

      irr_long <- irr_net %>%
        tidyr::pivot_longer(
          cols      = c(forward_reaction, reverse_reaction),
          names_to  = "direction",
          values_to = "reaction") %>%
        dplyr::filter(!is.na(reaction)) %>%
        dplyr::mutate(
          lb = dplyr::if_else(direction == "forward_reaction", lb_f, lb_r),
          ub = dplyr::if_else(direction == "forward_reaction", ub_f, ub_r)
        ) %>%
        dplyr::select(reaction, lb, ub)

      irr_out <- file.path(out_dir,
                           sprintf("irreversible_reactions_%s.csv",
                                   unique(.y$abbr)))
      readr::write_csv(irr_long, irr_out)
      cat(sprintf("  ✓ irreversible saved for %s\n",
                  unique(.y$abbr)))
    })

  invisible(list(
    shared_reactions    = purrr::set_names(shared_list$rxns,
                                           shared_list$abbr),
    common_reactions    = common_rxns,
    exclusive_reactions = purrr::set_names(exclusive_list$excl,
                                           exclusive_list$abbr),
    bounds              = bounds_tbl
  ))
}

