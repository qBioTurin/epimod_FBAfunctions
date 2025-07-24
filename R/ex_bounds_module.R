# ============================================================
#  File: ex_bounds_module_fixed.R      (2-CSV workflow)
# ============================================================

# --------------------------------------------------------------------
# Utility: strip directory + ".mat", leaving the clean FBA-model name
# --------------------------------------------------------------------
#' Strip the directory and “.mat” extension from a path.
#' @keywords internal
#' @export
model_dir <- function(path) {
  tools::file_path_sans_ext(fs::path_file(path))
}

# --------------------------------------------------------------------
# Helper: split one reaction into the directions actually needed
# --------------------------------------------------------------------
split_reaction <- function(row) {
  lb   <- row$lowbnd
  ub   <- row$uppbnd
  base <- row$abbreviation

  out <- list()

  ## 1) reversible → both dirs
  if (lb < 0 && ub > 0) {
    out$f <- data.frame(
      reaction     = base,
      abbreviation = paste0(base, "_f"),
      name         = paste0(row$name, " (forward)"),
      equation     = paste("EXPORT:", row$equation),
      subtype      = row$subtype,
      FBAmodel     = row$FBAmodel,
      txt_file     = row$txt_file,
      upper_bound  = ub,
      stringsAsFactors = FALSE
    )
    out$r <- data.frame(
      reaction     = base,
      abbreviation = paste0(base, "_r"),
      name         = paste0(row$name, " (reverse)"),
      equation     = paste("IMPORT:", row$equation),
      subtype      = row$subtype,
      FBAmodel     = row$FBAmodel,
      txt_file     = row$txt_file,
      upper_bound  = abs(lb),
      stringsAsFactors = FALSE
    )

  ## 2) irreversible reverse
  } else if (lb < 0 && ub <= 0) {
    out$r <- data.frame(
      reaction     = base,
      abbreviation = paste0(base, "_r"),
      name         = paste0(row$name, " (reverse)"),
      equation     = paste("IMPORT:", row$equation),
      subtype      = row$subtype,
      FBAmodel     = row$FBAmodel,
      txt_file     = row$txt_file,
      upper_bound  = abs(lb),
      stringsAsFactors = FALSE
    )

  ## 3) irreversible forward
  } else if (lb >= 0 && ub > 0) {
    out$f <- data.frame(
      reaction     = base,
      abbreviation = paste0(base, "_f"),
      name         = paste0(row$name, " (forward)"),
      equation     = paste("EXPORT:", row$equation),
      subtype      = row$subtype,
      FBAmodel     = row$FBAmodel,
      txt_file     = row$txt_file,
      upper_bound  = ub,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(out)
}

# --------------------------------------------------------------------
# Expand and write bound adjustments for boundary reactions
# --------------------------------------------------------------------
process_boundary_reactions <- function(
    hypernode_name,
    biounit_models,
    output_dir,
    projected_lower_bound,
    projected_upper_bound,
    not_projected_lower_bound,
    not_projected_upper_bound,
    volume,                    # unico per tutti
    projected_reactions_df,
    base_dir = getwd()
) {

  ## 1) metadati ----------------------------------------------------------
  metadata_all <- do.call(rbind, lapply(biounit_models, function(model) {
    model_name <- model_dir(model$FBAmodel)
    meta_path  <- fs::path(base_dir, "hypernodes", hypernode_name,
                           "biounits", model_name, "reactions_metadata.csv")
    if (!file.exists(meta_path))
      stop("Missing reactions_metadata.csv at: ", meta_path)

    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    meta$FBAmodel <- model_name
    meta$txt_file <- model$txt_file
    meta
  }))
  metadata_all <- subset(metadata_all, type == "boundary")

  ## 2) espansione dir ----------------------------------------------------
  expanded <- purrr::map_dfr(seq_len(nrow(metadata_all)),
                             ~ split_reaction(metadata_all[.x, ]))

  ## 3) flag projected ----------------------------------------------------
  expanded$projected <- mapply(function(rxn, org) {
    any(projected_reactions_df$reaction  == rxn &
        projected_reactions_df$FBAmodel == org)
  }, expanded$abbreviation, expanded$FBAmodel)

  ## 4) build due data-frame  --------------------------------------------
  df_proj  <- tibble::tibble()
  df_nproj <- tibble::tibble()

  for (i in seq_len(nrow(expanded))) {
    row        <- expanded[i, ]
    is_reverse <- endsWith(row$abbreviation, "_r")

    ## ---- PROJECTED ------------------------------------------------------
    if (row$projected) {
      bound <- if (is_reverse) abs(projected_lower_bound) else projected_upper_bound

      df_proj <- dplyr::bind_rows(df_proj, data.frame(
        reaction        = row$abbreviation,
        FBAmodel        = row$FBAmodel,
        background_conc = bound,
        system_volume   = volume,
        stringsAsFactors = FALSE
      ))

    ## ---- NON-PROJECTED --------------------------------------------------
    } else {
      bound <- if (is_reverse) abs(not_projected_lower_bound) else not_projected_upper_bound

      df_nproj <- dplyr::bind_rows(df_nproj, data.frame(
        reaction        = row$abbreviation,
        FBAmodel        = row$FBAmodel,
        background_conc = bound,
        system_volume   = volume,
        stringsAsFactors = FALSE
      ))
    }
  }

  ## 5) write the two CSV --------------------------------------------------
  utils::write.table(
    df_proj,
    file.path(output_dir, "ub_bounds_projected.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  utils::write.table(
    df_nproj,
    file.path(output_dir, "non_projected_bounds.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  message(
    "✔ Bounds generated:\n",
    "  • ub_bounds_projected.csv\n",
    "  • non_projected_bounds.csv"
  )
}

# --------------------------------------------------------------------
# Master wrapper
# --------------------------------------------------------------------
#' @export
run_full_ex_bounds <- function(
    hypernode_name,
    biounit_models,
    projected_lower_bound,
    projected_upper_bound,
    not_projected_lower_bound,
    not_projected_upper_bound,
    volume,
    base_dir = getwd()
) {
  output_dir <- fs::path(base_dir, "hypernodes", hypernode_name, "output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ## -- reaction_bounds.csv ------------------------------------------------
  reaction_bounds_path <- fs::path(output_dir, "reaction_bounds.csv")
  if (!file.exists(reaction_bounds_path))
    stop("Missing: ", reaction_bounds_path)

  reaction_bounds_df <- read.csv(reaction_bounds_path, stringsAsFactors = FALSE)

  ## -- abbreviazione → nome modello completo ------------------------------
  abbr_map <- setNames(
    vapply(biounit_models, function(x) model_dir(x$FBAmodel), character(1)),
    vapply(biounit_models, function(x) x$abbreviation[2],     character(1))
  )

  ## -- tabella split per proiezioni ---------------------------------------
  projected_reactions_df <- purrr::map_dfr(seq_len(nrow(reaction_bounds_df)), function(i) {
    row        <- reaction_bounds_df[i, ]
    full_model <- abbr_map[[row$organism]]
    if (is.null(full_model)) return(NULL)

    dummy <- data.frame(
      abbreviation = row$reaction,
      name         = "",
      equation     = "",
      subtype      = "",
      lowbnd       = row$lower_bound,
      uppbnd       = row$upper_bound,
      FBAmodel     = full_model,
      txt_file     = NA,
      stringsAsFactors = FALSE
    )

    split_reaction(dummy) %>%
      dplyr::transmute(
        reaction = abbreviation,
        FBAmodel = FBAmodel
      )
  })

  ## -- dispatch -----------------------------------------------------------
  process_boundary_reactions(
    hypernode_name             = hypernode_name,
    biounit_models             = biounit_models,
    output_dir                 = output_dir,
    projected_lower_bound      = projected_lower_bound,
    projected_upper_bound      = projected_upper_bound,
    not_projected_lower_bound  = not_projected_lower_bound,
    not_projected_upper_bound  = not_projected_upper_bound,
    volume                     = volume,
    projected_reactions_df     = projected_reactions_df,
    base_dir                   = base_dir
  )
}

