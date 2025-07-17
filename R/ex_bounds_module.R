# ============================================================
#  File: ex_bounds_module_fixed.R  (path-handling fixes applied)
# ============================================================

# strip directory + “.mat”, leaving the clean FBA-model folder name
model_dir <- function(path) {
  tools::file_path_sans_ext(fs::path_file(path))
}

# ------------------------------------------------------------
# Expand and write bound adjustments for boundary reactions
# Added base_dir argument and absolute paths for metadata
process_boundary_reactions <- function(
    hypernode_name,
    biounit_models,
    output_dir,
    bacteria_counts,
    projected_base_lb,
    projected_base_ub,
    not_projected_base_lb,
    not_projected_base_ub,
    projected_reactions_df,
    base_dir = getwd()
) {
  ## ── 1) Load and expand boundary reactions ───────────────────────
  metadata_all <- do.call(rbind, lapply(biounit_models, function(model) {
    model_name <- model_dir(model$FBAmodel)
    # build absolute path to reactions metadata
    meta_path <- fs::path(
      base_dir,
      "hypernodes", hypernode_name,
      "biounits", model_name,
      "reactions_metadata.csv"
    )
    if (!file.exists(meta_path)) {
      stop("Missing reactions_metadata.csv at: ", meta_path)
    }
    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    meta$FBAmodel <- model_name
    meta$txt_file <- model$txt_file
    meta
  }))
  metadata_all <- subset(metadata_all, type == "boundary")

  expanded <- do.call(rbind, lapply(seq_len(nrow(metadata_all)), function(i) {
    row      <- metadata_all[i, , drop = FALSE]
    base_rxn <- row$abbreviation
    data.frame(
      reaction     = rep(base_rxn, 2),
      abbreviation = c(paste0(base_rxn, "_r"), paste0(base_rxn, "_f")),
      name         = c(paste0(row$name, " (reverse)"),
                       paste0(row$name, " (forward)")),
      equation     = c(paste("IMPORT:", row$equation),
                       paste("EXPORT:", row$equation)),
      subtype      = rep(row$subtype, 2),
      FBAmodel     = rep(row$FBAmodel, 2),
      txt_file     = rep(row$txt_file, 2),
      upper_bound  = c(abs(row$lowbnd), row$uppbnd),
      stringsAsFactors = FALSE
    )
  }))

  expanded$projected <- mapply(function(rxn, org) {
    any(projected_reactions_df$reaction  == rxn &
        projected_reactions_df$FBAmodel == org)
  }, expanded$abbreviation, expanded$FBAmodel)

  cat("Total reactions in expanded:", nrow(expanded), "\n")
  cat("Marked as projected:", sum(expanded$projected), "\n")
  print(utils::head(expanded[expanded$projected, ], 5))

  ## ── 2) Compute bounds ───────────────────────────────────────────
  projected     <- list()
  non_projected <- list()

  for (i in seq_len(nrow(expanded))) {
    row <- expanded[i, ]
    if (grepl("EX_biomass_e", row$abbreviation)) next

    rxn_type <- if (grepl("_f$", row$abbreviation)) "f" else "r"

    if (row$projected) {
      ub <- if (rxn_type == "r") abs(projected_base_lb) else projected_base_ub
      projected[[length(projected) + 1]] <-
        data.frame(
          reaction    = row$abbreviation,
          FBAmodel    = row$FBAmodel,
          upper_bound = ub,
          stringsAsFactors = FALSE
        )
    } else {
      if (rxn_type == "r") {
        if (row$subtype == "exchange") {
          match_rows  <- expanded[expanded$abbreviation == row$abbreviation &
                                  !expanded$projected, ]
          orgs        <- unique(match_rows$FBAmodel)
          total_count <- sum(bacteria_counts[orgs])
          ub          <- abs(not_projected_base_lb / total_count)
        } else {
          cell_cnt <- bacteria_counts[row$FBAmodel]
          ub       <- abs(not_projected_base_lb / cell_cnt)
        }
      } else {
        ub <- not_projected_base_ub
      }

      non_projected[[length(non_projected) + 1]] <-
        data.frame(
          reaction    = row$abbreviation,
          FBAmodel    = row$FBAmodel,
          upper_bound = ub,
          stringsAsFactors = FALSE
        )
    }
  }

  # ensure empty data.frames have the correct columns
  df_proj <- if (length(projected)) {
    do.call(rbind, projected)
  } else {
    data.frame(
      reaction    = character(),
      FBAmodel    = character(),
      upper_bound = double(),
      stringsAsFactors = FALSE
    )
  }

  df_nonproj <- if (length(non_projected)) {
    do.call(rbind, non_projected)
  } else {
    data.frame(
      reaction    = character(),
      FBAmodel    = character(),
      upper_bound = double(),
      stringsAsFactors = FALSE
    )
  }

  # write CSVs with headers
  utils::write.table(
    df_proj,
    file.path(output_dir, "ub_bounds_projected.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  utils::write.table(
    df_nonproj,
    file.path(output_dir, "ub_bounds_not_projected.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  message("✔ Bounds generated:\n - ub_bounds_projected.csv\n - ub_bounds_not_projected.csv")
}

# ------------------------------------------------------------------
# Master function
# ------------------------------------------------------------------
#' @export
run_full_ex_bounds <- function(
    hypernode_name,
    biounit_models,
    projected_base_lb,
    projected_base_ub,
    not_projected_base_lb,
    not_projected_base_ub,
    base_dir = getwd()
) {
  # build absolute output directory
  output_dir <- fs::path(base_dir, "hypernodes", hypernode_name, "output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # counts keyed by clean FBA-model name
  bacteria_counts <- setNames(
    vapply(biounit_models, function(x) x$initial_count, numeric(1)),
    vapply(biounit_models, function(x) model_dir(x$FBAmodel), character(1))
  )

  # Load projection info ------------------------------------------------
  reaction_bounds_path <- fs::path(output_dir, "reaction_bounds.csv")
  if (!file.exists(reaction_bounds_path))
    stop("Missing: ", reaction_bounds_path)

  reaction_bounds_df <- read.csv(reaction_bounds_path, stringsAsFactors = FALSE)

  # Abbreviation ➜ model map
  abbr_map <- setNames(
    vapply(biounit_models, function(x) model_dir(x$FBAmodel), character(1)),
    vapply(biounit_models, function(x) x$abbreviation[2],     character(1))
  )

  # Build _r/_f projection table
  projected_reactions_df <- do.call(rbind, lapply(seq_len(nrow(reaction_bounds_df)), function(i) {
    row <- reaction_bounds_df[i, ]
    full_model <- abbr_map[[row$organism]]
    if (is.null(full_model)) {
      warning("No FBAmodel found for abbreviation: ", row$organism)
      return(NULL)
    }
    data.frame(
      reaction = c(paste0(row$reaction, "_r"), paste0(row$reaction, "_f")),
      type     = c("r", "f"),
      FBAmodel = rep(full_model, 2),
      bound    = c(row$lower_bound, row$upper_bound),
      stringsAsFactors = FALSE
    )
  }))

  # Dispatch ------------------------------------------------------------
  process_boundary_reactions(
    hypernode_name          = hypernode_name,
    biounit_models          = biounit_models,
    output_dir              = output_dir,
    bacteria_counts         = bacteria_counts,
    projected_base_lb       = projected_base_lb,
    projected_base_ub       = projected_base_ub,
    not_projected_base_lb   = not_projected_base_lb,
    not_projected_base_ub   = not_projected_base_ub,
    projected_reactions_df  = projected_reactions_df,
    base_dir                = base_dir
  )
}

