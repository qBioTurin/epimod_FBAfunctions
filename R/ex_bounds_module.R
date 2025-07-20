# ============================================================
#  File: ex_bounds_module_fixed.R  (path-handling fixes applied)
# ============================================================

# --------------------------------------------------------------------
# Utility: strip directory + ".mat", leaving the clean FBA‑model folder
# --------------------------------------------------------------------
model_dir <- function(path) {
  tools::file_path_sans_ext(fs::path_file(path))
}

# --------------------------------------------------------------------
# Helper: split one reaction into the directions that are actually needed
# --------------------------------------------------------------------
split_reaction <- function(row) {
  lb   <- row$lowbnd
  ub   <- row$uppbnd
  base <- row$abbreviation

  out <- list()

  # CASE 1 — reversible (lb < 0 & ub > 0) → _r + _f
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
  # CASE 2 — irreversible reverse (lb < 0 & ub ≤ 0) → _r only
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
  # CASE 3 — irreversible forward (lb ≥ 0 & ub > 0) → _f only
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
  # CASE 4 — blocked (lb == ub == 0) → no variables

  dplyr::bind_rows(out)
}

# --------------------------------------------------------------------
# Expand and write bound adjustments for boundary reactions
# --------------------------------------------------------------------
process_boundary_reactions <- function(
    hypernode_name,
    biounit_models,
    output_dir,
    bacteria_counts,
    projected_base_lb,
    projected_base_ub,
    background_met_base,
    volume,
    not_projected_base_ub,
    projected_reactions_df,
    base_dir = getwd()
) {
  # 1) Load boundary‑reaction metadata -------------------------------
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

  # 1b) Expand each row according to canonical splitting rules --------
  expanded <- purrr::map_dfr(seq_len(nrow(metadata_all)),
                             ~ split_reaction(metadata_all[.x, ]))

  # 1c) Mark the reactions that are projected ------------------------
  expanded$projected <- mapply(function(rxn, org) {
    any(projected_reactions_df$reaction  == rxn &
        projected_reactions_df$FBAmodel == org)
  }, expanded$abbreviation, expanded$FBAmodel)

  cat("Total expanded boundary reactions:", nrow(expanded), "\n")
  cat("  ↳ projected:", sum(expanded$projected), "\n")

  # 2) Compute new upper bounds -------------------------------------
  projected     <- list()
  nonproj_f     <- list()
  nonproj_r     <- list()

  for (i in seq_len(nrow(expanded))) {
    row <- expanded[i, ]
    if (grepl("EX_biomass_e", row$abbreviation, fixed = TRUE)) next

    rxn_type <- if (endsWith(row$abbreviation, "_f")) "f" else "r"

    if (row$projected) {
      ub <- if (rxn_type == "r") abs(projected_base_lb) else projected_base_ub
      projected[[length(projected) + 1]] <-
        data.frame(reaction     = row$abbreviation,
                   FBAmodel     = row$FBAmodel,
                   upper_bound  = ub,
                   stringsAsFactors = FALSE)
    } else {
      # Non-projected: forward uses base ub, reverse uses background_met_base
      if (rxn_type == "r") {
        ub <- abs(background_met_base)
        nonproj_r[[length(nonproj_r) + 1]] <-
          data.frame(reaction       = row$abbreviation,
                     FBAmodel       = row$FBAmodel,
                     background_met = ub,
                     volume         = volume,
                     stringsAsFactors = FALSE)
      } else {
        ub <- not_projected_base_ub
        nonproj_f[[length(nonproj_f) + 1]] <-
          data.frame(reaction       = row$abbreviation,
                     FBAmodel       = row$FBAmodel,
                     background_met = ub,
                     stringsAsFactors = FALSE)
      }
    }
  }

  df_proj      <- if (length(projected))     dplyr::bind_rows(projected)     else data.frame(reaction=character(),FBAmodel=character(),upper_bound=double())
  df_nonproj_f <- if (length(nonproj_f))     dplyr::bind_rows(nonproj_f)    else data.frame(reaction=character(),FBAmodel=character(),background_met=double())
  df_nonproj_r <- if (length(nonproj_r))     dplyr::bind_rows(nonproj_r)    else data.frame(reaction=character(),FBAmodel=character(),background_met=double(),volume=double())

  # 3) Write projected bounds (unchanged)
  utils::write.table(df_proj,
    file.path(output_dir, "ub_bounds_projected.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

  # 4) Write non-projected background_met split files
  utils::write.table(df_nonproj_f,
    file.path(output_dir, "non_projected_reverse_bounds.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  utils::write.table(df_nonproj_r,
    file.path(output_dir, "non_projected_reverse_background_met.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

  message("✔ Bounds generated:\n - ub_bounds_projected.csv\n",
          " - background_met_not_projected_f.csv\n",
          " - background_met_not_projected_r.csv")
}

# --------------------------------------------------------------------
# Master wrapper
# --------------------------------------------------------------------
#' @export
run_full_ex_bounds <- function(
    hypernode_name,
    biounit_models,
    projected_base_lb,
    projected_base_ub,
    background_met_base,
    volume,
    not_projected_base_ub,
    base_dir = getwd()
) {
  output_dir <- fs::path(base_dir, "hypernodes", hypernode_name, "output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # counts keyed by clean FBA model name -----------------------------
  bacteria_counts <- setNames(
    vapply(biounit_models, function(x) x$initial_count, numeric(1)),
    vapply(biounit_models, function(x) model_dir(x$FBAmodel), character(1))
  )

  # Load projection info --------------------------------------------
  reaction_bounds_path <- fs::path(output_dir, "reaction_bounds.csv")
  if (!file.exists(reaction_bounds_path))
    stop("Missing: ", reaction_bounds_path)

  reaction_bounds_df <- read.csv(reaction_bounds_path, stringsAsFactors = FALSE)

  # Abbreviation ➜ model map
  abbr_map <- setNames(
    vapply(biounit_models, function(x) model_dir(x$FBAmodel), character(1)),
    vapply(biounit_models, function(x) x$abbreviation[2],     character(1))
  )

  # Build projection table with the **same** splitting logic --------
  projected_reactions_df <- purrr::map_dfr(seq_len(nrow(reaction_bounds_df)), function(i) {
    row <- reaction_bounds_df[i, ]
    full_model <- abbr_map[[row$organism]]
    if (is.null(full_model)) {
      warning("No FBAmodel found for abbreviation: ", row$organism)
      return(NULL)
    }
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
    split_reaction(dummy) |>  # returns only required directions
      dplyr::transmute(reaction = abbreviation,
                       type     = ifelse(endsWith(reaction, "_f"), "f", "r"),
                       FBAmodel = FBAmodel,
                       bound    = upper_bound)
  })

  # Dispatch ---------------------------------------------------------
  process_boundary_reactions(
    hypernode_name         = hypernode_name,
    biounit_models         = biounit_models,
    output_dir             = output_dir,
    bacteria_counts        = bacteria_counts,
    projected_base_lb      = projected_base_lb,
    projected_base_ub      = projected_base_ub,
    background_met_base    = background_met_base,
    volume                 = volume,
    not_projected_base_ub  = not_projected_base_ub,
    projected_reactions_df = projected_reactions_df,
    base_dir               = base_dir
  )
}

