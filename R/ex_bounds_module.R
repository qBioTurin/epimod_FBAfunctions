# ============================================================
#  File: ex_bounds_module_fixed.R  (unificato in 2 CSV)
# ============================================================

# … funzioni ausiliarie model_dir()  e split_reaction() INVARIATE …


# --------------------------------------------------------------------
# Crea **due soli** file:
#   • ub_bounds_projected.csv
#   • non_projected_bounds.csv
# --------------------------------------------------------------------
process_boundary_reactions <- function(
    hypernode_name,
    biounit_models,
    output_dir,
    projected_lower_bound,
    projected_upper_bound,
    not_projected_lower_bound,
    not_projected_upper_bound,
    volume,                                  # unico per tutti
    projected_reactions_df,
    base_dir = getwd()
) {

  # --- 1) Metadati delle reazioni di confine ------------------------------
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

  # --- 1b) Espansione in forward/reverse ----------------------------------
  expanded <- purrr::map_dfr(seq_len(nrow(metadata_all)),
                             ~ split_reaction(metadata_all[.x, ]))

  # --- 1c) Flag “projected / non-projected” -------------------------------
  expanded$projected <- mapply(function(rxn, org) {
    any(projected_reactions_df$reaction  == rxn &
        projected_reactions_df$FBAmodel == org)
  }, expanded$abbreviation, expanded$FBAmodel)

  # --- 2) Calcolo dei nuovi bound + build tabelle -------------------------
  df_proj <- data.frame()        # bounds per metaboliti PROJECTED
  df_nproj<- data.frame()        # bounds per metaboliti NON-PROJECTED

  for (i in seq_len(nrow(expanded))) {
    row <- expanded[i, ]
    if (grepl("EX_biomass_e", row$abbreviation, fixed = TRUE)) next

    is_reverse <- endsWith(row$abbreviation, "_r")

    ## PROJECTED ------------------------------------------------------------
    if (row$projected) {

      bound <- if (is_reverse) abs(projected_lower_bound) else projected_upper_bound

      df_proj <- rbind(df_proj, data.frame(
        reaction        = row$abbreviation,
        FBAmodel        = row$FBAmodel,
        background_conc = bound,
        system_volume   = volume,
        stringsAsFactors = FALSE
      ))

    ## NON-PROJECTED --------------------------------------------------------
    } else {

      bound <- if (is_reverse) abs(not_projected_lower_bound) else not_projected_upper_bound

      df_nproj <- rbind(df_nproj, data.frame(
        reaction        = row$abbreviation,
        FBAmodel        = row$FBAmodel,
        background_conc = bound,
        system_volume   = volume,
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- 3) Scrittura file CSV ----------------------------------------------
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
    " - ub_bounds_projected.csv\n",
    " - non_projected_bounds.csv"
  )
}

# --------------------------------------------------------------------
# Master wrapper  (aggiornato con i nuovi argomenti)
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

  # … lettura reaction_bounds.csv e costruzione projected_reactions_df  INVARIATA …

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

