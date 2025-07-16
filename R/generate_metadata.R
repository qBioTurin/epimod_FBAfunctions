#' Generate metadata CSV files for a folder of *.mat models (with debug)
#'
#' @param mat_dir    Character.  Folder that contains one or more `*.mat` files.
#' @param output_dir Destination root for the generated folders (default = mat_dir).
#' @param biomass    Named list with components `max`, `mean`, `min` (defaults to 1, 1, 0).
#'                   Passed to `setBiomassParameters()`.
#' @param overwrite  Logical.  Overwrite existing model folders (FALSE).
#' @param quiet      Logical.  Suppress progress messages (FALSE).
#' @param progress   Function(idx, total) for progress updates (NULL).
#' @return A tibble (invisibly) with one row per processed model, including the
#'         path where the metadata files were written.
#' @export
generate_metadata <- function(mat_dir,
                              output_dir = mat_dir,
                              biomass    = list(max = 1, mean = 1, min = 0),
                              overwrite  = FALSE,
                              quiet      = FALSE,
                              progress   = NULL) {
  # — normalize paths to avoid any “~/” confusion —
  mat_dir    <- normalizePath(mat_dir,    mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  stopifnot(dir.exists(mat_dir))

  mat_files <- list.files(mat_dir, pattern = "\\.mat$", full.names = TRUE)
  if (length(mat_files) == 0) {
    stop("No '.mat' files found in ", mat_dir)
  }

  say <- function(...) if (!quiet) message(...)

  res <- purrr::map_dfr(seq_along(mat_files), function(idx) {
    mat_path <- mat_files[idx]
    if (!is.null(progress)) progress(idx, length(mat_files))

    model_name <- tools::file_path_sans_ext(basename(mat_path))
    model_dir  <- fs::path(output_dir, model_name)

    if (fs::dir_exists(model_dir)) {
      if (overwrite) {
        fs::dir_delete(model_dir)
      } else {
        stop("Destination folder exists (set overwrite = TRUE): ", model_dir)
      }
    }
    fs::dir_create(model_dir, recurse = TRUE)

    # 1) copy original .mat file
    fs::file_copy(
      mat_path,
      fs::path(model_dir, paste0(model_name, ".mat")),
      overwrite = TRUE
    )

    # 2) conversion pipeline
    say("▶ [", model_name, "] FBA4Greatmod.generation()")
    model_obj <- FBA4Greatmod.generation(
      fba_mat   = fs::path(model_dir, paste0(model_name, ".mat")),
      input_dir = model_dir,
      GUI       = TRUE
    )
    model_obj <- setBiomassParameters(
      model_obj,
      bioMax  = biomass$max,
      bioMean = biomass$mean,
      bioMin  = biomass$min
    )

    # derive a short base name for the .txt
    abbrs   <- derive_abbrs(model_name)
    abbr    <- if (length(abbrs) >= 2) abbrs[2] else abbrs[1]
    fba_base <- paste0(abbr, "_model")

    say("▶ [", model_name, "] writeFBAfile()")
    # — DEBUG: show exactly what dest_dir and fname are —
    real_dest <- paste0(fs::path_norm(model_dir), .Platform$file.sep)
    message(
      "   → about to write FBA .txt:\n",
      "       dest_dir  = ", real_dest, "\n",
      "       fba_fname = ", fba_base, ".txt"
    )
    writeFBAfile(
      model_obj,
      fba_fname = fba_base,
      dest_dir  = real_dest
    )

    # 3) extract boundary metabolites
    meta_path   <- fs::path(model_dir, "metabolites_metadata.csv")
    boundary_ok <- FALSE
    if (fs::file_exists(meta_path)) {
      meta_df      <- readr::read_csv(meta_path, show_col_types = FALSE)
      boundary_df  <- dplyr::filter(meta_df, .data$is_boundary)
      readr::write_csv(boundary_df, fs::path(model_dir, "boundary_metabolites.csv"))
      boundary_ok <- nrow(boundary_df) > 0
    }

    meta_ok <- all(
      fs::file_exists(meta_path),
      fs::file_exists(fs::path(model_dir, "reactions_metadata.csv"))
    )

    tibble::tibble(
      model          = model_name,
      directory      = model_dir,
      txt_file       = fs::path(model_dir, paste0(fba_base, ".txt")),
      meta_ready     = meta_ok,
      boundary_ready = boundary_ok
    )
  })

  invisible(res)
}

