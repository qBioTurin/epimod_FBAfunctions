# ---------------------------------------------------------------------------
# Generate a shortlist of plausible abbreviations for one model name ─ DEBUG
# ---------------------------------------------------------------------------
#' @export
derive_abbrs <- function(model_name) {

  ## ─── debug ────────────────────────────────────────────────────────
  message("derive_abbrs(): raw model_name = ", model_name)

  # drop any trailing “_model”
  core  <- stringr::str_remove(model_name, "_model$")
  parts <- strsplit(core, "_")[[1]]
  abbrs <- character()

  # candidate 1: everything before the first underscore
  abbrs <- c(abbrs, parts[1])

  # candidate 2: first letters of each part
  abbrs <- c(abbrs, paste0(substr(parts, 1, 1), collapse = ""))

  # candidate 3: if there are exactly two parts, 1st letter + 2nd letter
  if (length(parts) == 2)
    abbrs <- c(abbrs,
               paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1)))

  abbrs <- unique(tolower(abbrs))

  ## ─── debug ────────────────────────────────────────────────────────
  message("derive_abbrs(): candidates -> ",
          paste(abbrs, collapse = ", "))

  abbrs
}

# ---------------------------------------------------------------------------
# Build a general-purpose list of FBA-compatible biounits ─ DEBUG + FIX
# ---------------------------------------------------------------------------
#' @export
make_biounit_models <- function(model_names,
                                biomass_params,
                                population_params,
                                initial_counts) {

  stopifnot(length(model_names) == length(biomass_params),
            length(model_names) == length(population_params),
            length(model_names) == length(initial_counts))

  lapply(seq_along(model_names), function(i) {

    raw_path <- model_names[i]

    ## ─── NEW: strip directory & “.mat” ──────────────────────────────
    short    <- tools::file_path_sans_ext(fs::path_file(raw_path))
    # e.g. “/…/Escherichia_coli_SE11.mat” → “Escherichia_coli_SE11”

    abbrs <- derive_abbrs(short)

    ## ─── debug ──────────────────────────────────────────────────────
    message(sprintf("make_biounit_models(): [%02d] %s", i, short))
    message("  › abbr[1] = ", abbrs[1],
            "; abbr[2] = ", abbrs[2])

    list(
      # keep the full path so downstream functions can still read the file
      FBAmodel            = raw_path,
      unit                = gsub("_", " ", short),
      abbreviation        = abbrs,
      txt_file            = paste0(abbrs[2], "_model.txt"),
      biomass             = biomass_params[[i]],
      population_settings = population_params[[i]],
      initial_count       = initial_counts[i]
    )
  })
}

# ---------------------------------------------------------------------------
# Write population dynamics parameters (unchanged, just with a tiny ping)
# ---------------------------------------------------------------------------
#' @export
write_population_params <- function(biounit_models, path) {

  message("write_population_params(): writing → ", path)

  # Turn each population_settings list into a 1 × 3 data.frame
  df <- do.call(rbind, lapply(biounit_models, function(unit) {
    as.data.frame(as.list(unit$population_settings), stringsAsFactors = FALSE)
  }))

  # Save as plain CSV (no headers, no row names)
  write.table(df, path,
              row.names = FALSE,
              col.names = FALSE,
              sep       = ",",
              quote     = FALSE)
}

