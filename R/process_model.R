#' Process a single biounit model
#'
#' @param m             list returned by make_biounit_models()
#' @param hypernode_name character, parent hyper-node label
#' @param base_dir      where hypernodes/, metabolic_networks_library/ …
#'                      live (default = getwd()).
#' @return list with status + paths (invisibly)
#' @export
process_model <- function(m,
                          hypernode_name,
                          base_dir = getwd()) {

  cat("\n========== PROCESSING BIOUNIT ==========\n")

  ## -----------------------------------------------------------------
  ##  local helpers & paths
  ## -----------------------------------------------------------------
  abs_path <- function(...) fs::path_abs(fs::path(...))

  unit      <- m$unit
  abbr      <- m$abbreviation[2]
  FBAmodel  <- ifelse(is.null(m$FBAmodel), unit, m$FBAmodel)

  # where the original .mat sits
  mat_file  <- abs_path(base_dir, "metabolic_networks_library",
                        paste0(FBAmodel, ".mat"))

  # where model-specific artefacts go
  input_dir <- abs_path(base_dir, "hypernodes", hypernode_name,
                        "biounits", FBAmodel)

  output_file <- fs::path(input_dir, paste0(abbr, "_model.txt"))
  fs::dir_create(input_dir, recurse = TRUE)

  ## -----------------------------------------------------------------
  ##  copy *.mat + generate .txt with writeFBAfile()
  ## -----------------------------------------------------------------
  fs::file_copy(mat_file,
                fs::path(input_dir, paste0(FBAmodel, ".mat")),
                overwrite = TRUE)

  model_obj <- FBA4Greatmod.generation(
    fba_mat   = fs::path(input_dir, paste0(FBAmodel, ".mat")),
    input_dir = input_dir
  )

  model_obj <- setBiomassParameters(
    model_obj,
    bioMax  = m$biomass$max,
    bioMean = m$biomass$mean,
    bioMin  = m$biomass$min
  )

  cat("📝 Writing model with writeFBAfile...\n")
  before_files <- fs::dir_ls(fs::path_dir(input_dir),
                             glob = "*.txt")

  writeFBAfile(model_obj,
               fba_fname = paste0(abbr, "_model.txt"),
               dest_dir  = input_dir)

  after_files  <- fs::dir_ls(fs::path_dir(input_dir),
                             glob = "*.txt")
  new_txt      <- setdiff(after_files, before_files)

  if (length(new_txt) == 1) {
    fs::file_move(new_txt, output_file)
    cat("✅ Model saved to:", output_file, "\n")
    invisible(list(
      status      = "success",
      unit        = unit,
      abbr        = abbr,
      model_file  = output_file
    ))
  } else {
    cat("❌ Could not locate new .txt model output.\n")
    invisible(list(
      status      = "error",
      message     = "Model file not created",
      unit        = unit,
      abbr        = abbr
    ))
  }
}

