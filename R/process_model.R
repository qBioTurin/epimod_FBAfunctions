#' Process a single biounit model
#'
#' * Copies the chosen *.mat* file into the hyper-node workspace  
#' * Converts it to a text FBA model with `writeFBAfile()`  
#' * Returns the path to that single *.txt* file
#'
#' @param m              A list for one cellular unit (from
#'                       `make_biounit_models()`).
#' @param hypernode_name Character. Parent hyper-node label.
#' @param mat_dir        Directory that contains the user-supplied
#'                       *.mat* files (same value you pass to
#'                       `build_hypernode()`).
#' @param base_dir       Root working directory (default = `getwd()`).
#'
#' @return (Invisibly) a list with status and created paths.
#' @export
process_model <- function(m,
                          hypernode_name,
                          mat_dir,
                          base_dir = getwd()) {

  cat("\n========== PROCESSING BIOUNIT ==========\n")

  ## -----------------------------------------------------------------
  ##  helpers & metadata
  ## -----------------------------------------------------------------
  abs_path <- function(...) fs::path_abs(fs::path(...))

  unit     <- m$unit
  abbr     <- m$abbreviation[2]
  FBAmodel <- ifelse(is.null(m$FBAmodel), unit, m$FBAmodel)

  ## -----------------------------------------------------------------
  ##  locate the .mat file
  ## -----------------------------------------------------------------
  if (fs::file_exists(FBAmodel)) {
    ## caller supplied a full path
    mat_path <- fs::path_abs(FBAmodel)
    FBAmodel <- tools::file_path_sans_ext(fs::path_file(mat_path))
  } else {
    ## plain name → look inside mat_dir/
    mat_candidate <- fs::path(mat_dir, paste0(FBAmodel, ".mat"))
    if (!fs::file_exists(mat_candidate))
      stop("MAT model not found: ", mat_candidate)
    mat_path <- fs::path_abs(mat_candidate)
  }

  ## -----------------------------------------------------------------
  ##  workspace within the hyper-node
  ## -----------------------------------------------------------------
  input_dir <- abs_path(base_dir, "hypernodes", hypernode_name,
                        "biounits", FBAmodel)
  fs::dir_create(input_dir, recurse = TRUE)

  ## -----------------------------------------------------------------
  ##  copy *.mat*  →  generate *.txt*
  ## -----------------------------------------------------------------
  fs::file_copy(mat_path,
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

  ## -----------------------------------------------------------------
  ##  write the FBA text file  (exactly one “.txt”)
  ## -----------------------------------------------------------------
  fba_base     <- paste0(abbr, "_model")                 # no extension
  output_file  <- fs::path(input_dir, paste0(fba_base, ".txt"))

  cat("📝 Writing model with writeFBAfile() …\n")

  writeFBAfile(
    model_obj,
    fba_fname = fba_base,        # writeFBAfile() appends “.txt”
    dest_dir  = input_dir
  )

  ## -----------------------------------------------------------------
  ##  confirm + return
  ## -----------------------------------------------------------------
  if (fs::file_exists(output_file)) {
    cat("✅ Model saved to:", output_file, "\n")
    invisible(list(
      status     = "success",
      unit       = unit,
      abbr       = abbr,
      model_file = output_file
    ))
  } else {
    cat("❌ writeFBAfile() did not create", output_file, "\n")
    invisible(list(
      status  = "error",
      message = "Model file not created",
      unit    = unit,
      abbr    = abbr
    ))
  }
}

