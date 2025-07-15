# R/model_analysis.R

#' General Wrapper for epimod::model.analysis (GUI)
#'
#' Provides a flexible interface to `epimod::model.analysis`, constructing
#' solver, parameters, function scripts, and auxiliary file paths from a
#' hypernode directory structure and user inputs.
#'
#' @param paths           Named character vector of hypernode paths:
#'                        "gen", "config", "src", "output".
#' @param hypernode_name  Name of the hypernode (used to build filenames).
#' @param debug_solver    Logical. Pass to `debug` argument in analysis.
#' @param i_time          Numeric. Initial time for simulation.
#' @param f_time          Numeric. Final time for simulation.
#' @param s_time          Numeric. Step size for simulation.
#' @param atol            Numeric. Absolute tolerance.
#' @param rtol            Numeric. Relative tolerance.
#' @param fba_fname       Character vector. FBA model `.txt` files (from biounits).
#' @param user_files      Character vector. Additional user data files.
#' @param ...             Additional args passed to `epimod::model.analysis`.
#'
#' @return Invisible list of results from `epimod::model.analysis`.
#' @export
model_analysis_GUI <- function(
  paths,
  hypernode_name,
  debug_solver = FALSE,
  i_time       = 0,
  f_time       = 10,
  s_time       = 1,
  atol         = 1e-6,
  rtol         = 1e-6,
  fba_fname,
  user_files = character(),
  volume     = getwd()
) {
  # Attach epimod if needed
  if (!"package:epimod" %in% search()) library(epimod)

  # Validate inputs
  stopifnot(
    is.character(paths), length(paths) >= 4,
    all(c("gen","config","src","output") %in% names(paths)),
    is.character(hypernode_name), nzchar(hypernode_name),
    is.logical(debug_solver), length(debug_solver)==1,
    is.numeric(i_time), is.numeric(f_time), is.numeric(s_time),
    is.numeric(atol), is.numeric(rtol),
    is.character(fba_fname), length(fba_fname)>=1,
    all(file.exists(fba_fname)),
    is.character(user_files)
  )

  # Build standard file paths
  solver_fname     <- fs::path(paths["gen"], paste0(hypernode_name, ".solver"))
  parameters_fname <- fs::path(paths["config"], "initial_data.csv")
  functions_fname  <- fs::path(paths["src"],   paste0("functions_", hypernode_name, ".R"))

  # Ensure files exist
  files_needed <- c(solver_fname, parameters_fname, functions_fname)
  if (!all(file.exists(files_needed))) {
    stop("Missing required analysis files: ",
      paste(basename(files_needed[!file.exists(files_needed)]), collapse = ", "))
  }

  # Call epimod analysis
  results <- epimod::model.analysis(
    solver_fname     = solver_fname,
    parameters_fname = parameters_fname,
    functions_fname  = functions_fname,
    debug            = debug_solver,
    i_time           = i_time,
    f_time           = f_time,
    s_time           = s_time,
    atol             = atol,
    rtol             = rtol,
    fba_fname        = fba_fname,
    user_files       = user_files,
    volume					 = volume
  )

  invisible(results)
}
