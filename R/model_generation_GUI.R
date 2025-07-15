#' General Wrapper for epimod::model.generation (GUI)
#'
#' Ensures 'epimod' is attached so that underlying calls to path.package()
#' within epimod::model.generation() succeed.
#'
#' @param net_fname         Path to your .PNPRO file
#' @param transitions_fname Path to your .cpp transitions file
#' @param fba_fname         One or more .txt FBA model files
#' @param output_dir        Where to write generated outputs
#' @param ...               Passed to epimod::model.generation()
#' @return Invisibly, the result list from epimod::model.generation()
#' @export
model_generation_GUI <- function(
  net_fname,
  transitions_fname,
  fba_fname,
  output_dir = getwd(),
  ...
) {
  # make sure 'epimod' is attached so path.package() won't fail
  if (!"package:epimod" %in% search()) {
    library(epimod)
  }

  # validate inputs
  stopifnot(
    is.character(net_fname), length(net_fname) == 1, file.exists(net_fname),
    is.character(transitions_fname), length(transitions_fname) == 1, file.exists(transitions_fname),
    is.character(fba_fname), length(fba_fname) >= 1, all(file.exists(fba_fname)),
    is.character(output_dir), length(output_dir) == 1
  )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  results <- epimod::model.generation(
    net_fname         = net_fname,
    transitions_fname = transitions_fname,
    fba_fname         = fba_fname,
    volume = output_dir
  )

  invisible(results)
}

