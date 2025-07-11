#' Install developer tools and optional containers
#'
#' *End-users never need this; it is for contributors setting up a fresh
#' machine.  Runtime dependencies are already installed automatically
#' because they are declared in DESCRIPTION.*
#'
#' @param with_containers Logical.  If `TRUE` run
#'   \code{downloadContainers()} after packages are installed.
#' @export
install_dev_env <- function(with_containers = TRUE) {

  # developer-only helpers ----------------------------------------------
  dev_pkgs <- c("devtools", "pkgload", "usethis")

  for (p in dev_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }

  ## optional: refresh devtools
  if (packageVersion("devtools") < "2.4.0") {
    message("Upgrading devtools …")
    utils::remove.packages("devtools")
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }

  # install epimod from GitHub if absent -------------------------------
  if (!requireNamespace("epimod", quietly = TRUE)) {
    devtools::install_github("qBioTurin/epimod", ref = "epimod_pFBA")
  }

  # download containers -------------------------------------------------
  if (isTRUE(with_containers)) {
    if (!exists("downloadContainers", mode = "function")) {
      stop("downloadContainers() not found; please source or import it.")
    }
    downloadContainers()
  }

  invisible(TRUE)
}

