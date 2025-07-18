#' @keywords internal
#' @export
get_python_bin <- function() {
  # 1) opzione
  opt <- getOption("epimodFBAfunctions.python", "")
  if (nzchar(opt)) return(opt)

  # 2) variabile d'ambiente
  env <- Sys.getenv("EPIMOD_PYTHON", "")
  if (nzchar(env)) return(env)

  # 3) ricerca nel PATH
  py <- Sys.which("python3")
  if (nzchar(py)) return(py)

  stop(
    "Non trovo un interprete Python.\n",
    "• Installa python3 + cobra\n",
    "• oppure imposta options(epimodFBAfunctions.python='...')\n",
    "• oppure esporta EPIMOD_PYTHON='...'"
  )
}

