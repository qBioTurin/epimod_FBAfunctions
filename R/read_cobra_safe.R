#' Lettura sicura di file MAT (C++ nativo se presente, altrimenti fallback JSON)
#'
#' @param path Percorso al file `.mat`
#' @return Lista con campi `S`, `lb`, `ub`, `c`, `rxns`, `mets`
#' @export
read_cobra_safe <- function(path) {
  if (requireNamespace("RcppArmadillo", quietly = TRUE) &&
      "read_cobra_mat" %in% getNamespaceExports("epimodFBAfunctions")) {
    read_cobra_mat(path)
  } else {
    message("⚠️  Lettore C++ non disponibile: uso Python/JSON")
    read_model_json(sub("\\.mat$", ".json", path))
  }
}

