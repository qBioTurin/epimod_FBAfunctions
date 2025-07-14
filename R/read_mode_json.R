# file: R/read_model_json.R

#' Read a COBRA model JSON into the same list structure as read_cobra_mat()
#'
#' @param json_path Path to the JSON file produced by mat2json.py
#' @return A named list with the *old* MAT field names:
#'   S, lb, ub, c, rxns, rxnNames, mets, metNames, subSystems, rev,
#'   rxnGeneMat, genes, grRules, rxnKEGGID, metKEGGID
#' @importFrom Matrix sparseMatrix
#' @importFrom jsonlite fromJSON
#' @export
read_model_json <- function(json_path) {
  j <- jsonlite::fromJSON(json_path)
  S <- Matrix::sparseMatrix(
    i        = unlist(j$S_i) + 1L,
    j        = unlist(j$S_j) + 1L,
    x        = unlist(j$S_x),
    dims     = j$dims,
    dimnames = list(j$mets, j$rxns)
  )
  # **Use the old slot names here** so FBAmat.read() keeps working
  list(
    S           = S,
    lb          = j$lb,
    ub          = j$ub,
    c           = j$c,
    rxns        = j$rxns,
    rxnNames    = j$rxnNames,
    mets        = j$mets,
    metNames    = j$metNames,
    subSystems  = j$subSystems,
    rev         = j$rev,
    rxnGeneMat  = j$rxnGeneMat,
    genes       = j$genes,
    grRules     = j$grRules,
    rxnKEGGID   = j$rxnKEGGID,
    metKEGGID   = j$metKEGGID
  )
}

