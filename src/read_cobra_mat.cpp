// ──────────────────────────────────────────────────────────────────────────────
//  read_cobra_mat.cpp   – lettore .mat nativo per modelli COBRA
//
//  Dipendenze: Rcpp, RcppArmadillo, libmatio ≥ 1.5
// ──────────────────────────────────────────────────────────────────────────────

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "matio.h"

using namespace Rcpp;

// ──────────────────────────────────────────────
// helper: prova più nomi (anche dentro struct 'model')
// ──────────────────────────────────────────────
inline matvar_t* pull(mat_t* mat,
                      const std::vector<std::string>& names) {

  /* 1) livello principale --------------------------------------------------- */
  for (const auto& n : names) {
    if (auto* v = Mat_VarRead(mat, n.c_str()))
      return v;                              // trovato
  }

  /* 2) dentro struct 'model' ------------------------------------------------ */
  if (matvar_t* mdl = Mat_VarRead(mat, "model")) {
    if (mdl->class_type == MAT_C_STRUCT) {
      for (const auto& n : names) {
        if (auto* v = Mat_VarGetStructFieldByName(mdl, n.c_str(), 0))
          return v;                          // (NON liberare qui)
      }
    }
    Mat_VarFree(mdl);                        // inutile → libera
  }

  stop("Variabile '%s' non trovata nel .mat",
       names.front().c_str());
  return nullptr;                            // mai eseguito
}

/* ──────────────────────────────────────────────
   Converter generico per:
     • char array   (MAT_C_CHAR)
     • cell array   (MAT_C_CELL) di char array
   ──────────────────────────────────────────── */
CharacterVector charvec(mat_t* mat,
                        const std::vector<std::string>& names) {

  matvar_t* v = pull(mat, names);
  CharacterVector out;

  if (v->class_type == MAT_C_CHAR) {         // char array  [nChar × nStrings]
    size_t nrow = v->dims[1], ncol = v->dims[0];
    const char* raw = static_cast<char*>(v->data);
    out = CharacterVector(nrow);
    for (size_t i = 0; i < nrow; ++i)
      out[i] = std::string(raw + i*ncol,
                           strnlen(raw + i*ncol, ncol));
  }
  else if (v->class_type == MAT_C_CELL) {    // cellstr 1×N
    size_t n = v->dims[1];
    out = CharacterVector(n);
    auto** cells = static_cast<matvar_t**>(v->data);

    for (size_t i = 0; i < n; ++i) {
      matvar_t* cell = cells[i];
      if (cell && cell->class_type == MAT_C_CHAR) {
        out[i] = std::string(static_cast<char*>(cell->data),
                             cell->nbytes / sizeof(char));
      } else {
        out[i] = NA_STRING;
      }
    }
  }
  else {
    Mat_VarFree(v);
    stop("Variabile %s: tipo non gestito", names.front().c_str());
  }

  Mat_VarFree(v);
  return out;
}

// [[Rcpp::export]]
List read_cobra_mat(const std::string& file) {

  /* apertura file ---------------------------------------------------------- */
  mat_t* mat = Mat_Open(file.c_str(), MAT_ACC_RDONLY);
  if (!mat) stop("Impossibile aprire %s", file.c_str());

  /* 1. matrice S ----------------------------------------------------------- */
  matvar_t* vS = pull(mat, {"S", "Stoichiometry"});
  if (vS->class_type != MAT_C_SPARSE)
    stop("'S' non è una matrice sparsa");

  mat_sparse_t* sp = static_cast<mat_sparse_t*>(vS->data);
  arma::sp_mat S( arma::uvec(sp->ir, sp->ndata, false, true),
                  arma::uvec(sp->jc, sp->njc, false, true),
                  arma::vec (static_cast<double*>(sp->data),
                             sp->ndata, false, true),
                  vS->dims[0], vS->dims[1] );
  Mat_VarFree(vS);

  /* 2. vettori numerici ---------------------------------------------------- */
  auto numvec = [&](const char* a, const char* b = nullptr) {
    matvar_t* v = pull(mat, {a, b ? std::string(b) : ""});
    NumericVector out(v->nbytes / sizeof(double));
    std::memcpy(out.begin(), v->data, v->nbytes);
    Mat_VarFree(v);
    return out;
  };
  NumericVector lb  = numvec("lb", "LBound");
  NumericVector ub  = numvec("ub", "UBound");
  NumericVector obj = numvec("c" , "obj");

  /* 3. vettori di stringhe -------------------------------------------------- */
  CharacterVector rxns       = charvec(mat, {"rxns", "reactionID"});
  CharacterVector mets       = charvec(mat, {"mets", "metaboliteID"});
  CharacterVector rxnNames   = charvec(mat, {"rxnNames", "reactionName"});
  CharacterVector metNames   = charvec(mat, {"metNames", "metaboliteName"});
  CharacterVector subSystems = charvec(mat, {"subSystems"});

  /* opzionali */
  CharacterVector rxnKEGG, metKEGG;
  try { rxnKEGG = charvec(mat, {"rxnKEGGID"}); }
  catch(...) { rxnKEGG = CharacterVector(rxns.size(), NA_STRING); }

  try { metKEGG = charvec(mat, {"metKEGGID"}); }
  catch(...) { metKEGG = CharacterVector(mets.size(), NA_STRING); }

  Mat_Close(mat);

  /* output ----------------------------------------------------------------- */
  return List::create(
    _["S"]          = S,
    _["lb"]         = lb,
    _["ub"]         = ub,
    _["c"]          = obj,
    _["rxns"]       = rxns,
    _["rxnNames"]   = rxnNames,
    _["rxnKEGGID"]  = rxnKEGG,
    _["mets"]       = mets,
    _["metNames"]   = metNames,
    _["metKEGGID"]  = metKEGG,
    _["subSystems"] = subSystems
  );
}

