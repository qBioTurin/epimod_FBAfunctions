// src/read_cobra_mat.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
#include <Rcpp.h>
#include <matio.h>
using namespace Rcpp;

// build a dgCMatrix from a MATLAB CSC sparse matvar
static S4 matSparseToDgC(matvar_t* v) {
  if (v->isComplex || v->class_type != MAT_C_SPARSE)
    stop("Expected a real sparse matrix for 'S'.");

  mat_sparse_t* s = (mat_sparse_t*)v->data;
  size_t nz = s->nzmax, nrow = v->dims[0], ncol = v->dims[1];

  IntegerVector p(ncol+1), i(nz);
  NumericVector x(nz);

  std::copy(s->jc, s->jc + ncol+1, p.begin());
  std::copy(s->ir, s->ir + nz,      i.begin());
  std::copy((double*)s->data, (double*)s->data + nz, x.begin());

  S4 M("dgCMatrix");
  M.slot("i")   = i;
  M.slot("p")   = p;
  M.slot("x")   = x;
  M.slot("Dim") = IntegerVector::create(nrow, ncol);
  M.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
  return M;
}

// read a MATLAB cell of uint8 strings into CharacterVector
static CharacterVector getStrings(matvar_t* v) {
  if (!v || v->class_type != MAT_C_CELL) return CharacterVector();
  size_t n = 1;
  for (int d=0; d<v->rank; ++d) n *= v->dims[d];
  CharacterVector out(n);
  matvar_t** c = (matvar_t**)v->data;
  for (size_t k=0; k<n; ++k) {
    if (c[k] && c[k]->data_type == MAT_T_UINT8) {
      // build std::string of exact length then let Rcpp convert
      out[k] = std::string((char*)c[k]->data, c[k]->nbytes);
    } else {
      out[k] = NA_STRING;
    }
  }
  return out;
}

// read a MATLAB double array into NumericVector
static NumericVector getNumeric(matvar_t* v) {
  if (!v || v->class_type != MAT_C_DOUBLE) return NumericVector();
  size_t n = 1;
  for (int d=0; d<v->rank; ++d) n *= v->dims[d];
  return NumericVector((double*)v->data, (double*)v->data + n);
}

// [[Rcpp::export]]
List read_cobra_mat(std::string path) {
  mat_t* m = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
  if (!m) stop("Cannot open MAT file '%s'.", path);

  // get directory of variables (libmatio v1.5 signature)
  size_t nv;
  char** dir = Mat_GetDir(m, &nv);
  if (nv < 1) {
    free(dir);
    Mat_Close(m);
    stop("No variables in MAT file.");
  }

  // first var is the COBRA struct
  matvar_t* model = Mat_VarRead(m, dir[0]);
  free(dir);
  if (!model || model->class_type != MAT_C_STRUCT) {
    Mat_VarFree(model);
    Mat_Close(m);
    stop("First variable is not a struct; is this a COBRA model?");
  }

  auto fld = [&](const char* n){ 
    return Mat_VarGetStructFieldByName(model, n, 0); 
  };

  // mandatory fields
  matvar_t *vS    = fld("S"),
           *vlb   = fld("lb"),
           *vub   = fld("ub"),
           *vC    = fld("c"),
           *vrxns = fld("rxns"),
           *vmets = fld("mets");
  if (!vS || !vlb || !vub || !vC || !vrxns || !vmets) {
    Mat_VarFree(model);
    Mat_Close(m);
    stop("Missing one of S / lb / ub / c / rxns / mets in the model.");
  }

  // convert to R objects
  S4            S        = matSparseToDgC(vS);
  NumericVector lb       = getNumeric(vlb),
                ub       = getNumeric(vub),
                obj      = getNumeric(vC);
  CharacterVector rxns   = getStrings(vrxns),
                  mets   = getStrings(vmets);

  // optional
  CharacterVector rxnNames   = getStrings(fld("rxnNames")),
                  metNames   = getStrings(fld("metNames")),
                  subSystems = getStrings(fld("subSystems")),
                  genes      = getStrings(fld("genes")),
                  grRules    = getStrings(fld("grRules")),
                  rxnKEGGID  = getStrings(fld("rxnKEGGID")),
                  metKEGGID  = getStrings(fld("metKEGGID"));
  NumericVector rev        = getNumeric(fld("rev")),
                rxnGeneMat = getNumeric(fld("rxnGeneMat"));

  Mat_VarFree(model);
  Mat_Close(m);

  return List::create(
    _["S"]            = S,
    _["lb"]           = lb,
    _["ub"]           = ub,
    _["c"]            = obj,
    _["rxns"]         = rxns,
    _["rxnNames"]     = rxnNames,
    _["mets"]         = mets,
    _["metNames"]     = metNames,
    _["subSystems"]   = subSystems,
    _["rev"]          = rev,
    _["rxnGeneMat"]   = rxnGeneMat,
    _["genes"]        = genes,
    _["grRules"]      = grRules,
    _["rxnKEGGID"]    = rxnKEGGID,
    _["metKEGGID"]    = metKEGGID
  );
}

