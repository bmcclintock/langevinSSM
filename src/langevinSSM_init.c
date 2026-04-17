#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton(getwd(),,,FALSE) (R 4.0)

/* .Call calls */
extern SEXP _langevinSSM_measurementError_LS_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _langevinSSM_measurementError_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _langevinSSM_simulate_langevin_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _langevinSSM_simulate_regionprob_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _langevinSSM_simulate_ud_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_langevinSSM_measurementError_LS_rcpp", (DL_FUNC) &_langevinSSM_measurementError_LS_rcpp, 7},
  {"_langevinSSM_measurementError_rcpp",    (DL_FUNC) &_langevinSSM_measurementError_rcpp,    6},
  {"_langevinSSM_simulate_langevin_cpp",    (DL_FUNC) &_langevinSSM_simulate_langevin_cpp,    9},
  {"_langevinSSM_simulate_regionprob_cpp",  (DL_FUNC) &_langevinSSM_simulate_regionprob_cpp,  8},
  {"_langevinSSM_simulate_ud_cpp",          (DL_FUNC) &_langevinSSM_simulate_ud_cpp,          7},
  {NULL, NULL, 0}
};

void R_init_langevinSSM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
