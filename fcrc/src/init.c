#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _fcrc_admm_grplasso_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fcrc_admm_grplasso_path_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fcrc_alm_cgl_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _fcrc_alm_cgl_path_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_fcrc_admm_grplasso_int",      (DL_FUNC) &_fcrc_admm_grplasso_int,      11},
    {"_fcrc_admm_grplasso_path_int", (DL_FUNC) &_fcrc_admm_grplasso_path_int, 11},
    {"_fcrc_alm_cgl_int",            (DL_FUNC) &_fcrc_alm_cgl_int,            17},
    {"_fcrc_alm_cgl_path_int",       (DL_FUNC) &_fcrc_alm_cgl_path_int,       17},
    {NULL, NULL, 0}
};

void R_init_fcrc(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

