#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _openpopscr_C_calc_D(SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_llk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_move_llk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_pdet(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_move_pdet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_pr_capture(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_openpopscr_C_calc_D",          (DL_FUNC) &_openpopscr_C_calc_D,          4},
    {"_openpopscr_C_calc_llk",        (DL_FUNC) &_openpopscr_C_calc_llk,        8},
    {"_openpopscr_C_calc_move_llk",   (DL_FUNC) &_openpopscr_C_calc_move_llk,   12},
    {"_openpopscr_C_calc_pdet",       (DL_FUNC) &_openpopscr_C_calc_pdet,       5},
    {"_openpopscr_C_calc_move_pdet",  (DL_FUNC) &_openpopscr_C_calc_move_pdet,  10},
    {"_openpopscr_C_calc_pr_capture", (DL_FUNC) &_openpopscr_C_calc_pr_capture, 15},
    {NULL, NULL, 0}
};

void R_init_openpopscr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
