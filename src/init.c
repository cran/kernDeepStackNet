#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>  // optional

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP kernDeepStackNet_crossprodRcpp(SEXP);
extern SEXP kernDeepStackNet_getEigenValuesRcpp(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"kernDeepStackNet_crossprodRcpp",      (DL_FUNC) &kernDeepStackNet_crossprodRcpp,      1},
    {"kernDeepStackNet_getEigenValuesRcpp", (DL_FUNC) &kernDeepStackNet_getEigenValuesRcpp, 1},
    {NULL, NULL, 0}
};

void R_init_kernDeepStackNet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
