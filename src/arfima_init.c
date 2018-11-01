#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/
 
/* .C calls */
extern void durlevsim(void *, void *, void *, void *, void *, void *);
extern void shift_C(void *, void *, void *, void *);
extern void tacfFGN_C(void *, void *, void *);
extern void tacfHD_C(void *, void *, void *);
extern void tacvfARMA_C(void *, void *, void *, void *, void *, void *);
extern void tacvfFDWN_C(void *, void *, void *);
extern void tfcalc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(integd)(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"durlevsim",   (DL_FUNC) &durlevsim,    6},
    {"shift_C",     (DL_FUNC) &shift_C,      4},
    {"tacfFGN_C",   (DL_FUNC) &tacfFGN_C,    3},
    {"tacfHD_C",    (DL_FUNC) &tacfHD_C,     3},
    {"tacvfARMA_C", (DL_FUNC) &tacvfARMA_C,  6},
    {"tacvfFDWN_C", (DL_FUNC) &tacvfFDWN_C,  3},
    {"tfcalc",      (DL_FUNC) &tfcalc,      10},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"integd", (DL_FUNC) &F77_NAME(integd), 8},
    {NULL, NULL, 0}
};

void R_init_arfima(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
