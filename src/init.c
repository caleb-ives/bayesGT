#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(sample)(double *p, int *Yt_mat, int *Z_mat, int *N,
                     int *Ycols, int *Zrows, int *Zcols, double *U,
                     int *Ztil, int *assay_vec, int *L,
                     double *se_in, double *sp_in,
                     int *se_counts, int *sp_counts);

static const R_FortranMethodDef fortranMethods[] = {
  {"sample", (DL_FUNC) &F77_NAME(sample), 15},
  {NULL, NULL, 0}
};

void R_init_bayesGT(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
