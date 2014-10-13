#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

void massdist(double *x, double *xmass, int *nx, double *xlow, double *xhigh, double *y, int *ny);
void nncg(int *Rn, double *x0, double *A, double *b, double *fmin, int *fail, int *Rtype, int *Rtrace, int *objfcount, int *gradcount);

static const R_CMethodDef cMethods[] = {
  {"massdist", (DL_FUNC) &massdist, 7},
  {"nncg", (DL_FUNC) &nncg, 10},
  NULL
};

void R_init_SSPA(DllInfo *info)
{
  R_registerRoutines(info, cMethods,  NULL, NULL, NULL);
};
