#include <R.h>
#include "Rmath.h"

extern "C" {

  void zhpevxC(char **jobz, char **range, char **uplo, int *n, Rcomplex *ap,
	      double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
	      double *w, Rcomplex *z, int *ldz, Rcomplex *work, double *rwork,
	      int *iwork, int *ifail, int *info) {

          extern void F77_NAME(zhpevx)(const void*, const void*, const void*, void*, void*, void*,
				       void*, void*, void*, void*, void*, void*,
				       void*, void*, void*, void*, void*, void*,
				       void*);

	  const char *JOBZ = jobz[0];
	  const char *RANGE = range[0];
	  const char *UPLO = uplo[0];
	  
	        F77_CALL(zhpevx)(JOBZ, RANGE, UPLO, &n[0], ap, &vl[0], &vu[0], &il[0], &iu[0], &abstol[0], &m[0],
				 w, z, &ldz[0], work, rwork, iwork, ifail, &info[0]);

  }

}
