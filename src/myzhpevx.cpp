#include <Rconfig.h>       // for FC_LEN_T
#include <R_ext/Complex.h> // for Rcomplex
#include <R_ext/RS.h>      // for F77_CALL, F77_NAME



extern "C" {
void F77_NAME(zhpevx)(char *jobz, char *range, char *uplo,
                      int *n, Rcomplex *ap, double *vl, double *vu,
                      int *il, int *iu, double *abstol, int *m, double *w,
                      Rcomplex *z, int *ldz, Rcomplex *work, double *rwork,
                      int *iwork, int *ifail, int *info
#ifdef FC_LEN_T
                      , FC_LEN_T, FC_LEN_T, FC_LEN_T
#endif
);
}

extern "C" {
void zhpevxC(char **jobz, char **range, char **uplo, int *n, Rcomplex *ap,
             double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
             double *w, Rcomplex *z, int *ldz, Rcomplex *work, double *rwork,
             int *iwork, int *ifail, int *info)
{
    // Pass single characters explicitly
  char jobz_char = jobz[0][0];
  char range_char = range[0][0];
  char uplo_char = uplo[0][0];


    F77_CALL(zhpevx)(&jobz_char, &range_char, &uplo_char,
                     n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work,
                     rwork, iwork, ifail, info
#ifdef FC_LEN_T
                     , (FC_LEN_T)1, (FC_LEN_T)1, (FC_LEN_T)1
#endif
    );
}
}
