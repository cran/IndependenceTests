#include <stddef.h>

#include <R_ext/Complex.h>
#include <R_ext/RS.h>

#if defined(__clang__)
# include <complex>
typedef std::complex<double> fcomplex;
#else
typedef __complex__ double fcomplex;
#endif

extern "C" {

void F77_NAME(zhpevx)(char *jobz, char *range, char *uplo,
                      int *n, fcomplex *ap,
                      double *vl, double *vu,
                      int *il, int *iu, double *abstol,
                      int *m, double *w,
                      fcomplex *z, int *ldz,
                      fcomplex *work,
                      double *rwork, int *iwork, int *ifail,
                      int *info,
                      size_t jobz_len, size_t range_len, size_t uplo_len);

void zhpevxC(char **jobz, char **range, char **uplo, int *n, Rcomplex *ap,
             double *vl, double *vu, int *il, int *iu, double *abstol,
             int *m, double *w, Rcomplex *z, int *ldz, Rcomplex *work,
             double *rwork, int *iwork, int *ifail, int *info)
{
    char jobz_char = jobz[0][0];
    char range_char = range[0][0];
    char uplo_char = uplo[0][0];

    F77_CALL(zhpevx)(&jobz_char, &range_char, &uplo_char,
                     n,
                     reinterpret_cast<fcomplex *>(ap),
                     vl, vu, il, iu, abstol, m, w,
                     reinterpret_cast<fcomplex *>(z),
                     ldz,
                     reinterpret_cast<fcomplex *>(work),
                     rwork, iwork, ifail, info,
                     (size_t) 1, (size_t) 1, (size_t) 1);
}

}

