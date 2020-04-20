#include <R.h>
#include "Rmath.h"
#include <complex.h>

extern"C" {

  void phinhatC(double *vect, double *X, int *q, int *n, double _Complex *res) {

    int j, col;
    double tmp = 0.0;
    double _Complex somme = 0.0 + 0.0*_Complex_I;
    
    for (j = 1; j<= n[0]; j++) {
      tmp = 0.0;
      for (col = 1; col <= q[0]; col++) {
	tmp = tmp + vect[col-1] * X[(col-1) * n[0] + (j-1)];
      }
      somme = somme + cexp(_Complex_I * tmp);
    }
    res[0] = somme / (double)n[0];
  }

}



