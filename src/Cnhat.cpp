#include <R.h>
#include "Rmath.h"
#include <complex.h>
#include <iostream>
using namespace std;

extern"C" {

  void CnhatC(double *vecs, double *vect, double *X, int *n, int *q, int *p, int *vecd, double _Complex *res) {

    void phinhatReturn(double *vect1, double *vect2, double *vect3, double *X, int *q, int *n, double _Complex *res1, double _Complex *res2, double _Complex *res3);
    double _Complex tmp2, prod1 = 1.0 + 0.0*_Complex_I, prod2 = 1.0 + 0.0*_Complex_I, somme = 0.0 + 0.0*_Complex_I;
    int indbeginblocl, l, i, *vecdl, k, j;
    double *vecsl, *mvectl, *diffvecsvectl, *Xl;
    double _Complex *res1, *res2, *res3;
    vecdl = new int[1];
    res1 =  new _Complex double[1];
    res2 =  new _Complex double[1];
    res3 =  new _Complex double[1];

    for (l = 0; l <= (p[0]-1); l++) {
      indbeginblocl = 1;
      if (l != 0) {
	for (i = 0; i <= (l-1); i++) {
	  indbeginblocl = indbeginblocl + vecd[i];
	}
      }
      vecdl[0] = vecd[l];
      vecsl = new double[vecdl[0]];
      for (k = 0; k <= (vecdl[0]-1); k++) vecsl[k] = vecs[indbeginblocl+k-1];
      mvectl = new double[vecdl[0]];
      diffvecsvectl = new double[vecdl[0]];
      for (k = 0; k <= (vecdl[0]-1); k++) {
	mvectl[k] = - vect[indbeginblocl+k-1];
	diffvecsvectl[k] = vecsl[k] + mvectl[k];
      }
      Xl = new double[n[0]*vecdl[0]];
      for (j = 0; j <= (n[0]-1); j++) {
	for (k = 0; k <= (vecdl[0]-1); k++) {
	  Xl[k*n[0] + j] = X[(indbeginblocl+k-1)*n[0] + j];
	}
      }
      phinhatReturn(diffvecsvectl,vecsl,mvectl,Xl,vecdl,n,res1,res2,res3);

      prod1 = prod1 * res1[0];
      tmp2 = res2[0] * res3[0];
      prod2 = prod2 * tmp2;
      somme = somme + res1[0] / tmp2;

      delete[] vecsl;
      delete[] mvectl;
      delete[] diffvecsvectl;
      delete[] Xl;
    }
    delete[] vecdl;
    delete[] res1;
    delete[] res2;
    delete[] res3;
    res[0] = prod1 - prod2 * (1 - p[0] + somme);

  }

  void CnhatmatC(double *yMat, int *N, double *X, int *n, int *q, int *p, int *vecd, double _Complex *res) {

    int i, j, k, qm = q[0] - 1, Nm = N[0] - 1;
    double *vect, *vecs;
    double _Complex *restmp;
    restmp = new _Complex double[1];
    vecs = new double[q[0]];
    vect = new double[q[0]];

    for (i = 0; i <= Nm; i++) {
      for (j = 0; j <= i; j++) {
	for (k = 0; k <= qm; k++) {
	  vecs[k] = yMat[k*N[0] + i];
	  vect[k] = yMat[k*N[0] + j];
	}
	CnhatC(vecs,vect,X,n,q,p,vecd,restmp);
	res[j*N[0] + i] = restmp[0];
      }

    }

    for (i = 0; i <= (N[0]-2); i++) {
      for (j = i+1; j <= Nm; j++) {
	res[j*N[0] + i] = conj(res[i*N[0] + j]);
      }
    }

    delete[] restmp;
    delete[] vecs;
    delete[] vect;
  }

  void CnhatmatClower(double *yMat, int *N, double *X, int *n, int *q, int *p, int *vecd, double _Complex *res, double *sqrtweights) {

    int i, j, k, l, qm = q[0] - 1, Nm = N[0] - 1;
    double *vect, *vecs;
    double _Complex  *restmp;
    restmp = new _Complex double[1];
    vecs = new double[q[0]];
    vect = new double[q[0]];

    l = 0;
    for (j = 0; j <= Nm; j++) {
      for (k = 0; k <= qm; k++) {
	vect[k] = yMat[k*N[0] + j];
      }
      for (i = j; i <= Nm; i++) {
	for (k = 0; k <= qm; k++) {
	  vecs[k] = yMat[k*N[0] + i];
	}
	CnhatC(vecs,vect,X,n,q,p,vecd,restmp);
	res[l] = restmp[0] * sqrtweights[i] * sqrtweights[j];
	l = l + 1;
      }
    }

    delete[] restmp;
    delete[] vecs;
    delete[] vect;
  }
  
  void phinhatReturn(double *vect1, double *vect2, double *vect3, double *X, int *q, int *n, double _Complex *res1, double _Complex *res2, double _Complex *res3) {
    int j, col, qm = q[0]-1;
    double tmp1, tmp2, tmp3, tmp;
    res1[0] = 0.0 + 0.0*_Complex_I;
    res2[0] = 0.0 + 0.0*_Complex_I;
    res3[0] = 0.0 + 0.0*_Complex_I;

    for (j = 0; j<= (n[0]-1); j++) {
      tmp1 = 0.0;
      tmp2 = 0.0;
      tmp3 = 0.0;
      for (col = 0; col <= qm; col++) {
	tmp = X[col * n[0] + j];
	tmp1 = tmp1 + vect1[col] * tmp;
	tmp2 = tmp2 + vect2[col] * tmp;
	tmp3 = tmp3 + vect3[col] * tmp;
      }      
      res1[0] = res1[0] + cexp(_Complex_I * tmp1);
      res2[0] = res2[0] + cexp(_Complex_I * tmp2);
      res3[0] = res3[0] + cexp(_Complex_I * tmp3);
    }
    res1[0] = res1[0] / (double _Complex)n[0];
    res2[0] = res2[0] / (double _Complex)n[0];
    res3[0] = res3[0] / (double _Complex)n[0];
  }
  

}
