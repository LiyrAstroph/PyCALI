/*
 *  PIXON
 *  A Pixon-based method for reconstructing velocity-delay map in reverberation mapping.
 * 
 *  Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

/*!
 *  \file mathfun.h
 *  \brief header file for mathfun.c.
 */

#ifndef _MATHFUN_H

#define _MATHFUN_H


#ifdef __cplusplus
extern "C" {
#endif


/* matrix operations */
void inverse_mat(double *a, int n, int *info);
void inverse_symat(double *a, int n, int *info);
void inverse_pomat(double *a, int n, int *info);
double det_mat(double *a, int n, int *info);
double lndet_mat(double *a, int n, int *info);
void display_mat(double *a, int m, int n);
void multiply_mat(double * a, double *b, double *c, int n);
void multiply_mat_transposeA(double * a, double *b, double *c, int n);
void multiply_mat_transposeB(double * a, double *b, double *c, int n);
void multiply_matvec(double *a, double *x, int n, double *y);
void multiply_matvec_transposeA(double *a, double *x, int n, double *y);
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y);
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k);
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n);
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
void Chol_decomp_L(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);
void test_mathfun();
void compute_semiseparable_drw(double *t, int n, double a1, double c1, double *sigma, double syserr,  double *W, double *D, double *phi);
void multiply_matvec_semiseparable_drw(double *y, double  *W, double *D, double *phi, int n, double a1, double *z);
void multiply_mat_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z);
void multiply_mat_transposeB_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z);

void inverse_semiseparable_uv(double *t, int n, double a1, double c1, double *A);

void multiply_matvec_semiseparable_uv(double *y, double *u, double  *W, double *D, double *phi,
                                      int n, double *z);

void compute_inverse_semiseparable_plus_diag(double *t, int n, double a1, double c1, 
                double *sigma, double syserr, double *u, double *v,
                double *W, double *D, double *phi, double *work);

int compare(const void* a, const void* b);

#ifdef __cplusplus
}
#endif

#endif