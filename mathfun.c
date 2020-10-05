/*
 *  PIXON
 *  A Pixon-based method for reconstructing velocity-delay map in reverberation mapping.
 * 
 *  Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

/*!
 *  \file mathfun.c
 *  \brief mathematic functions.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <lapacke.h>
#include <cblas.h>

#include "mathfun.h"

#define EPS (DBL_MIN)

/*!
 * This function calculates matrix multiply C(nxn) = A(nxn) * B(nxn).
 */
void multiply_mat(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

/*!
 * This function calculates matrix multiply C(nxn) = A^T(nxn) * B(nxn).
 */
void multiply_mat_transposeA(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

/*!
 * This function calculates matrix multiply C(nxn) = A(nxn) * B^T(nxn).
 */
void multiply_mat_transposeB(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

/*!
 * This function calculates matrix multiply Y(n) = A^T(nxn) * X(n).
 */
void multiply_matvec_transposeA(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/*!
 * This function calculates matrix multiply Y(n) = A(nxn) * X(n).
 */
void multiply_matvec(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/*!
 * This function calculates  Y(m) = A(m, n) * X(n).
 */
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/* C(m*n) = A(m*k) * B(k*n) */
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, k, b, n, 0.0f, c, n);
}

/* C(m*n) = A^T(m*k) * B(k*n) */
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, m, b, n, 0.0f, c, n);
}
/* C(m*n) = A(m*k) * B^T(k*n) */
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0f
                             , a, k, b, k, 0.0f, c, n);
}

/* A(mxm)^-1 * B(mxn), store the output in B
 * note that A will be changed on exit. */
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n)
{
  int * ipiv, info;
  ipiv=malloc(m*sizeof(int));

  info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, a, m, ipiv);
  if(info!=0)
  {
    printf("multiply_mat_MN_inverseA 1.\n this usually caused by improper nc.\n increase the low limit of nc");
    exit(0);
    return info;
  }
  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', m, n, a, m, ipiv, b, n);
  if(info!=0)
  {
    printf("multiply_mat_MN_inverseA 2\n this usually caused by improper nc.\n increase the low limit of nc");
    exit(0);
    return info;
  }

  free(ipiv);
  return info;
}

/*!
 * This functions calculate A^-1(nxn).
 * A is a generic matrix.
 */
void inverse_mat(double * a, int n, int *info)
{
  int * ipiv;
  ipiv=malloc(n*sizeof(int));

//  dgetrf_(&n, &n, a, &n, ipiv, info);
//  dgetri_(&n, a, &n, ipiv, work, &lwork, info);

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);

  free(ipiv);
  return;
}

/*!
 * This functions calculate A^-1(nxn).
 * A is a symmetric matrix.
 */
void inverse_symat(double * a, int n, int *info)
{
  int * ipiv, i, j;
  ipiv=malloc(n*sizeof(int));

  LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);
  LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);

  /* fill up the lower triangle */
  for(i=0; i<n; i++)
    for(j=0; j<i; j++)
      a[i*n+j] = a[j*n+i];

  free(ipiv);
  return;
}

/*!
 * This functions calculate A^-1(nxn).
 * A is a postive-definite symmetrix matrix.
 */
void inverse_pomat(double * a, int n, int *info)
{
  int i, j;

  LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, a, n);
  LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', n, a, n);

  /* fill up the lower triangle */
  for(i=0; i<n; i++)
    for(j=0; j<i; j++)
      a[i*n+j] = a[j*n+i];
  return;
}

/*!
 * This function calculate eigenvectors and eigenvalues of matrix.
 */
void eigen_sym_mat(double *a, int n, double *val, int *info)
{
    char jobz='V', uplo='U';

/* store the eigenvectors  in a by rows.
 * store the eigenvalues in val in ascending order.
 */
//    dsyev_(&jobz, &uplo, &n, a, &n, val, work, &lwork, info);

/* store the eigenvectors  in a by columns.
 * store the eigenvalues in val in ascending order.
 */
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, val);
    return;
}

/*!
 * This function calculate A(nxn) = X^T(1xn)*X(1xn)
 */
void multiply_vec2mat(double * x, double * a, int n)
{
//  cblas_dsyr(CblasRowMajor, CblasUpper, n, 1.0f, x, 1, a, n);
  int i, j;
  for(i=0; i<n; i++)
    for(j=0; j<=i; j++)
    {
      a[i*n+j] = a[j*n+i] = x[i]*x[j];
    }
}

/*!
 * This function calculates determinant of matrix A. \n
 * There are two versions in the internet, the main difference lies at dealing with the sign of det. \n
 * The version II is verifed to be \b INCORRECT. \n
 * Note that LAPACK is written in Fortran, the indix diffes by 1 with that in C.
 * LAPACK version 3.5.0
 * ********************************
 * Version I:
 * \code{.sh}
 *  det = 1.0; 
 *  for(i=0; i<n; i++)
 *  {
 *   det *= a[i*n+i];
 *   if (ipiv[i] != i+1) 
 *   {
 *     det = -det;
 *   }
 *  }
 * \endcode
 * ********************************
 * Version II:
 * \code{.sh}
 *  det = 1.0;
 *  for(i=0; i<n; i++)
 *  {
 *   det *= a[i*n+i];
 *   if (ipiv[i] != i+1) 
 *   {
 *     ipiv[ipiv[i]-1] = ipiv[i];
 *     det = -det;
 *   }
 *  }
 * \endcode
 * ********************************
 */
double det_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double det;
  ipiv=malloc(n*sizeof(int));
  
  /* LU decomposition */
//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    printf("# Error, Wrong in det_mat!\n");
    exit(-1);
  }

  det = 1.0;
  for(i=0; i<n; i++)
  {
    printf("%d\n", ipiv[i]);
    det *= a[i*n+i];
    if (ipiv[i] != i+1) // note that LAPACK is written in C, the indix diffes by 1 with C.
    {
      det = -det;
    }
  }
  //det=fabs(det);
  free(ipiv);
  return det;
}

/*!
 * This function calculates logarithm determinant of matrix.
 */
double lndet_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double lndet;
  ipiv=malloc(n*sizeof(int));

  /* LU factorization */
//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    printf("Wrong!\n");
    exit(-1);
  }

  lndet = 0.0;
  for(i=0; i<n; i++)
  {
    lndet += log(a[i*n+i]);
  }
  free(ipiv);
  return lndet;
}

/*!
 * This function performs Cholesky decomposition of matrix into upper triangle matrixs
 */
void Chol_decomp_U(double *a, int n, int *info)
{
  int i,j;
  char uplo = 'U';   // decomposite as A = U^T*U, i.e., upper triangle matrix.
//  dpotrf_(&uplo, &n, a, &n, info);
  *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  if(*info<0)
  {
    fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
//    exit(-1);
    return;
  }
  else if (*info>0)
  {
    fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
//    exit(-1);
    return;
  }

  // only the upper triangle is referenced by dpotrf output, 
  // so the strictly lower triangle must be set to zero
  for(i=0;i<n;i++)
    for(j=0;j<i;j++)
      a[i*n+j] = 0.0;
  return;
}

/*!
 * This function performs Cholesky decomposition of matrix into lower triangle matrixs
 */
void Chol_decomp_L(double *a, int n, int *info)
{
  int i,j;
  char uplo = 'L';   // decomposite as A = L*L^T, i.e., lower triangle matrix.
//  dpotrf_(&uplo, &n, a, &n, info);
  *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  if(*info<0)
  {
    fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
//    exit(-1);
    return;
  }
  else if (*info>0)
  {
    fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
//    exit(-1);
    return;
  }

  // only the lower triangle is referenced by dpotrf output, 
  // so the strictly upper triangle must be set to zero
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      a[i*n+j] = 0.0;
  return;
}

/*
 * semiseparable matrix
 */
void compute_semiseparable_drw(double *t, int n, double a1, double c1, double *sigma, double syserr, double *W, double *D, double *phi)
{
  int i;
  double S, A;
  phi[0] = 0.0;
  for(i=1; i<n; i++)
  {
    phi[i] = exp(-c1 * (t[i] - t[i-1]));
  }

  S = 0.0;
  A = sigma[0]*sigma[0] + syserr*syserr + a1;
  D[0] = A;
  W[0] = 1.0/D[0];
  for(i=1; i<n; i++)
  {
    S = phi[i]*phi[i] * (S + D[i-1]*W[i-1]*W[i-1]);
    A = sigma[i]*sigma[i] + syserr*syserr + a1;
    D[i] = A - a1 * a1 * S;
    W[i] = 1.0/D[i] * (1.0 - a1*S);
  }
}
/*
 * z = C^-1 x y
 *
 * y is a vector
 */
void multiply_matvec_semiseparable_drw(double *y, double  *W, double *D, double *phi, int n, double a1, double *z)
{
  int i;
  double f, g;

  // forward substitution
  f = 0.0;
  z[0] = y[0];
  for(i=1; i<n;i++)
  {
    f = phi[i] * (f + W[i-1] * z[i-1]);
    z[i] = y[i] - a1*f;
  }

  //backward substitution
  g = 0.0;
  z[n-1] = z[n-1]/D[n-1];
  for(i=n-2; i>=0; i--)
  {
    g = phi[i+1] *(g + a1*z[i+1]);
    z[i] = z[i]/D[i] - W[i]*g;
  }
}
/*
 * Z = C^-1 x Y
 * 
 * Y is an (nxm) matrix. 
 * Note that Y is row-major
 */
void multiply_mat_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z)
{
  int i, j;
  double f, g;

  // forward substitution
  for(j=0; j<m; j++)
  {
    f = 0.0;
    Z[0*m+j] = Y[0*m+j];
    for(i=1; i<n;i++)
    {
      f = phi[i] * (f + W[i-1] * Z[(i-1)*m + j]);
      Z[i*m+j] = Y[i*m+j] - a1*f;
    }
  }

  //backward substitution
  for(j=0; j<m; j++)
  {
    g = 0.0;
    Z[(n-1)*m+j] = Z[(n-1)*m+j]/D[n-1];
    for(i=n-2; i>=0; i--)
    {
      g = phi[i+1] *(g + a1*Z[(i+1)*m+j]);
      Z[i*m+j] = Z[i*m+j]/D[i] - W[i]*g;
    }
  }
}

/*
 * Z = C^-1 x Y^T
 * 
 * Y is an (mxn) matrix. 
 * Note that Y is row-major
 */
void multiply_mat_transposeB_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z)
{
  int i, j;
  double f, g;

  // forward substitution
  for(j=0; j<m; j++)
  {
    f = 0.0;
    Z[0*m+j] = Y[0+j*n];
    for(i=1; i<n;i++)
    {
      f = phi[i] * (f + W[i-1] * Z[(i-1)*m + j]);
      Z[i*m+j] = Y[i+j*n] - a1*f;
    }
  }

  //backward substitution
  for(j=0; j<m; j++)
  {
    g = 0.0;
    Z[(n-1)*m+j] = Z[(n-1)*m+j]/D[n-1];
    for(i=n-2; i>=0; i--)
    {
      g = phi[i+1] *(g + a1*Z[(i+1)*m+j]);
      Z[i*m+j] = Z[i*m+j]/D[i] - W[i]*g;
    }
  }
}

/**
 *  calculate A^-1.
 * 
 *  A = LxDxL^T, L = I + tril(UxW^T), D is a diagonal matrix.
 * 
 *  M = LxD^1/2,  A = MxM^T,  A^-1 = (M^T)^-1xM^-1.
 */
void inverse_semiseparable(double *t, int n, double a1, double c1, double *sigma, 
                           double syserr, double *W, double *D, double *phi,
                           double *A, double *lndet)
{
  int i, j;

  compute_semiseparable_drw(t, n, a1, c1, sigma, syserr, W, D, phi);
  
  *lndet = 0.0;
  for(i=0; i<n; i++)
  {
    A[i*n + i] = 1.0 * sqrt(D[i]);
    for(j=0; j<i; j++)
    {
      A[i*n + j] = a1 * (exp(-c1*(t[i]-t[j]))*W[j]) * sqrt(D[j]);
    }

    *lndet += log(D[i]);
  }
  LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'N', n, A, n);

  LAPACKE_dlauum(LAPACK_ROW_MAJOR, 'L', n, A, n);

  /* fill up upper triangle */
  for(i=0; i<n; i++)
    for(j=i+1; j<n; j++)
      A[i*n+j] = A[j*n+i];
  return;
}

/**
 *  caclulate inverse of a DRW semiseparable matrix, which 
 *  is a tridiagonal matrix.
 */
void inverse_semiseparable_uv(double *t, int n, double a1, double c1, double *A)
{
  int i, j;
  double b1, b2, dt;

  dt = t[1] - t[0];
  b1 = 1.0/( a1 * (exp(-c1*dt) - exp(c1*dt)) );
  A[0] = - b1 * exp(c1*dt);
  A[1] = A[1*n+0] = b1;
  for(i=2; i<n; i++)
    A[i] = 0.0;

  for(i=1; i<n-1; i++)
  {
    for(j=0; j<i-1; j++)
      A[i*n+j] = A[j*n+i] = 0.0;
    
    for(j=i+1; j<n; j++)
      A[i*n+j] = A[j*n+i] = 0.0;
    
    dt = t[i+1] - t[i];
    b2 = 1.0/( a1 * (exp(-c1*dt) - exp(c1*dt)) );
    dt = t[i+1] - t[i-1];
    A[i*n+i] = -b1 * b2 *  a1 * (exp(-c1*dt) - exp(c1*dt));
    
    A[i*n+(i+1)] = A[(i+1)*n + i] = b2;

    b1 = b2;
  }

  i = n-1;
  dt = t[i] - t[i-1];
  A[i*n+i] = - b1 * exp(c1 * dt);
  
  for(j=0; j<n-2; j++)
    A[i*n+j] = 0.0;
  return;
}

/**
 *  calculate Q = [S^-1+N^-1]^-1, where S is a DRW symmetric semiseparable matrix 
 *  and N is a diagonal matrix.
 *  
 *  S^-1 + N^-1 is a tridiagonal symmetric matrix.
 * 
 *  [S^-1 + N^-1]^-1 is a semiseparable symmetric matrix.
 * 
 */
void compute_inverse_semiseparable_plus_diag(double *t, int n, double a1, double c1, 
                double *sigma, double syserr, double *u, double *v,
                double *W, double *D, double *phi, double *work)
{
  int i;
  double *a, *b, f1, f2;

  a = work;
  b = work + n;

  phi[0] = 0.0;
  for(i=1; i<n; i++)
  {
    phi[i] = exp(-c1 * (t[i] - t[i-1]));
  }

  /* first S^-1 + N^-1 */
  f1 = (phi[1]*phi[1] - 1.0) + EPS;
  b[0] = phi[1]/( a1 * f1);
  a[0] = - 1.0/( a1 * f1) + 1.0/(sigma[0]*sigma[0] + syserr*syserr);

  for(i=1; i<n-1; i++)
  { 
    f2 = (phi[i+1]*phi[i+1] -1.0) + EPS;
    b[i] = phi[i+1]/( a1 *  f2);
    a[i] = -1.0/a1 * (1.0 + 1.0/f1 + 1.0/f2)  + 1.0/(sigma[i]*sigma[i] + syserr*syserr);

    f1 = f2;
  }
  i = n-1;
  a[i] = - 1.0 / (a1 * f1) + 1.0/(sigma[i]*sigma[i] + syserr*syserr);;
  
  /* now inverse [S^-1+N^-1]^-1*/
  v[0] = 1.0 ;
  v[1] = -a[0]/b[0] * phi[1];
  for(i=2; i<n; i++)
  {
    v[i] = -(a[i-1]*v[i-1] * phi[i]  + b[i-2]*v[i-2] * phi[i]*phi[i-1] )/b[i-1];
  }

  u[n-1] = 1.0/(b[n-2]*v[n-2] * phi[n-1] + a[n-1]*v[n-1]);
  for(i=n-2; i>0; i--)
  {
    u[i] = (1.0 - b[i]*v[i]*u[i+1] * phi[i+1] )/(a[i]*v[i]+b[i-1]*v[i-1]*phi[i]);
  }
  u[0] = (1.0-b[0]*v[0]*u[1] * phi[1] )/(a[0]*v[0]);

  /* calculate W, D */
  D[0] = u[0]*v[0];
  W[0] = 1.0/D[0];
  for(i=1; i<n; i++)
  {
    /*S = phi[i]*phi[i]*(S + D[i-1]*W[i-1]*W[i-1]);
    A = u[i]*v[i];
    D[i] = A - u[i]*u[i] * S;
    W[i] = 1.0/D[i] * (v[i] - u[i]*S);*/

    D[i] = -u[i]/u[i-1]/b[i-1] * phi[i];
    W[i] = 1.0/u[i];
    //printf("%e %e %e %e %e %e\n", u[i], v[i], D[i], W[i], a[i], b[i]);
  }
  return;
}

/*
 * z = C^1/2 x y
 *
 * C is a semiseparable matrix with (u,v) representation, y is a vector
 * 
 * C = LxDxL^T,   L=I + tril(UxW^T)
 * C^1/2 = LxD^1/2
 */
void multiply_matvec_semiseparable_uv(double *y, double *u, double  *W, double *D, double *phi,
                                      int n, double *z)
{
  int i;
  double f;

  f = 0.0;
  z[0] = sqrt(D[0]) * y[0];
  for(i=1; i<n; i++)
  {
    f = phi[i] * (f + W[i-1] * sqrt(D[i-1]) * y[i-1]);
    z[i] = sqrt(D[i]) * y[i] + u[i] * f;
  }
  return;
}


/*!
 * This function display matrix on the screen.
 */
void display_mat(double *a, int m, int n)
{
    int i, j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%e\t", a[i*n+j]);
        }
        printf("\n");
    }
}

/*!
 * This function allocates memory for matrix.
 */
double ** matrix_malloc(int n1, int n2)
{
  double ** mat;
  int i;

  if(!(mat = malloc(n1*sizeof(double*))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    exit(-1);
  }

  for(i=0; i<n1; i++)
  {
    if(!(mat[i] = malloc(n2*sizeof(double))))
    {
      fprintf(stderr, "Unable to allocate the matrix!\n");
      exit(-1);
    }
  }
  return mat;
}

/*!
 * This function allocates memory for array.
 */
double * array_malloc(int n)
{
  double *array;

  if(!(array = malloc(n*sizeof(double))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    exit(-1);
  }

  return array;
}

/*!
 * This function is to test the functions defined in mathfun.c.
 */
void test_mathfun()
{
  double *A, det;
  int n = 3, info;
  A = malloc(n*n*sizeof(double));

  A[0*n+0] = 1.0;
  A[0*n+1] = 10.3;
  A[0*n+2] = 6.3;

  A[1*n+0] = -8.3;
  A[1*n+1] = -8.0;
  A[1*n+2] = 5.3;

  A[2*n+0] = 0.3;
  A[2*n+1] = -9.3;
  A[2*n+2] = -3.0;
  
  det = det_mat(A, n, &info);
  printf("%f\n", det);
}

/* used for sorting */
int compare(const void* a, const void* b)
{
	if (*(double *)a > *(double *)b)
		return 1;
	else if (*(double *)a > *(double *)b)
		return -1;
	return 0;
}