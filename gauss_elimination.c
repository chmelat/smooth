/*
 * Reseni soustavy linearnich rovnic s matici koeficientu MA a matici
 * pravych stran MB Gaussovou eliminaci s vyberem nejvetsiho pivotu
 * Pro vypocet inversni matice je tedy treba zadat B = I, vysledkem bude
 * inversni matice v B a I v A.
 * V1/14.5.2018
 */

#include <stdlib.h>
#include <stdio.h> 
#include <math.h> 

#include "matrix.h"

#define D 0.01  /* Chyba jednicky v Identity matrix */

/* Local declare function */
static double pivot(double **M, int m, int n, int *pi, int *pj);
static void row_divide(double *V, int n, double k);
static void row_subtract(double **M, double **B, int n, int nb, int m, int pj);
static void row_swap(double **M, int i, int j);


/* Compute inversion mantrix by Gauss elimination */
void gauss_elimination(MATRIX MA, MATRIX MB)
{
  double **A, **B;
  int m,n,nb;
  int pi,pj; /* Pivots coordinate */
  double p; /* Pivot */

  if (MA.m != MA.n) {
    fprintf(stderr,"gauss_eleimination: Matrix A isn't square (%d,%d)!\n",MA.m,MA.n);
    exit(EXIT_FAILURE);
  }

  if (MA.n != MB.m) {
    fprintf(stderr,"gauss_eleimination: Number of rows B isn't number of colums A!\n");
    exit(EXIT_FAILURE);
  }

  n = MA.m; /* Dimension A */
  nb = MB.n; /* Number of right sides B */
  A = MA.M;
  B = MB.M;

  for (m=n; m>0; m--) {
    p = pivot(A,m,n,&pi,&pj); /* Max pivot */

    row_divide(A[pi],n,p); /* Deleni radku 'pi' pivotem 'p' */
    row_divide(B[pi],nb,p);

    row_swap(A,m-1,pi); /* Radku s pivotem zaradit na konec */
    row_swap(B,m-1,pi);

    row_subtract(A,B,n,nb,m,pj);
  }

/* Prerovnani radek I => Inverse matrix */
  {
    int i,j;

    for (i=0; i<n-1; i++) {
      if (fabs(A[i][i]-1) < D) {
        continue;
      }
      for (j=i+1; j<n; j++) {            // prochazej sloupec
        if (fabs(A[j][i]-1) < D) {
          row_swap(A,j,i);
          row_swap(B,j,i);
          break;
        }
      }
    }
  }
}  

/* 
 * Local functions 
 */

/* Coordinates of Maximal absolute element of matrix M[m][n] */
static double pivot(double **M, int m, int n, int *pi, int *pj)
{
  int i,j;
  double max = 0;

  for (i=0; i<m; i++)
   for (j=0; j<n; j++)
     if (fabs(M[i][j]) > max) {
       max = fabs(M[i][j]);
       *pi = i;
       *pj = j;  
     }
  return M[*pi][*pj];
}

/* deleni radku 'n' matice 'V' konstantou k */
static void row_divide(double *V, int n, double k)
{
  int i;

  if (k == 0) {
    fprintf(stderr,"row_divide: Divide by zero!\n");
    exit(EXIT_FAILURE);
  }

  for (i=0; i<n; i++) 
    *V++ /= k;

}

/* Odecet A[i][pj] nasobku posledniho radku (m-1) ode vsech ostatnich  */
static void row_subtract(double **M, double **B, int n, int nb, int m, int pj)
{
  int i,j;
  double c;

  for (i=0; i<n; i++) {
    if (i != (m-1)) {
      c = M[i][pj];
      for (j=0; j<n; j++) {
        M[i][j] -= c*M[m-1][j];
        if (j<nb)
          B[i][j] -= c*B[m-1][j];
      }
    }
  }
}

/* Prohodi i-ty a j-ty radek matice  */
static void row_swap(double **M, int i, int j)
{
  double *p_V = M[i];  

  M[i] = M[j];
  M[j] = p_V;
}
