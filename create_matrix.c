/* Dynamicka alokace matice dimense M[m][n] */

#include <stdlib.h>
#include <stdio.h> 

#include "matrix.h"

MATRIX create_matrix(int m, int n)
{
  MATRIX M;
  double **p_p_M;
  int i;

  if ( (p_p_M = (double **) malloc(m*sizeof(double *))) == NULL) {
    fprintf(stderr,"create_matrix: No memory for row pointers allocation!\n");
    exit (EXIT_FAILURE);
  }
  for (i=0; i<m; i++) {
    if ( (p_p_M[i] = (double *) malloc(n*sizeof(double))) == NULL) {
      fprintf(stderr,"create_matrix: No memory for row allocation!\n");
      exit (EXIT_FAILURE);
    }
  }
  M.m = m;
  M.n = n;
  M.M = p_p_M;
  return M;
}
