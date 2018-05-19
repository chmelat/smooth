/*
 * Free matrix[m][n]
 *
 * V0/14.5.2018/
 */

#include <stdlib.h>
#include "matrix.h"

void free_matrix(MATRIX M)
{
  int i;
  double **p_p_M = M.M;

  for (i=0; i<M.m; i++)
     free(p_p_M[i]);
  free(p_p_M);
}
