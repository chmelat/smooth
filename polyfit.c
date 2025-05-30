/*  Polynomial fitting for data smoothing
 *  Implementation of least squares polynomial approximation
 *  V2.0/2025-05-28/ Updated to use LAPACK/BLAS instead of lib_matrix
 *  V1.0/2025-05-27/ Extracted from smooth.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "polyfit.h"

#ifndef DPMAX  
#define DPMAX 12  /* Maximum degree of approximation polynomial - configurable */
#endif

/* LAPACK function declarations */
extern void dposv_(char *uplo, int *n, int *nrhs, double *a, int *lda, 
                   double *b, int *ldb, int *info);

/* Local function declarations */
static double power(double x, int n);

/* Power function x^n */
static double power(double x, int n)
{
    double p = 1;
    
    while (n-- > 0)
        p *= x;
    
    return p;
}

/* Main polynomial fitting function */
PolyfitResult* polyfit_smooth(double *x, double *y, int n, int window_size, int poly_degree)
{
    PolyfitResult *result;
    double *C;  /* Matrix of normal equations (stored column-major for LAPACK) */
    double *B;  /* Right-hand side vector */
    double sx[DPMAX+3], sy[DPMAX+1];
    int i, j, k;
    int offset;
    int matrix_size;
    int info;
    int nrhs = 1;
    char uplo = 'U';  /* Upper triangular part */
    
    /* Input validation */
    if (x == NULL || y == NULL || n < 1) {
        fprintf(stderr, "Error: Invalid input parameters to polyfit_smooth\n");
        return NULL;
    }
    
    if (window_size < 3 || !(window_size % 2)) {
        fprintf(stderr, "Error: Window size must be odd and >= 3\n");
        return NULL;
    }
    
    if (poly_degree < 0 || poly_degree > DPMAX) {
        fprintf(stderr, "Error: Polynomial degree must be between 0 and %d\n", DPMAX);
        return NULL;
    }

    if (poly_degree > 6) {
    fprintf(stderr, "Warning: High polynomial degree (%d) may cause numerical instability\n", poly_degree);
}
    
    if (n < window_size) {
        fprintf(stderr, "Error: Not enough data points (n=%d < window_size=%d)\n", n, window_size);
        return NULL;
    }
    
    /* Check that x is monotonic */
    for (i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "Error: x array must be strictly increasing\n");
            return NULL;
        }
    }
    
    offset = window_size / 2;
    matrix_size = poly_degree + 1;
    
    /* Allocate result structure */
    result = (PolyfitResult *)malloc(sizeof(PolyfitResult));
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for PolyfitResult\n");
        return NULL;
    }
    
    result->n = n;
    result->poly_degree = poly_degree;
    result->window_size = window_size;
    result->y_smooth = (double *)calloc(n, sizeof(double));
    result->y_deriv = (double *)calloc(n, sizeof(double));
    result->coeffs = NULL;  /* Not used in current implementation */
    
    if (result->y_smooth == NULL || result->y_deriv == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for result arrays\n");
        free_polyfit_result(result);
        return NULL;
    }
    
    /* Allocate matrices for least squares */
    C = (double *)malloc(matrix_size * matrix_size * sizeof(double));
    B = (double *)malloc(matrix_size * sizeof(double));
    
    if (C == NULL || B == NULL) {
        fprintf(stderr, "Error: Matrix allocation failed\n");
        free(C);
        free(B);
        free_polyfit_result(result);
        return NULL;
    }
    
    /* Process each point in the interior */
    for (i = offset; i < n - offset; i++) {
        /* Initialize sums */
        for (k = 0; k < poly_degree + 3; k++)
            sx[k] = 0;
        for (k = 0; k < poly_degree + 1; k++)
            sy[k] = 0;
        
        /* Compute sums for normal equations */
        for (j = i - offset; j <= i + offset; j++) {
            double dx = x[j] - x[i];
            
            for (k = 0; k < poly_degree + 3; k++)
                sx[k] += power(dx, k);
            
            for (k = 0; k < poly_degree + 1; k++)
                sy[k] += y[j] * power(dx, k);
        }
        
        /* Build matrix of coefficients and right-hand side */
        /* Matrix C is symmetric, so we only need to fill upper triangle */
        /* LAPACK uses column-major storage */
        {
            int i2, j2;
            for (j2 = 0; j2 <= poly_degree; j2++) {
                for (i2 = 0; i2 <= j2; i2++) {
                    /* Column-major: C[i2,j2] = C[i2 + j2*matrix_size] */
                    C[i2 + j2*matrix_size] = sx[i2 + j2];
                }
            }
            
            /* Copy right-hand side */
            for (i2 = 0; i2 <= poly_degree; i2++)
                B[i2] = sy[i2];
        }
        
        /* Solve linear system using LAPACK */
        /* dposv solves A*X = B for symmetric positive definite matrix A */
        dposv_(&uplo, &matrix_size, &nrhs, C, &matrix_size, B, &matrix_size, &info);
        
        if (info != 0) {
            fprintf(stderr, "Error: LAPACK dposv failed with info = %d at point %d\n", info, i);
            if (info > 0) {
                fprintf(stderr, "Matrix is not positive definite (leading minor %d)\n", info);
            }
            free(C);
            free(B);
            free_polyfit_result(result);
            return NULL;
        }
        
        /* Store smoothed value and derivative at center point */
        result->y_smooth[i] = B[0];
        result->y_deriv[i] = (poly_degree > 0) ? B[1] : 0.0;
        
        /* Handle left boundary on first interior point */
        if (i == offset) {
            for (k = 0; k < offset; k++) {
                double dx = x[k] - x[offset];
                double fi = 0, dfi = 0;
                int m;
                
                for (m = 0; m <= poly_degree; m++)
                    fi += B[m] * power(dx, m);
                
                for (m = 1; m <= poly_degree; m++)
                    dfi += m * B[m] * power(dx, m - 1);
                
                result->y_smooth[k] = fi;
                result->y_deriv[k] = dfi;
            }
        }
        
        /* Handle right boundary on last interior point */
        if (i == n - offset - 1) {
            for (k = 0; k < offset; k++) {
                double dx = x[n - offset + k] - x[n - offset - 1];
                double fi = 0, dfi = 0;
                int m;
                
                for (m = 0; m <= poly_degree; m++)
                    fi += B[m] * power(dx, m);
                
                for (m = 1; m <= poly_degree; m++)
                    dfi += m * B[m] * power(dx, m - 1);
                
                result->y_smooth[n - offset + k] = fi;
                result->y_deriv[n - offset + k] = dfi;
            }
        }
    }
    
    /* Clean up */
    free(C);
    free(B);
    
    return result;
}

/* Free allocated memory */
void free_polyfit_result(PolyfitResult *result)
{
    if (result != NULL) {
        free(result->y_smooth);
        free(result->y_deriv);
        free(result->coeffs);
        free(result);
    }
}
