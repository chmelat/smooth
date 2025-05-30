/*  Savitzky-Golay filter for data smoothing
 *  Implementation
 *  V2.0/2025-05-28/ Updated to use LAPACK/BLAS instead of lib_matrix
 *  V1.1/2025-05-27/ Added grid analysis support
 *  V1.0/2025-05-27/ Extracted from smooth.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "savgol.h"
#include "grid_analysis.h"

#ifndef DPMAX
#define DPMAX 12  /* Maximum degree of approximation polynomial - configurable */  
#endif

/* LAPACK function declarations */
extern void dposv_(char *uplo, int *n, int *nrhs, double *a, int *lda, 
                   double *b, int *ldb, int *info);

/* Local function declarations */
static double power(double x, int n);
static double factorial(int n);
static double apply_savgol_filter(double *x, double *y, int center_idx, 
                                  int nl, int nr, int poly_degree, 
                                  int deriv_order, int data_size);

/* Power function x^n */
static double power(double x, int n)
{
    double p = 1;
    
    while (n-- > 0)
        p *= x;
    
    return p;
}

/* Calculate factorial */
static double factorial(int n)
{
    double result = 1.0;
    
    for (int i = 2; i <= n; i++)
        result *= i;
    
    return result;
}

/* Calculate Savitzky-Golay filter coefficients */
void savgol_coefficients(int nl, int nr, int poly_degree, int deriv_order, double *c)
{
    int i, j;
    double *a;
    double *A;  /* Matrix for normal equations (column-major storage) */
    double *B;  /* Right-hand side vector */
    int matrix_size;
    int info;
    int nrhs = 1;
    char uplo = 'U';  /* Upper triangular part */
    
    /* Input validation */
    if (c == NULL || poly_degree < 0 || deriv_order < 0 || deriv_order > poly_degree) {
        fprintf(stderr, "Error: Invalid parameters for savgol_coefficients\n");
        return;
    }
    
    /* Check that we have enough points for the polynomial degree */
    if (nl + nr < poly_degree) {
        fprintf(stderr, "Error: Not enough points for polynomial degree\n");
        return;
    }
    
    matrix_size = poly_degree + 1;
    
    /* Allocate matrices and arrays */
    A = (double*)malloc(matrix_size * matrix_size * sizeof(double));
    B = (double*)malloc(matrix_size * sizeof(double));
    a = (double*)calloc(2 * poly_degree + 1, sizeof(double));
    
    if (A == NULL || B == NULL || a == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in savgol_coefficients\n");
        free(A);
        free(B);
        free(a);
        return;
    }
    
    /* Fill 'a' array with the moments of the data positions */
    for (i = 0; i <= 2 * poly_degree; i++) {
        for (j = -nl; j <= nr; j++)
            a[i] += power(j, i);
    }
    
    /* Set up the normal equations for the desired polynomial fit */
    /* Matrix A is symmetric, so we only need to fill upper triangle */
    /* LAPACK uses column-major storage */
    for (j = 0; j <= poly_degree; j++) {
        for (i = 0; i <= j; i++) {
            /* Column-major: A[i,j] = A[i + j*matrix_size] */
            A[i + j*matrix_size] = a[i + j];
        }
        
        /* Right-hand side for the derivatives */
        if (j == deriv_order)
            B[j] = factorial(deriv_order);
        else
            B[j] = 0.0;
    }
    
    /* Solve the linear system using LAPACK */
    dposv_(&uplo, &matrix_size, &nrhs, A, &matrix_size, B, &matrix_size, &info);
    
    if (info != 0) {
        fprintf(stderr, "Error: LAPACK dposv failed with info = %d in savgol_coefficients\n", info);
        if (info > 0) {
            fprintf(stderr, "Matrix is not positive definite (leading minor %d)\n", info);
        }
        free(A);
        free(B);
        free(a);
        return;
    }
    
    /* Compute the filter coefficients using the solution */
    for (i = 0; i <= nl + nr; i++) {
        double sum = B[0];
        double pos = i - nl;  /* Relative position to center */
        
        for (j = 1; j <= poly_degree; j++)
            sum += B[j] * power(pos, j);
        
        c[i] = sum;
    }
    
    /* Clean up */
    free(a);
    free(A);
    free(B);
}

/* Apply Savitzky-Golay filter at a specific point */
static double apply_savgol_filter(double *x, double *y, int center_idx, 
                                  int nl, int nr, int poly_degree, 
                                  int deriv_order, int data_size)
{
    int i;
    double result = 0.0;
    double *c;
    double h = 1.0;  /* Normalized spacing */
    
    /* Allocate coefficients array */
    c = (double*)calloc(nl + nr + 1, sizeof(double));
    if (c == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for filter coefficients\n");
        return 0.0;
    }
    
    /* Calculate filter coefficients */
    savgol_coefficients(nl, nr, poly_degree, deriv_order, c);
    
    /* Calculate the average spacing for non-uniform grids */
    if (deriv_order > 0) {
        /* For derivatives, we need to account for non-uniform spacing */
        h = 0.0;
        int count = 0;
        for (i = center_idx - nl + 1; i <= center_idx + nr; i++) {
            if (i > center_idx - nl && i <= center_idx + nr && i > 0 && i < data_size) {
                h += x[i] - x[i-1];
                count++;
            }
        }
        /* Avoid division by zero */
        if (count > 0) {
            h /= count;  /* Average spacing */
        } else {
            h = 1.0;  /* Default if no valid intervals */
        }
    }
    
    /* Apply the filter */
    for (i = 0; i <= nl + nr; i++) {
        int idx = center_idx - nl + i;
        
        /* Boundary check */
        if (idx >= 0 && idx < data_size) {
            result += c[i] * y[idx];
        }
    }
    
    /* Scale for derivatives - for non-uniform x values */
    if (deriv_order > 0) {
        /* Convert from normalized derivatives to actual derivatives */
        double scale_factor = 1.0;
        for (i = 0; i < deriv_order; i++)
            scale_factor /= h;
        
        result *= scale_factor;
    }
    
    free(c);
    return result;
}

/* Main Savitzky-Golay smoothing function */
SavgolResult* savgol_smooth(double *x, double *y, int n, int window_size, int poly_degree)
{
    SavgolResult *result;
    GridAnalysis *grid_info;
    int i;
    int offset;
    double fi, dfi;
    
    /* Input validation */
    if (x == NULL || y == NULL || n < 1) {
        fprintf(stderr, "Error: Invalid input parameters to savgol_smooth\n");
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

    if (poly_degree >= window_size) {
        fprintf(stderr, "Error: Polynomial degree must be less than window size\n");
        return NULL;
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
    
    /* Analyze grid uniformity */
    grid_info = analyze_grid(x, n, 0);
    if (grid_info == NULL) {
        fprintf(stderr, "Warning: Grid analysis failed, proceeding with standard method\n");
    } else {
        /* Print grid analysis if significant non-uniformity detected */
        if (grid_info->reliability_warning) {
            printf("# Savitzky-Golay: Grid analysis detected non-uniformity\n");
            print_grid_analysis(grid_info, 0, "# ");
            printf("# %s\n", get_grid_recommendation(grid_info));
            
            /* Suggest adaptive window size for highly non-uniform grids */
            if (grid_info->ratio_max_min > 20.0) {
                int suggested_window = optimal_window_size(grid_info, 3, window_size);
                if (suggested_window < window_size) {
                    printf("# WARNING: Consider using smaller window size (%d) for this non-uniform grid\n", 
                           suggested_window);
                }
            }
        }
        free_grid_analysis(grid_info);
    }
    
    offset = window_size / 2;
    
    /* Allocate result structure */
    result = (SavgolResult *)malloc(sizeof(SavgolResult));
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for SavgolResult\n");
        return NULL;
    }
    
    result->n = n;
    result->poly_degree = poly_degree;
    result->window_size = window_size;
    result->y_smooth = (double *)calloc(n, sizeof(double));
    result->y_deriv = (double *)calloc(n, sizeof(double));
    
    if (result->y_smooth == NULL || result->y_deriv == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for result arrays\n");
        free_savgol_result(result);
        return NULL;
    }
    
    /* Left boundary handling - special case with reduced window */
    for (i = 0; i < offset; i++) {
        /* For left boundary points, use an asymmetric window */
        int left_pts = i;
        int right_pts = window_size - 1 - left_pts;
        
        /* Ensure we don't exceed polynomial degree constraints */
        if (left_pts + right_pts < poly_degree) {
            right_pts = poly_degree - left_pts;
            if (i + right_pts >= n) {
                /* Can't fit polynomial, use nearest valid point */
                result->y_smooth[i] = y[i];
                result->y_deriv[i] = 0.0;
                continue;
            }
        }
        
        /* Calculate smoothed value and derivative */
        fi = apply_savgol_filter(x, y, i, left_pts, right_pts, poly_degree, 0, n);
        dfi = apply_savgol_filter(x, y, i, left_pts, right_pts, poly_degree, 1, n);
        
        result->y_smooth[i] = fi;
        result->y_deriv[i] = dfi;
    }
    
    /* Central points - use full symmetric window */
    for (i = offset; i < n - offset; i++) {
        /* Calculate smoothed value and derivative */
        fi = apply_savgol_filter(x, y, i, offset, offset, poly_degree, 0, n);
        dfi = apply_savgol_filter(x, y, i, offset, offset, poly_degree, 1, n);
        
        result->y_smooth[i] = fi;
        result->y_deriv[i] = dfi;
    }
    
    /* Right boundary handling - special case with reduced window */
    for (i = n - offset; i < n; i++) {
        /* For right boundary points, use an asymmetric window */
        int right_pts = n - 1 - i;
        int left_pts = window_size - 1 - right_pts;
        
        /* Ensure we don't exceed polynomial degree constraints */
        if (left_pts + right_pts < poly_degree) {
            left_pts = poly_degree - right_pts;
            if (i - left_pts < 0) {
                /* Can't fit polynomial, use nearest valid point */
                result->y_smooth[i] = y[i];
                result->y_deriv[i] = 0.0;
                continue;
            }
        }
        
        /* Calculate smoothed value and derivative */
        fi = apply_savgol_filter(x, y, i, left_pts, right_pts, poly_degree, 0, n);
        dfi = apply_savgol_filter(x, y, i, left_pts, right_pts, poly_degree, 1, n);
        
        result->y_smooth[i] = fi;
        result->y_deriv[i] = dfi;
    }
    
    return result;
}

/* Free allocated memory */
void free_savgol_result(SavgolResult *result)
{
    if (result != NULL) {
        free(result->y_smooth);
        free(result->y_deriv);
        free(result);
    }
}
