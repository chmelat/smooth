/*  Savitzky-Golay filter for data smoothing
 *  Implementation
 *  V2.1/2025-10-04/ Added uniform grid check (SG not suitable for non-uniform grids)
 *  V2.0/2025-05-28/ Updated to use LAPACK/BLAS instead of lib_matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "savgol.h"
#include "grid_analysis.h"

#ifndef DPMAX
#define DPMAX 12
#endif

/* Threshold for grid uniformity - above this CV, grid is considered non-uniform */
#define UNIFORMITY_CV_THRESHOLD 0.05

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
    double *A;
    double *B;
    int matrix_size;
    int info;
    int nrhs = 1;
    char uplo = 'U';
    
    /* Input validation */
    if (c == NULL || poly_degree < 0 || deriv_order < 0 || deriv_order > poly_degree) {
        fprintf(stderr, "Error: Invalid parameters for savgol_coefficients\n");
        return;
    }
    
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
    for (j = 0; j <= poly_degree; j++) {
        for (i = 0; i <= j; i++) {
            A[i + j*matrix_size] = a[i + j];
        }
        
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
        double pos = i - nl;
        
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
    double h = 1.0;
    
    /* Allocate coefficients array */
    c = (double*)calloc(nl + nr + 1, sizeof(double));
    if (c == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for filter coefficients\n");
        return 0.0;
    }
    
    /* Calculate filter coefficients */
    savgol_coefficients(nl, nr, poly_degree, deriv_order, c);
    
    /* For derivatives, calculate average spacing */
    if (deriv_order > 0) {
        h = 0.0;
        int count = 0;
        for (i = center_idx - nl + 1; i <= center_idx + nr; i++) {
            if (i > center_idx - nl && i <= center_idx + nr && i > 0 && i < data_size) {
                h += x[i] - x[i-1];
                count++;
            }
        }
        if (count > 0) {
            h /= count;
        } else {
            h = 1.0;
        }
    }
    
    /* Apply the filter */
    for (i = 0; i <= nl + nr; i++) {
        int idx = center_idx - nl + i;
        
        if (idx >= 0 && idx < data_size) {
            result += c[i] * y[idx];
        }
    }
    
    /* Scale for derivatives */
    if (deriv_order > 0) {
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
    
    /* CRITICAL: Check grid uniformity - SG method requires uniform grid */
    grid_info = analyze_grid(x, n, 0);
    if (grid_info == NULL) {
        fprintf(stderr, "Error: Grid analysis failed\n");
        return NULL;
    }
    
    /* Check if grid is sufficiently uniform for Savitzky-Golay method */
    if (grid_info->cv > UNIFORMITY_CV_THRESHOLD) {
        fprintf(stderr, "\n");
        fprintf(stderr, "========================================\n");
        fprintf(stderr, "ERROR: Savitzky-Golay method not suitable for non-uniform grid!\n");
        fprintf(stderr, "========================================\n");
        fprintf(stderr, "Grid analysis:\n");
        fprintf(stderr, "  Coefficient of variation (CV) = %.4f\n", grid_info->cv);
        fprintf(stderr, "  Threshold for uniformity = %.4f\n", UNIFORMITY_CV_THRESHOLD);
        fprintf(stderr, "  h_min = %.6e, h_max = %.6e, h_avg = %.6e\n", 
                grid_info->h_min, grid_info->h_max, grid_info->h_avg);
        fprintf(stderr, "  Ratio h_max/h_min = %.2f\n", grid_info->ratio_max_min);
        fprintf(stderr, "\n");
        fprintf(stderr, "The Savitzky-Golay filter assumes uniformly spaced data points.\n");
        fprintf(stderr, "Your data has significant spacing variation.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "RECOMMENDED ALTERNATIVES:\n");
        fprintf(stderr, "  1. Use Tikhonov method: -m 2 -l auto\n");
        fprintf(stderr, "     (Works correctly with non-uniform grids)\n");
        fprintf(stderr, "  2. Use Polyfit method: -m 0 -n %d -p %d\n", window_size, poly_degree);
        fprintf(stderr, "     (Local fitting, less sensitive to spacing)\n");
        fprintf(stderr, "  3. Resample your data to uniform grid before smoothing\n");
        fprintf(stderr, "\n");
        
        free_grid_analysis(grid_info);
        return NULL;
    }
    
    /* Grid is sufficiently uniform - proceed with SG filtering */
    if (grid_info->cv > 0.01) {
        /* Grid is uniform enough but not perfectly uniform - just warn */
        printf("# Savitzky-Golay: Grid is nearly uniform (CV=%.4f)\n", grid_info->cv);
    }
    
    free_grid_analysis(grid_info);
    
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
    
    /* Left boundary handling */
    for (i = 0; i < offset; i++) {
        int left_pts = i;
        int right_pts = window_size - 1 - left_pts;
        
        if (left_pts + right_pts < poly_degree) {
            right_pts = poly_degree - left_pts;
            if (i + right_pts >= n) {
                result->y_smooth[i] = y[i];
                result->y_deriv[i] = 0.0;
                continue;
            }
        }
        
        fi = apply_savgol_filter(x, y, i, left_pts, right_pts, poly_degree, 0, n);
        dfi = apply_savgol_filter(x, y, i, left_pts, right_pts, poly_degree, 1, n);
        
        result->y_smooth[i] = fi;
        result->y_deriv[i] = dfi;
    }
    
    /* Central points */
    for (i = offset; i < n - offset; i++) {
        fi = apply_savgol_filter(x, y, i, offset, offset, poly_degree, 0, n);
        dfi = apply_savgol_filter(x, y, i, offset, offset, poly_degree, 1, n);
        
        result->y_smooth[i] = fi;
        result->y_deriv[i] = dfi;
    }
    
    /* Right boundary handling */
    for (i = n - offset; i < n; i++) {
        int right_pts = n - 1 - i;
        int left_pts = window_size - 1 - right_pts;
        
        if (left_pts + right_pts < poly_degree) {
            left_pts = poly_degree - right_pts;
            if (i - left_pts < 0) {
                result->y_smooth[i] = y[i];
                result->y_deriv[i] = 0.0;
                continue;
            }
        }
        
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
