/*  Savitzky-Golay filter for data smoothing
 *  OPTIMIZED Implementation with pre-computed coefficients
 *  V2.2/2025-10-14/ FIXED: Pre-compute coefficients once, massive speedup!
 *  V2.1/2025-10-04/ Added uniform grid check
 *  V2.0/2025-05-28/ Updated to use LAPACK/BLAS
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

/* Threshold for grid uniformity */
#define UNIFORMITY_CV_THRESHOLD 0.05

/* LAPACK function declarations */
extern void dposv_(char *uplo, int *n, int *nrhs, double *a, int *lda, 
                   double *b, int *ldb, int *info);

/* Local function declarations */
static double power(double x, int n);
static double factorial(int n);

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

/* Main Savitzky-Golay smoothing function - OPTIMIZED */
SavgolResult* savgol_smooth(double *x, double *y, int n, int window_size, int poly_degree,
                            GridAnalysis *grid_info)
{
    SavgolResult *result;
    int i, k;
    int offset;
    double h_avg;
    
    /* Pre-computed coefficient arrays */
    double *c_func = NULL;    /* Coefficients for function (deriv_order=0) */
    double *c_deriv = NULL;   /* Coefficients for derivative (deriv_order=1) */
    
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

    /* Check grid_info is available */
    if (grid_info == NULL) {
        fprintf(stderr, "Error: Grid info not available\n");
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

        return NULL;
    }

    /* Grid is sufficiently uniform - proceed with SG filtering */
    if (grid_info->cv > 0.01) {
        /* Grid is uniform enough but not perfectly uniform - just warn */
        printf("# Savitzky-Golay: Grid is nearly uniform (CV=%.4f)\n", grid_info->cv);
    }

    h_avg = grid_info->h_avg;
    
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
    
    /* ========================================================================
     * KEY OPTIMIZATION: Pre-compute coefficients ONCE for central points
     * ======================================================================== */
    
    c_func = (double *)calloc(window_size, sizeof(double));
    c_deriv = (double *)calloc(window_size, sizeof(double));
    
    if (c_func == NULL || c_deriv == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for coefficient arrays\n");
        free(c_func);
        free(c_deriv);
        free_savgol_result(result);
        return NULL;
    }
    
    /* Compute coefficients ONCE for symmetric window (central points) */
    savgol_coefficients(offset, offset, poly_degree, 0, c_func);   /* Function */
    savgol_coefficients(offset, offset, poly_degree, 1, c_deriv); /* Derivative */
    
    /* ========================================================================
     * Process central points using pre-computed coefficients - FAST!
     * ======================================================================== */
    
    for (i = offset; i < n - offset; i++) {
        double val = 0.0;
        double deriv = 0.0;
        
        /* Apply convolution with pre-computed coefficients */
        for (k = 0; k < window_size; k++) {
            int idx = i - offset + k;
            val += c_func[k] * y[idx];
            deriv += c_deriv[k] * y[idx];
        }
        
        result->y_smooth[i] = val;
        
        /* Scale derivative by average spacing */
        result->y_deriv[i] = deriv / h_avg;
    }
    
    /* ========================================================================
     * Boundary handling - compute coefficients for asymmetric windows
     * (Only a few points, so per-point computation is acceptable)
     * ======================================================================== */
    
    /* Left boundary */
    for (i = 0; i < offset; i++) {
        int left_pts = i;
        int right_pts = window_size - 1 - left_pts;
        double *c_bound_func, *c_bound_deriv;
        int n_coeff;
        
        /* Ensure enough points for polynomial degree */
        if (left_pts + right_pts < poly_degree) {
            right_pts = poly_degree - left_pts;
            if (i + right_pts >= n) {
                result->y_smooth[i] = y[i];
                result->y_deriv[i] = 0.0;
                continue;
            }
        }
        
        n_coeff = left_pts + right_pts + 1;
        c_bound_func = (double *)calloc(n_coeff, sizeof(double));
        c_bound_deriv = (double *)calloc(n_coeff, sizeof(double));
        
        if (c_bound_func && c_bound_deriv) {
            savgol_coefficients(left_pts, right_pts, poly_degree, 0, c_bound_func);
            savgol_coefficients(left_pts, right_pts, poly_degree, 1, c_bound_deriv);
            
            double val = 0.0;
            double deriv = 0.0;
            
            for (k = 0; k < n_coeff; k++) {
                int idx = i - left_pts + k;
                if (idx >= 0 && idx < n) {
                    val += c_bound_func[k] * y[idx];
                    deriv += c_bound_deriv[k] * y[idx];
                }
            }
            
            result->y_smooth[i] = val;
            result->y_deriv[i] = deriv / h_avg;
        } else {
            result->y_smooth[i] = y[i];
            result->y_deriv[i] = 0.0;
        }
        
        free(c_bound_func);
        free(c_bound_deriv);
    }
    
    /* Right boundary */
    for (i = n - offset; i < n; i++) {
        int right_pts = n - 1 - i;
        int left_pts = window_size - 1 - right_pts;
        double *c_bound_func, *c_bound_deriv;
        int n_coeff;
        
        /* Ensure enough points for polynomial degree */
        if (left_pts + right_pts < poly_degree) {
            left_pts = poly_degree - right_pts;
            if (i - left_pts < 0) {
                result->y_smooth[i] = y[i];
                result->y_deriv[i] = 0.0;
                continue;
            }
        }
        
        n_coeff = left_pts + right_pts + 1;
        c_bound_func = (double *)calloc(n_coeff, sizeof(double));
        c_bound_deriv = (double *)calloc(n_coeff, sizeof(double));
        
        if (c_bound_func && c_bound_deriv) {
            savgol_coefficients(left_pts, right_pts, poly_degree, 0, c_bound_func);
            savgol_coefficients(left_pts, right_pts, poly_degree, 1, c_bound_deriv);
            
            double val = 0.0;
            double deriv = 0.0;
            
            for (k = 0; k < n_coeff; k++) {
                int idx = i - left_pts + k;
                if (idx >= 0 && idx < n) {
                    val += c_bound_func[k] * y[idx];
                    deriv += c_bound_deriv[k] * y[idx];
                }
            }
            
            result->y_smooth[i] = val;
            result->y_deriv[i] = deriv / h_avg;
        } else {
            result->y_smooth[i] = y[i];
            result->y_deriv[i] = 0.0;
        }
        
        free(c_bound_func);
        free(c_bound_deriv);
    }
    
    /* Clean up pre-computed coefficients */
    free(c_func);
    free(c_deriv);
    
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
