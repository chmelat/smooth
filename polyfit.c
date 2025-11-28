/* Polynomial fitting for data smoothing
 * Implementation of least squares polynomial approximation
 * V3.1/2025-11-28/ SVD solver instead of QR for maximum robustness
 * V3.0/2025-11-28/ Major: QR decomposition instead of normal equations, condition check
 * V2.3/2025-11-28/ Removed unused coeffs field
 * V2.2/2025-11-28/ Fixes: boundary fallback on LAPACK failure, numerically stable derivative
 * V2.1/2025-11-23/ Fixes: sx array size, performance opt, robust error handling
 * V2.0/2025-05-28/ Updated to use LAPACK/BLAS instead of lib_matrix
 * V1.0/2025-05-27/ Extracted from smooth.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "polyfit.h"

#ifndef DPMAX  
#define DPMAX 12  /* Maximum degree of approximation polynomial */
#endif

/* SVD truncation threshold: singular values below rcond * s_max are treated as zero.
 * This provides implicit regularization for ill-conditioned systems.
 * 
 * Value 1e-10 means: ignore directions with singular values < 1e-10 * largest
 * This is conservative - machine epsilon is ~1e-16 for double precision.
 */
#define SVD_RCOND 1e-10

/* Condition number threshold for warning message */
#define CONDITION_WARNING_THRESHOLD 1e8

/* LAPACK function declarations */

/* SVD-based least squares solver
 * Solves min||b - A*x||_2 using Singular Value Decomposition
 * Most stable method for ill-conditioned problems
 * 
 * Singular values s[i] < rcond * s[0] are treated as zero,
 * providing automatic regularization.
 */
extern void dgelss_(int *m, int *n, int *nrhs,
                    double *a, int *lda, double *b, int *ldb,
                    double *s, double *rcond, int *rank,
                    double *work, int *lwork, int *info);

/* ============================================================================
 * Helper function: Evaluate polynomial and its derivative at a point
 * 
 * Uses Horner's scheme for numerical stability.
 * ============================================================================
 */
static void evaluate_polynomial(const double *coeffs, int degree, double dx,
                                double *value, double *derivative)
{
    double f = coeffs[degree];
    double df = 0.0;
    
    for (int m = degree - 1; m >= 0; m--) {
        df = df * dx + (m + 1) * coeffs[m + 1];
        f = f * dx + coeffs[m];
    }
    
    *value = f;
    *derivative = (degree > 0) ? df : 0.0;
}

/* ============================================================================
 * Helper function: Build Vandermonde matrix for least squares
 * 
 * V[i,j] = (x[i] - x_center)^j
 * 
 * Matrix is stored column-major for LAPACK.
 * ============================================================================
 */
static void build_vandermonde(const double *x, int x_start, int x_end, 
                              double x_center, int degree, 
                              double *V, int ldv)
{
    int nrows = x_end - x_start + 1;
    int ncols = degree + 1;
    
    for (int i = 0; i < nrows; i++) {
        double dx = x[x_start + i] - x_center;
        double p_dx = 1.0;
        
        for (int j = 0; j < ncols; j++) {
            V[i + j * ldv] = p_dx;  /* Column-major storage */
            p_dx *= dx;
        }
    }
}

/* ============================================================================
 * Helper function: Fill boundary points with fallback values
 * ============================================================================
 */
static void fill_boundary_fallback_left(double *x, double *y, int n, int offset,
                                        double *y_smooth, double *y_deriv)
{
    int i0 = offset;
    int i1 = offset + 1;
    
    if (i1 < n - offset) {
        double slope = (y_smooth[i1] - y_smooth[i0]) / (x[i1] - x[i0]);
        
        for (int k = 0; k < offset; k++) {
            y_smooth[k] = y_smooth[i0] + slope * (x[k] - x[i0]);
            y_deriv[k] = slope;
        }
    } else {
        for (int k = 0; k < offset; k++) {
            y_smooth[k] = y[k];
            y_deriv[k] = 0.0;
        }
    }
}

static void fill_boundary_fallback_right(double *x, double *y, int n, int offset,
                                         double *y_smooth, double *y_deriv)
{
    int i0 = n - offset - 1;
    int i1 = n - offset - 2;
    
    if (i1 >= offset) {
        double slope = (y_smooth[i0] - y_smooth[i1]) / (x[i0] - x[i1]);
        
        for (int k = 0; k < offset; k++) {
            int idx = n - offset + k;
            y_smooth[idx] = y_smooth[i0] + slope * (x[idx] - x[i0]);
            y_deriv[idx] = slope;
        }
    } else {
        for (int k = 0; k < offset; k++) {
            int idx = n - offset + k;
            y_smooth[idx] = y[idx];
            y_deriv[idx] = 0.0;
        }
    }
}

/* ============================================================================
 * Main polynomial fitting function
 * 
 * V3.1 changes:
 * - Uses SVD decomposition (dgelss) for maximum numerical stability
 * - Automatic regularization via rcond parameter
 * - Reports effective rank and condition number from singular values
 * ============================================================================
 */
PolyfitResult* polyfit_smooth(double *x, double *y, int n, int window_size, int poly_degree)
{
    PolyfitResult *result = NULL;
    double *V = NULL;           /* Vandermonde matrix */
    double *rhs = NULL;         /* Right-hand side / solution vector */
    double *work = NULL;        /* LAPACK workspace */
    double *sing_vals = NULL;   /* Singular values */
    
    int i, j, k;
    int offset;
    int info;
    int lwork;
    int effective_rank;
    double work_query;
    double rcond = SVD_RCOND;
    
    int left_boundary_done = 0;
    int right_boundary_done = 0;
    int condition_warning_printed = 0;
    int rank_warning_printed = 0;
    
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
    
    offset = window_size / 2;
    int matrix_cols = poly_degree + 1;
    
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
    
    if (result->y_smooth == NULL || result->y_deriv == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for result arrays\n");
        free_polyfit_result(result);
        return NULL;
    }
    
    /* Allocate Vandermonde matrix, RHS vector, and singular values
     * rhs needs to be at least max(window_size, matrix_cols) for dgelss
     */
    int rhs_size = (window_size > matrix_cols) ? window_size : matrix_cols;
    
    V = (double *)malloc(window_size * matrix_cols * sizeof(double));
    rhs = (double *)malloc(rhs_size * sizeof(double));
    sing_vals = (double *)malloc(matrix_cols * sizeof(double));
    
    if (V == NULL || rhs == NULL || sing_vals == NULL) {
        fprintf(stderr, "Error: Matrix allocation failed\n");
        goto cleanup_error;
    }
    
    /* Query optimal workspace size for dgelss */
    int nrhs = 1;
    lwork = -1;
    dgelss_(&window_size, &matrix_cols, &nrhs, V, &window_size, rhs, &rhs_size,
            sing_vals, &rcond, &effective_rank, &work_query, &lwork, &info);
    
    lwork = (int)work_query + 1;
    work = (double *)malloc(lwork * sizeof(double));
    
    if (work == NULL) {
        fprintf(stderr, "Error: Workspace allocation failed\n");
        goto cleanup_error;
    }
    
    /* Process each point in the interior */
    for (i = offset; i < n - offset; i++) {
        
        /* Build Vandermonde matrix */
        build_vandermonde(x, i - offset, i + offset, x[i], poly_degree, V, window_size);
        
        /* Copy y values to RHS vector */
        for (j = 0; j < window_size; j++) {
            rhs[j] = y[i - offset + j];
        }
        
        /* Solve least squares using SVD decomposition */
        dgelss_(&window_size, &matrix_cols, &nrhs, V, &window_size, 
                rhs, &rhs_size, sing_vals, &rcond, &effective_rank,
                work, &lwork, &info);
        
        if (info != 0) {
            /* SVD failed - fallback to original value */
            result->y_smooth[i] = y[i];
            result->y_deriv[i] = 0.0;
            continue;
        }
        
        /* Check condition number and rank on first window */
        if (i == offset) {
            /* Condition number = s_max / s_min (for non-zero singular values) */
            double s_max = sing_vals[0];
            double s_min = sing_vals[effective_rank - 1];
            double cond = (s_min > 0) ? (s_max / s_min) : 1e16;
            
            if (cond > CONDITION_WARNING_THRESHOLD && !condition_warning_printed) {
                fprintf(stderr, "Note: Vandermonde matrix condition number = %.2e\n", cond);
                condition_warning_printed = 1;
            }
            
            if (effective_rank < matrix_cols && !rank_warning_printed) {
                fprintf(stderr, "Note: Matrix is rank-deficient (rank=%d of %d). "
                        "SVD regularization active.\n", effective_rank, matrix_cols);
                rank_warning_printed = 1;
            }
        }
        
        /* Store smoothed value and derivative at center point
         * rhs[0] = c0 (value at center), rhs[1] = c1 (first derivative)
         */
        result->y_smooth[i] = rhs[0];
        result->y_deriv[i] = (poly_degree > 0) ? rhs[1] : 0.0;
        
        /* Handle left boundary on first interior point */
        if (i == offset) {
            for (k = 0; k < offset; k++) {
                double dx = x[k] - x[offset];
                double fi, dfi;
                evaluate_polynomial(rhs, poly_degree, dx, &fi, &dfi);
                result->y_smooth[k] = fi;
                result->y_deriv[k] = dfi;
            }
            left_boundary_done = 1;
        }
        
        /* Handle right boundary on last interior point */
        if (i == n - offset - 1) {
            for (k = 0; k < offset; k++) {
                double dx = x[n - offset + k] - x[n - offset - 1];
                double fi, dfi;
                evaluate_polynomial(rhs, poly_degree, dx, &fi, &dfi);
                result->y_smooth[n - offset + k] = fi;
                result->y_deriv[n - offset + k] = dfi;
            }
            right_boundary_done = 1;
        }
    }
    
    /* Handle case where boundaries were not processed */
    if (!left_boundary_done && offset > 0) {
        fprintf(stderr, "Warning: Left boundary filled with fallback\n");
        fill_boundary_fallback_left(x, y, n, offset, result->y_smooth, result->y_deriv);
    }
    
    if (!right_boundary_done && offset > 0) {
        fprintf(stderr, "Warning: Right boundary filled with fallback\n");
        fill_boundary_fallback_right(x, y, n, offset, result->y_smooth, result->y_deriv);
    }
    
    /* Successful cleanup */
    free(V);
    free(rhs);
    free(work);
    free(sing_vals);
    
    return result;

cleanup_error:
    free(V);
    free(rhs);
    free(work);
    free(sing_vals);
    free_polyfit_result(result);
    return NULL;
}

/* Free allocated memory */
void free_polyfit_result(PolyfitResult *result)
{
    if (result != NULL) {
        free(result->y_smooth);
        free(result->y_deriv);
        free(result);
    }
}
