/* Tikhonov regularization for data smoothing
 * Second derivative penalty via (D²)ᵀ W D² (pentadiagonal Gram matrix)
 * V5.0/2026-02-07/ Corrected penalty from 1st to 2nd order: (D²)ᵀWD² pentadiagonal matrix
 * V4.7/2025-11-28/ Fixed boundary condition asymmetry in Local Spacing Method
 * V4.6/2025-11-28/ Fixed critical bugs: discretization consistency, functional computation
 * V4.5/2025-11-21/ Refactored memory management (goto pattern)
 * V4.4/2025-10-13/ Fixed discretization for non-uniform grids
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tikhonov.h"
#include "grid_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* LAPACK function declarations */
extern void dpbsv_(char *uplo, int *n, int *kd, int *nrhs, 
                   double *ab, int *ldab, double *b, int *ldb, int *info);

/* BLAS function declarations */
extern void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

/* CV threshold for selecting discretization method */
#define CV_THRESHOLD 0.15

/* ============================================================================
 * FIX 1.1: Unified discretization method selection
 * 
 * This enum and helper function ensure consistent method selection
 * across build_band_matrix() and compute_functional().
 * ============================================================================
 */
typedef enum {
    DISCRETIZATION_AVERAGE,  /* Average coefficient method - for near-uniform grids */
    DISCRETIZATION_LOCAL     /* Local spacing method - for non-uniform grids */
} DiscretizationMethod;

/* Select discretization method based on grid properties
 * 
 * Parameters:
 *   grid_info - Grid analysis results (can be NULL)
 *   x         - Array of x-coordinates (used as fallback if grid_info is NULL)
 *   n         - Number of points
 * 
 * Returns:
 *   DISCRETIZATION_AVERAGE for near-uniform grids (CV < 0.15)
 *   DISCRETIZATION_LOCAL for non-uniform grids
 */
static DiscretizationMethod select_discretization_method(GridAnalysis *grid_info, 
                                                          const double *x, int n)
{
    double cv;
    
    if (grid_info != NULL) {
        cv = grid_info->cv;
    } else {
        /* Fallback: compute CV from x array */
        if (n < 2) {
            return DISCRETIZATION_AVERAGE;
        }
        
        double h_sum = 0.0;
        double h_sq_sum = 0.0;
        
        for (int i = 1; i < n; i++) {
            double h = x[i] - x[i-1];
            h_sum += h;
            h_sq_sum += h * h;
        }
        
        double h_avg = h_sum / (n - 1);
        double h_var = h_sq_sum / (n - 1) - h_avg * h_avg;
        
        if (h_var < 0.0) h_var = 0.0;  /* Numerical safety */
        double h_std = sqrt(h_var);
        
        cv = (h_avg > 1e-15) ? (h_std / h_avg) : 0.0;
    }
    
    return (cv < CV_THRESHOLD) ? DISCRETIZATION_AVERAGE : DISCRETIZATION_LOCAL;
}

/* Build band matrix with hybrid discretization */
static void build_band_matrix(double *x, int n, double lambda, double *AB, int ldab, int kd,
                              GridAnalysis *grid_info)
{
    int j;

    /* Clear matrix memory explicitly just to be safe, 
     * mostly handled by calloc in caller but useful if reused */
    memset(AB, 0, ldab * n * sizeof(double));

    /* Identity matrix (Data fidelity term) */
    for (j = 0; j < n; j++) {
        /* Diagonal element is at row index 'kd' in LAPACK band storage */
        AB[kd + j*ldab] = 1.0;
    }

    if (lambda <= 0.0 || n < 3) return;

    /* FIX 1.1: Use unified method selection */
    DiscretizationMethod method = select_discretization_method(grid_info, x, n);

    if (method == DISCRETIZATION_AVERAGE) {
        /* Average Coefficient Method: (D²)ᵀ D² with uniform stencil
         * D² is (n-2)×n matrix, row k: [1, -2, 1] / h_avg²
         * (D²)ᵀ D² = sum_k d_k^T d_k / h_avg^4
         * Resulting pentadiagonal matrix with stencil [1, -4, 6, -4, 1]/h^4
         * Natural BCs (D²u = 0 at boundaries) are implicit — no boundary rows in D²
         */
        double h_avg = (x[n-1] - x[0]) / (n - 1);
        double h4 = h_avg * h_avg * h_avg * h_avg;
        double coeff = lambda / h4;

        for (int k = 1; k <= n-2; k++) {
            /* d_k = [1, -2, 1] at positions k-1, k, k+1 */
            /* Accumulate rank-1 contribution d_k^T d_k (upper triangle only) */

            /* Diagonal: a²=1, b²=4, c²=1 */
            AB[kd + (k-1)*ldab] += coeff * 1.0;
            AB[kd + k*ldab]     += coeff * 4.0;
            AB[kd + (k+1)*ldab] += coeff * 1.0;

            /* 1st superdiagonal: a*b = -2, b*c = -2 */
            AB[(kd-1) + k*ldab]     += coeff * (-2.0);
            AB[(kd-1) + (k+1)*ldab] += coeff * (-2.0);

            /* 2nd superdiagonal: a*c = 1 */
            AB[(kd-2) + (k+1)*ldab] += coeff * 1.0;
        }

    } else {
        /* Local Spacing Method: (D²)ᵀ W D² for non-uniform grid
         * Row k of D²: (D²u)_k = (2/(h_l+h_r)) * [u_{k-1}/h_l - u_k*(1/h_l+1/h_r) + u_{k+1}/h_r]
         * Integration weight: w_k = (h_l + h_r) / 2
         * Matrix = sum_k w_k * d_k^T d_k — automatically symmetric (Gram matrix)
         * Natural BCs implicit — no boundary rows in D²
         */
        for (int k = 1; k <= n-2; k++) {
            double h_l = x[k] - x[k-1];
            double h_r = x[k+1] - x[k];
            double h_sum = h_l + h_r;
            double w_k = h_sum / 2.0;

            double a = 2.0 / (h_sum * h_l);        /* coeff of u_{k-1} */
            double b = -2.0 / (h_l * h_r);         /* coeff of u_k */
            double c = 2.0 / (h_sum * h_r);        /* coeff of u_{k+1} */

            double lw = lambda * w_k;

            /* Accumulate d_k^T * w_k * d_k (upper triangle only) */
            /* Diagonal */
            AB[kd + (k-1)*ldab] += lw * a * a;
            AB[kd + k*ldab]     += lw * b * b;
            AB[kd + (k+1)*ldab] += lw * c * c;

            /* 1st superdiagonal */
            AB[(kd-1) + k*ldab]     += lw * a * b;
            AB[(kd-1) + (k+1)*ldab] += lw * b * c;

            /* 2nd superdiagonal */
            AB[(kd-2) + (k+1)*ldab] += lw * a * c;
        }
    }
}

static void compute_derivatives(double *x, double *y_smooth, int n, double *y_deriv)
{
    if (n < 2) {
        if (n == 1) y_deriv[0] = 0.0;
        return;
    }
    
    if (n == 2) {
        double slope = (y_smooth[1] - y_smooth[0]) / (x[1] - x[0]);
        y_deriv[0] = slope;
        y_deriv[1] = slope;
        return;
    }
    
    /* Forward difference at start */
    y_deriv[0] = (y_smooth[1] - y_smooth[0]) / (x[1] - x[0]);
    
    /* Central differences in interior */
    for (int i = 1; i < n-1; i++) {
        y_deriv[i] = (y_smooth[i+1] - y_smooth[i-1]) / (x[i+1] - x[i-1]);
    }
    
    /* Backward difference at end */
    y_deriv[n-1] = (y_smooth[n-1] - y_smooth[n-2]) / (x[n-1] - x[n-2]);
}

/* ============================================================================
 * FIX 1.2 & 1.3: Corrected functional computation
 * 
 * Changes:
 * - Uses unified select_discretization_method() for consistency with build_band_matrix
 * - Always computes reg_term when lambda > 0 and n >= 3 (no dead code paths)
 * - Fixed boundary second derivative computation using proper one-sided formulas
 * ============================================================================
 */
static void compute_functional(double *x, double *y, double *y_smooth, int n, double lambda,
                              GridAnalysis *grid_info,
                              double *data_term, double *reg_term, double *total_functional)
{
    /* Data fidelity term: ||y - u||² */
    *data_term = 0.0;
    for (int i = 0; i < n; i++) {
        double residual = y[i] - y_smooth[i];
        *data_term += residual * residual;
    }
    
    /* Regularization term: λ||D²u||² */
    *reg_term = 0.0;
    
    if (lambda > 0.0 && n >= 3) {
        /* FIX 1.1: Use unified method selection - same as build_band_matrix */
        DiscretizationMethod method = select_discretization_method(grid_info, x, n);
        
        if (method == DISCRETIZATION_AVERAGE) {
            /* Average coefficient method: ||D²u||² with uniform stencil
             * Interior points only — natural BCs (D²u = 0 at boundaries) are implicit */
            double h_avg = (x[n-1] - x[0]) / (n - 1);
            double h_avg_sq = h_avg * h_avg;
            for (int i = 1; i < n-1; i++) {
                double d2u = (y_smooth[i-1] - 2.0*y_smooth[i] + y_smooth[i+1]) / h_avg_sq;
                *reg_term += d2u * d2u;
            }

        } else {
            /* Local spacing method: ||D²u||²_W with local weights
             * Interior points only — natural BCs implicit */
            for (int i = 1; i < n-1; i++) {
                double h_left = x[i] - x[i-1];
                double h_right = x[i+1] - x[i];
                double h_sum = h_left + h_right;

                /* Second derivative with non-uniform spacing */
                double d2u = (2.0 / h_sum) * (
                    y_smooth[i-1] / h_left -
                    y_smooth[i] * (1.0/h_left + 1.0/h_right) +
                    y_smooth[i+1] / h_right
                );
                /* Weight by local interval length for integration */
                *reg_term += d2u * d2u * h_sum / 2.0;
            }
        }
        
        *reg_term *= lambda;
    }
    
    *total_functional = *data_term + *reg_term;
}

/* Main function with IMPROVED Memory Management */
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda,
                                GridAnalysis *grid_info)
{
    /* Initialize pointers to NULL for safe cleanup */
    TikhonovResult *result = NULL;
    double *AB = NULL;
    double *b = NULL;

    /* Input Validation */
    if (x == NULL || y == NULL || n < 1 || lambda < 0) {
        fprintf(stderr, "Error: Invalid input parameters\n");
        return NULL;
    }

    /* Sanity check on x monotonicity */
    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "Error: x array must be strictly increasing\n");
            return NULL;
        }
    }

    /* --- ALLOCATION --- */

    /* 1. Result Structure */
    result = (TikhonovResult *)malloc(sizeof(TikhonovResult));
    if (!result) {
        fprintf(stderr, "Error: Memory allocation failed (struct)\n");
        goto error;
    }
    result->y_smooth = NULL;
    result->y_deriv = NULL;

    /* 2. Data Arrays */
    result->n = n;
    result->lambda = lambda;
    result->y_smooth = (double *)malloc(n * sizeof(double));
    result->y_deriv = (double *)malloc(n * sizeof(double));

    if (!result->y_smooth || !result->y_deriv) {
        fprintf(stderr, "Error: Memory allocation failed (arrays)\n");
        goto error;
    }

    /* 3. Solver Arrays */
    /* Band matrix storage: 
     * LDA = KD + 1 = 2 (1 superdiagonal + 1 diagonal) 
     * Dimensions: ldab * n 
     */
    int kd = 2;
    int ldab = kd + 1;
    
    AB = (double *)calloc(ldab * n, sizeof(double));
    b = (double *)malloc(n * sizeof(double));

    if (!AB || !b) {
        fprintf(stderr, "Error: Memory allocation failed (solver buffers)\n");
        goto error;
    }

    /* --- CALCULATION --- */

    /* Prepare RHS vector (copy y to b) */
    int inc = 1;
    dcopy_(&n, y, &inc, b, &inc);

    /* Build System Matrix */
    build_band_matrix(x, n, lambda, AB, ldab, kd, grid_info);
    
    /* Solve using LAPACK dpbsv */
    char uplo = 'U'; /* Upper triangle of symmetric band matrix */
    int nrhs = 1;
    int info;
    
    dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
    
    if (info != 0) {
        fprintf(stderr, "Error: LAPACK dpbsv failed (info=%d)\n", info);
        if (info > 0) {
            fprintf(stderr, "Matrix not positive definite. Try larger lambda.\n");
        }
        goto error;
    }
    
    /* Copy result back to struct */
    dcopy_(&n, b, &inc, result->y_smooth, &inc);
    
    /* Post-processing */
    compute_derivatives(x, result->y_smooth, n, result->y_deriv);
    
    /* FIX 1.1: Pass grid_info to compute_functional for consistent method selection */
    compute_functional(x, y, result->y_smooth, n, lambda, grid_info,
                      &result->data_term, &result->regularization_term, 
                      &result->total_functional);
    
    /* Cleanup temporary buffers (Success path) */
    free(AB);
    free(b);

    return result;

    /* --- ERROR HANDLER --- */
error:
    if (AB) free(AB);
    if (b) free(b);
    free_tikhonov_result(result); /* Handles internal NULLs safely */
    return NULL;
}

/* Improved GCV with trace(H) penalty to avoid over-fitting */
static double compute_gcv_score_robust(double *x, double *y, int n, double lambda, 
                                       GridAnalysis *grid_info, int verbose)
{
    TikhonovResult *result;
    double rss = 0.0;
    double trace_H;
    double gcv_score;
    
    /* FIX: Pass grid_info for consistent discretization */
    result = tikhonov_smooth(x, y, n, lambda, grid_info);
    
    if (result == NULL) {
        return 1e20;
    }
    
    for (int i = 0; i < n; i++) {
        double residual = y[i] - result->y_smooth[i];
        rss += residual * residual;
    }
    
    /* Check grid uniformity for trace approximation validity */
    double h_min = x[1] - x[0];
    double h_max = h_min;
    for (int i = 1; i < n; i++) {
        double h = x[i] - x[i-1];
        if (h < h_min) h_min = h;
        if (h > h_max) h_max = h;
    }
    double ratio = h_max / h_min;
    
    double h_avg = (x[n-1] - x[0]) / (n-1);
    double h4 = h_avg * h_avg * h_avg * h_avg;
    double scale = lambda / h4;

    if (n <= 5000) {
        /* Analytical trace - NOTE: approximate for non-uniform grids
         * Eigenvalues of (D²)ᵀD² = (eigenvalues of D¹ᵀD¹)²
         * Null space of D² is 2-dimensional (constants + linear), so trace starts at 2.0 */
        trace_H = 2.0;
        for (int k = 1; k <= n-2; k++) {
            double theta = M_PI * k / n;
            double sin_half = sin(theta / 2.0);
            double ev1 = 4.0 * sin_half * sin_half / (h_avg * h_avg);
            double eigenval = ev1 * ev1;
            trace_H += 1.0 / (1.0 + lambda * eigenval);
        }
        
        if (ratio > 2.0 && verbose) {
            printf("# Note: Trace(H) approximation less accurate for non-uniform grid (ratio=%.2f)\n", ratio);
        }
    } else {
        /* Fast approximation for large datasets */
        trace_H = (double)n / (1.0 + sqrt(scale));
        if (trace_H < 2.0) trace_H = 2.0;
        if (trace_H > n-2) trace_H = n-2;
    }
    
    /* Standard GCV */
    double denom = 1.0 - trace_H / n;
    if (denom > 1e-8) {
        gcv_score = (rss / n) / (denom * denom);
    } else {
        gcv_score = 1e20;
    }
    
    /* Add penalty for excessive trace(H) - penalize over-fitting */
    double trace_ratio = trace_H / n;
    if (trace_ratio > 0.7) {
        double penalty = exp(10.0 * (trace_ratio - 0.7));
        gcv_score *= penalty;
    }
    
    if (verbose) {
        printf("# λ=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f (%.2f), GCV=%9.3e\n",
               lambda, result->total_functional, rss, trace_H, trace_ratio, gcv_score);
    }
    
    free_tikhonov_result(result);
    return gcv_score;
}

/* L-curve method: find corner of L-curve (RSS vs Regularization) */
static double find_lambda_lcurve(double *x, double *y, int n, double *lambda_range, int n_lambda,
                                 GridAnalysis *grid_info)
{
    /* Safe allocation with goto cleanup */
    double *rss_vals = NULL;
    double *reg_vals = NULL;
    double *curv_vals = NULL;
    double best_lambda = 0.01;

    rss_vals = (double *)malloc(n_lambda * sizeof(double));
    reg_vals = (double *)malloc(n_lambda * sizeof(double));
    curv_vals = (double *)malloc(n_lambda * sizeof(double));
    
    if (!rss_vals || !reg_vals || !curv_vals) {
        goto lcurve_cleanup;
    }
    
    /* Compute L-curve points */
    for (int i = 0; i < n_lambda; i++) {
        /* FIX: Pass grid_info for consistent discretization */
        TikhonovResult *result = tikhonov_smooth(x, y, n, lambda_range[i], grid_info);
        if (result) {
            /* FIX: Guard against log(0) - use small epsilon */
            double dt = (result->data_term > 1e-300) ? result->data_term : 1e-300;
            double rt = (result->regularization_term > 1e-300) ? result->regularization_term : 1e-300;
            rss_vals[i] = log(dt);
            reg_vals[i] = log(rt);
            free_tikhonov_result(result);
        } else {
            /* Mark invalid points */
            rss_vals[i] = 0.0;
            reg_vals[i] = 0.0;
        }
    }
    
    /* Compute curvature */
    double max_curv = -1e20;
    int best_idx = n_lambda / 2;
    
    for (int i = 1; i < n_lambda - 1; i++) {
        /* Skip invalid points */
        if (rss_vals[i-1] == 0.0 || rss_vals[i] == 0.0 || rss_vals[i+1] == 0.0) {
            continue;
        }
        
        double dx1 = rss_vals[i] - rss_vals[i-1];
        double dx2 = rss_vals[i+1] - rss_vals[i];
        double dy1 = reg_vals[i] - reg_vals[i-1];
        double dy2 = reg_vals[i+1] - reg_vals[i];
        
        double dx = (dx1 + dx2) / 2.0;
        double dy = (dy1 + dy2) / 2.0;
        double ddx = dx2 - dx1;
        double ddy = dy2 - dy1;
        
        double numer = fabs(dx * ddy - dy * ddx);
        double denom = pow(dx*dx + dy*dy, 1.5);
        
        if (denom > 1e-10) {
            curv_vals[i] = numer / denom;
            if (curv_vals[i] > max_curv) {
                max_curv = curv_vals[i];
                best_idx = i;
            }
        }
    }
    
    best_lambda = lambda_range[best_idx];

lcurve_cleanup:
    if (rss_vals) free(rss_vals);
    if (reg_vals) free(reg_vals);
    if (curv_vals) free(curv_vals);
    
    return best_lambda;
}

/* Enhanced lambda selection with multiple methods */
double find_optimal_lambda_gcv(double *x, double *y, int n, GridAnalysis *grid_info)
{
    double best_lambda = 0.01;
    double best_gcv = 1e20;

    if (n < 3) {
        fprintf(stderr, "Warning: Too few points for GCV (n=%d)\n", n);
        return best_lambda;
    }

    printf("# GCV optimization for n=%d points\n", n);
    if (grid_info) {
        printf("# Grid CV = %.3f, using %s method\n", 
               grid_info->cv,
               (grid_info->cv < CV_THRESHOLD) ? "AVERAGE" : "LOCAL");
        if (grid_info->cv > 0.2) {
            printf("# WARNING: Highly non-uniform grid detected. Trace approximation less accurate.\n");
        }
    }
    
    if (n > 20000) {
        /* Large dataset: conservative range + robust GCV */
        double lambdas[] = {1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0};
        int n_lambdas = sizeof(lambdas) / sizeof(lambdas[0]);
        
        for (int i = 0; i < n_lambdas; i++) {
            double gcv = compute_gcv_score_robust(x, y, n, lambdas[i], grid_info, 1);
            if (gcv < best_gcv) {
                best_gcv = gcv;
                best_lambda = lambdas[i];
            }
        }
        
        double lambda_lcurve = find_lambda_lcurve(x, y, n, lambdas, n_lambdas, grid_info);
        printf("# L-curve suggests λ = %.6e\n", lambda_lcurve);
        
        if (fabs(log10(lambda_lcurve) - log10(best_lambda)) > 0.5) {
            printf("# GCV and L-curve disagree - using more conservative choice\n");
            best_lambda = (lambda_lcurve > best_lambda) ? lambda_lcurve : best_lambda;
        }
        
    } else {
        /* Standard GCV search */
        double lambda_min = 1e-8;
        double lambda_max = 1e0;
        int n_points = 13;
        
        for (int i = 0; i < n_points; i++) {
            double log_lambda = log10(lambda_min) + (log10(lambda_max) - log10(lambda_min)) * i / (n_points - 1);
            double lambda_test = pow(10.0, log_lambda);
            
            double gcv = compute_gcv_score_robust(x, y, n, lambda_test, grid_info, 1);
            
            if (gcv < best_gcv) {
                best_gcv = gcv;
                best_lambda = lambda_test;
            }
        }
        
        /* Refinement around best_lambda */
        if (n <= 5000) {
            printf("# Refinement around λ=%.6e\n", best_lambda);
            for (int i = 0; i < 8; i++) {
                double factor = 0.3 + 1.4 * i / 7.0;
                double lambda_test = best_lambda * factor;
                
                if (lambda_test > lambda_min && lambda_test < lambda_max) {
                    double gcv = compute_gcv_score_robust(x, y, n, lambda_test, grid_info, 1);
                    if (gcv < best_gcv) {
                        best_gcv = gcv;
                        best_lambda = lambda_test;
                    }
                }
            }
        }
    }
    
    printf("# Optimal λ: %.6e (GCV=%.3e)\n", best_lambda, best_gcv);
    return best_lambda;
}

void free_tikhonov_result(TikhonovResult *result)
{
    if (result != NULL) {
        if (result->y_smooth) free(result->y_smooth);
        if (result->y_deriv) free(result->y_deriv);
        free(result);
    }
}
