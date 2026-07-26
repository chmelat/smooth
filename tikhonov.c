/* Tikhonov regularization for data smoothing
 * Second derivative penalty via (D²)ᵀ W D² (pentadiagonal Gram matrix)
 * V5.4/2026-07-26/ Removed the L-curve method and the n>20000 branch it was
 *                  reachable from; one log-spaced GCV sweep for all n
 * V5.3/2026-06-10/ GCV: analytical trace for all n (removed n>5000 shortcut with
 *                  wrong asymptotic exponent); warn when optimal lambda is pinned
 *                  to the search-range edge
 * V5.2/2026-06-01/ L-curve: use explicit valid[] flag instead of a 0.0 sentinel
 *                  (log(data_term) can legitimately be 0.0)
 * V5.1/2026-05-30/ A1: single integral-measure discretization (removed AVERAGE/LOCAL branch + CV=0.15 switch)
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
extern void dcopy_(int *n, const double *x, int *incx, double *y, int *incy);

/* Build band matrix: A = I + lambda * (D2)^T W D2 (pentadiagonal Gram matrix) */
static void build_band_matrix(const double *x, int n, double lambda, double *AB, int ldab, int kd)
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

    /* Penalty discretizes the integral lambda * integral (u'')^2 dx via the
     * Gram matrix of the second-derivative operator D2:
     *   Matrix = sum_k w_k * d_k^T d_k   (symmetric, pentadiagonal, kd=2)
     * Row k of D2 (interior point k, 1 <= k <= n-2):
     *   (D2 u)_k = (2/(h_l+h_r)) * [u_{k-1}/h_l - u_k*(1/h_l+1/h_r) + u_{k+1}/h_r]
     * Integration weight: w_k = (h_l + h_r)/2.
     * On a uniform grid (h_l = h_r = h) this reduces to the [1,-4,6,-4,1]*lambda/h^3
     * stencil. Natural BCs (D2 u = 0 at the endpoints) are implicit: D2 has no
     * rows for boundary points, so the grid spacing carries through consistently
     * for uniform and non-uniform grids alike. */
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

static void compute_derivatives(const double *x, const double *y_smooth, int n, double *y_deriv)
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

/* Compute data fidelity, regularization, and total functional values.
 * Uses the same integral-measure discretization as build_band_matrix. */
static void compute_functional(const double *x, const double *y, const double *y_smooth, int n, double lambda,
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
        /* Regularization term lambda * integral (u'')^2 dx, discretized with the
         * same Gram-matrix weighting as build_band_matrix (interior points only;
         * natural BCs implicit). Uniform and non-uniform grids share one path. */
        for (int i = 1; i < n-1; i++) {
            double h_left = x[i] - x[i-1];
            double h_right = x[i+1] - x[i];
            double h_sum = h_left + h_right;

            /* Second derivative with (possibly non-uniform) spacing */
            double d2u = (2.0 / h_sum) * (
                y_smooth[i-1] / h_left -
                y_smooth[i] * (1.0/h_left + 1.0/h_right) +
                y_smooth[i+1] / h_right
            );
            /* Weight by local interval length for integration */
            *reg_term += d2u * d2u * h_sum / 2.0;
        }

        *reg_term *= lambda;
    }
    
    *total_functional = *data_term + *reg_term;
}

/* Main function with IMPROVED Memory Management */
TikhonovResult* tikhonov_smooth(const double *x, const double *y, int n, double lambda,
                                const GridAnalysis *grid_info)
{
    /* Initialize pointers to NULL for safe cleanup */
    TikhonovResult *result = NULL;
    double *AB = NULL;
    double *b = NULL;

    /* Input Validation */
    if (x == NULL || y == NULL || n < 1 || lambda < 0) {
        fprintf(stderr, "ERROR: Invalid input parameters\n");
        return NULL;
    }

    if (grid_info == NULL) {
        fprintf(stderr, "ERROR: Grid info not available\n");
        return NULL;
    }

    /* Sanity check on x monotonicity */
    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "ERROR: x array must be strictly increasing\n");
            return NULL;
        }
    }

    /* --- ALLOCATION --- */

    /* 1. Result Structure */
    result = (TikhonovResult *)malloc(sizeof(TikhonovResult));
    if (!result) {
        fprintf(stderr, "ERROR: Memory allocation failed (struct)\n");
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
        fprintf(stderr, "ERROR: Memory allocation failed (arrays)\n");
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
        fprintf(stderr, "ERROR: Memory allocation failed (solver buffers)\n");
        goto error;
    }

    /* --- CALCULATION --- */

    /* Prepare RHS vector (copy y to b) */
    int inc = 1;
    dcopy_(&n, y, &inc, b, &inc);

    /* Build System Matrix */
    build_band_matrix(x, n, lambda, AB, ldab, kd);
    
    /* Solve using LAPACK dpbsv */
    char uplo = 'U'; /* Upper triangle of symmetric band matrix */
    int nrhs = 1;
    int info;
    
    dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
    
    if (info != 0) {
        fprintf(stderr, "ERROR: LAPACK dpbsv failed (info=%d)\n", info);
        if (info > 0) {
            fprintf(stderr, "Leading minor of order %d not positive definite. "
                            "Unexpected: I + lambda*(D2)^T D2 is SPD for lambda>=0, "
                            "so this points to numerical ill-conditioning.\n", info);
        }
        goto error;
    }
    
    /* Copy result back to struct */
    dcopy_(&n, b, &inc, result->y_smooth, &inc);
    
    /* Post-processing */
    compute_derivatives(x, result->y_smooth, n, result->y_deriv);
    
    compute_functional(x, y, result->y_smooth, n, lambda,
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

/* Standard Generalized Cross Validation score for a single lambda */
static double compute_gcv_score_robust(const double *x, const double *y, int n, double lambda,
                                       const GridAnalysis *grid_info)
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
    
    /* Grid uniformity for trace-approximation validity (already in grid_info) */
    double ratio = grid_info->ratio_max_min;
    double h_avg = grid_info->h_avg;

    /* Analytical trace from the uniform-grid eigenvalue model of the penalty
     * (approximate for non-uniform grids). Penalty matrix K = sum_k w_k d_k^T d_k
     * uses the integral measure, so on a uniform grid its eigenvalues are
     * h_avg * (4 sin^2(theta/2)/h^2)^2. Null space of D2 is 2-dimensional
     * (constants + linear), so trace starts at 2.0. The O(n) sum is the same
     * order as the band solve, so it is used for every n. */
    trace_H = 2.0;
    for (int k = 1; k <= n-2; k++) {
        double theta = M_PI * k / n;
        double sin_half = sin(theta / 2.0);
        double ev1 = 4.0 * sin_half * sin_half / (h_avg * h_avg);
        double eigenval = ev1 * ev1 * h_avg;
        trace_H += 1.0 / (1.0 + lambda * eigenval);
    }

    if (ratio > 2.0) {
        printf("# Note: Trace(H) approximation less accurate for non-uniform grid (ratio=%.2f)\n", ratio);
    }
    
    /* Standard GCV */
    double denom = 1.0 - trace_H / n;
    if (denom > 1e-8) {
        gcv_score = (rss / n) / (denom * denom);
    } else {
        gcv_score = 1e20;
    }
    
    printf("# λ=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f (%.2f), GCV=%9.3e\n",
           lambda, result->total_functional, rss, trace_H, trace_H / n, gcv_score);


    free_tikhonov_result(result);
    return gcv_score;
}

/* Enhanced lambda selection with multiple methods */
double find_optimal_lambda_gcv(const double *x, const double *y, int n, const GridAnalysis *grid_info)
{
    double best_lambda = 0.01;
    double best_gcv = 1e20;
    /* Search range of the log-spaced GCV sweep */
    const double lambda_min = 1e-8;
    const double lambda_max = 1e0;

    if (grid_info == NULL) {
        fprintf(stderr, "ERROR: Grid info not available\n");
        return best_lambda;
    }

    if (n < 3) {
        fprintf(stderr, "Warning: Too few points for GCV (n=%d)\n", n);
        return best_lambda;
    }

    printf("# GCV optimization for n=%d points (see README: Generalized Cross Validation)\n", n);
    printf("# Grid CV = %.3f (integral-measure penalty discretization)\n",
           grid_info->cv);
    if (grid_info->cv > 0.2) {
        printf("# WARNING: Highly non-uniform grid detected. Trace approximation less accurate.\n");
    }
    
    /* Log-spaced GCV search over [lambda_min, lambda_max] */
    int n_points = 13;

    for (int i = 0; i < n_points; i++) {
        double log_lambda = log10(lambda_min) + (log10(lambda_max) - log10(lambda_min)) * i / (n_points - 1);
        double lambda_test = pow(10.0, log_lambda);

        double gcv = compute_gcv_score_robust(x, y, n, lambda_test, grid_info);

        if (gcv < best_gcv) {
            best_gcv = gcv;
            best_lambda = lambda_test;
        }
    }

    /* Refinement around best_lambda (skipped for large n: each score is a full
     * O(n) solve, so the extra 8 evaluations are not worth the sub-grid gain) */
    if (n <= 5000) {
        printf("# Refinement around λ=%.6e\n", best_lambda);
        for (int i = 0; i < 8; i++) {
            double factor = 0.3 + 1.4 * i / 7.0;
            double lambda_test = best_lambda * factor;

            if (lambda_test > lambda_min && lambda_test < lambda_max) {
                double gcv = compute_gcv_score_robust(x, y, n, lambda_test, grid_info);
                if (gcv < best_gcv) {
                    best_gcv = gcv;
                    best_lambda = lambda_test;
                }
            }
        }
    }


    /* Lambda is dimensional (scales with h^3 and the y amplitude), so a fixed
     * search range cannot fit every data scale. Flag a result pinned to the
     * range edge instead of returning it silently as "optimal". */
    if (best_lambda <= lambda_min * 1.01 || best_lambda >= lambda_max * 0.99) {
        printf("# WARNING: optimal λ = %.3e lies at the edge of the search range [%.0e, %.0e].\n",
               best_lambda, lambda_min, lambda_max);
        printf("#          The true optimum may lie outside this range; consider setting λ manually (-l <value>).\n");
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
