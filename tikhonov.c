/*  Tikhonov regularization for data smoothing
 *  Hybrid version: combines best of old and new approaches
 *  V4.4/2025-10-13/ Fixed discretization for non-uniform grids
 *  FIXES:
 *  - Corrected boundary superdiagonal in local spacing method
 *  - Fixed compute_functional to use proper non-uniform discretization
 *  - Harmonic mean for interval averaging (more accurate than arithmetic)
 *  - Added GCV warning for highly non-uniform grids
 *  - Consistent boundary conditions throughout
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tikhonov.h"
#include "grid_analysis.h"

/* LAPACK function declarations */
extern void dpbsv_(char *uplo, int *n, int *kd, int *nrhs, 
                   double *ab, int *ldab, double *b, int *ldb, int *info);

/* BLAS function declarations */
extern void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

/* CV threshold for selecting discretization method */
#define CV_THRESHOLD 0.15

/* Build band matrix with hybrid discretization */
static void build_band_matrix(double *x, int n, double lambda, double *AB, int ldab, int kd,
                              GridAnalysis *grid_info)
{
    int j;

    memset(AB, 0, ldab * n * sizeof(double));

    /* Identity matrix */
    for (j = 0; j < n; j++) {
        AB[kd + j*ldab] = 1.0;
    }

    if (lambda <= 0.0 || n < 3) return;

    /* Select discretization method based on CV */
    int use_average_coef = (grid_info && grid_info->cv < CV_THRESHOLD);

    if (use_average_coef) {
        /* Average Coefficient Method (from old version)
         * Better for uniform/mildly non-uniform grids */
        
        double h_avg_squared = 0.0;
        for (int i = 1; i < n; i++) {
            double h_local = x[i] - x[i-1];
            h_avg_squared += 1.0 / (h_local * h_local);
        }
        double c = lambda * h_avg_squared / (n-1);
        
        /* Interior points */
        for (j = 1; j < n-1; j++) {
            AB[kd + j*ldab] += c * 2.0;
        }
        
        /* Boundaries */
        AB[kd + 0*ldab] += c * 1.0;
        AB[kd + (n-1)*ldab] += c * 1.0;
        
        /* Superdiagonal */
        for (j = 1; j < n; j++) {
            AB[0 + j*ldab] += c * (-1.0);
        }

        if (grid_info && grid_info->cv > 0.01) {
            printf("# Using Average Coefficient Method (CV = %.3f)\n", grid_info->cv);
        }

    } else {
        /* Local Spacing Method (from new version)
         * Better for highly non-uniform grids */

        if (grid_info) {
            printf("# Using Local Spacing Method (CV = %.3f)\n", grid_info->cv);
        }
        
        /* Interior points */
        for (j = 1; j < n-1; j++) {
            double h_left = x[j] - x[j-1];
            double h_right = x[j+1] - x[j];
            double h_sum = h_left + h_right;
            double w = 2.0 * lambda / h_sum;
            
            AB[kd + j*ldab] += w * (1.0/h_left + 1.0/h_right);
            AB[0 + (j+1)*ldab] += -w / h_right;
        }
        
        /* FIXED: Boundary conditions with consistent discretization */
        double h_first = x[1] - x[0];
        AB[kd + 0*ldab] += lambda / (h_first * h_first);
        AB[0 + 1*ldab] += -lambda / (h_first * h_first);  // FIX: Was missing!
        
        double h_last = x[n-1] - x[n-2];
        AB[kd + (n-1)*ldab] += lambda / (h_last * h_last);
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

/* FIXED: Proper discretization for non-uniform grids */
static void compute_functional(double *x, double *y, double *y_smooth, int n, double lambda,
                              double *data_term, double *reg_term, double *total_functional)
{
    /* Data fidelity term */
    *data_term = 0.0;
    for (int i = 0; i < n; i++) {
        double residual = y[i] - y_smooth[i];
        *data_term += residual * residual;
    }
    
    /* Regularization term with correct non-uniform discretization */
    *reg_term = 0.0;
    
    if (lambda > 0.0 && n >= 3) {
        /* Determine which discretization method was used */
        double h_min = x[1] - x[0];
        double h_max = h_min;
        
        for (int i = 1; i < n; i++) {
            double h = x[i] - x[i-1];
            if (h < h_min) h_min = h;
            if (h > h_max) h_max = h;
        }
        
        double ratio = h_max / h_min;
        
        if (ratio < 2.5) {
            /* Average coefficient method - use simplified formula with harmonic mean */
            for (int i = 1; i < n-1; i++) {
                double h_left = x[i] - x[i-1];
                double h_right = x[i+1] - x[i];
                /* Harmonic mean is more appropriate for averaging intervals */
                double h_harm = 2.0 * h_left * h_right / (h_left + h_right);
                
                double d2u = (y_smooth[i-1] - 2.0*y_smooth[i] + y_smooth[i+1]) / (h_harm * h_harm);
                *reg_term += d2u * d2u;
            }
            
            /* Boundaries with simple discretization */
            double h_start = x[1] - x[0];
            double d2u_left = (y_smooth[1] - y_smooth[0]) / (h_start * h_start);
            *reg_term += 0.5 * d2u_left * d2u_left;
            
            double h_end = x[n-1] - x[n-2];
            double d2u_right = (y_smooth[n-1] - y_smooth[n-2]) / (h_end * h_end);
            *reg_term += 0.5 * d2u_right * d2u_right;
            
        } else {
            /* Local spacing method - use proper non-uniform discretization */
            
            /* Interior points: D²u_i = (2/(h_i+h_{i+1})) * (u_{i-1}/h_i - u_i*(1/h_i + 1/h_{i+1}) + u_{i+1}/h_{i+1}) */
            for (int i = 1; i < n-1; i++) {
                double h_left = x[i] - x[i-1];
                double h_right = x[i+1] - x[i];
                double h_sum = h_left + h_right;
                
                double d2u = (2.0 / h_sum) * (
                    y_smooth[i-1] / h_left - 
                    y_smooth[i] * (1.0/h_left + 1.0/h_right) + 
                    y_smooth[i+1] / h_right
                );
                
                /* Weight by interval length for proper integration */
                *reg_term += d2u * d2u * h_sum / 2.0;
            }
            
            /* Left boundary: D²u_0 = (u_1 - u_0) / h_0² */
            double h_start = x[1] - x[0];
            double d2u_left = (y_smooth[1] - y_smooth[0]) / (h_start * h_start);
            *reg_term += d2u_left * d2u_left * h_start / 2.0;
            
            /* Right boundary: D²u_{n-1} = (u_{n-1} - u_{n-2}) / h_{n-1}² */
            double h_end = x[n-1] - x[n-2];
            double d2u_right = (y_smooth[n-1] - y_smooth[n-2]) / (h_end * h_end);
            *reg_term += d2u_right * d2u_right * h_end / 2.0;
        }
        
        *reg_term *= lambda;
    }
    
    *total_functional = *data_term + *reg_term;
}

TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda,
                                GridAnalysis *grid_info)
{
    TikhonovResult *result;
    double *AB, *b;
    int kd = 1;
    int ldab = kd + 1;
    int nrhs = 1;
    int info;
    int inc = 1;
    char uplo = 'U';

    if (x == NULL || y == NULL || n < 1 || lambda < 0) {
        fprintf(stderr, "Error: Invalid input parameters\n");
        return NULL;
    }

    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "Error: x array must be strictly increasing\n");
            return NULL;
        }
    }

    result = (TikhonovResult *)malloc(sizeof(TikhonovResult));
    if (!result) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return NULL;
    }

    result->n = n;
    result->lambda = lambda;
    result->y_smooth = (double *)malloc(n * sizeof(double));
    result->y_deriv = (double *)malloc(n * sizeof(double));

    if (!result->y_smooth || !result->y_deriv) {
        free_tikhonov_result(result);
        return NULL;
    }

    AB = (double *)calloc(ldab * n, sizeof(double));
    b = (double *)malloc(n * sizeof(double));

    if (!AB || !b) {
        free(AB);
        free(b);
        free_tikhonov_result(result);
        return NULL;
    }

    dcopy_(&n, y, &inc, b, &inc);

    build_band_matrix(x, n, lambda, AB, ldab, kd, grid_info);
    
    dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
    
    if (info != 0) {
        fprintf(stderr, "Error: LAPACK failed (info=%d)\n", info);
        if (info > 0) {
            fprintf(stderr, "Matrix not positive definite. Try larger lambda.\n");
        }
        free(AB); 
        free(b);
        free_tikhonov_result(result);
        return NULL;
    }
    
    dcopy_(&n, b, &inc, result->y_smooth, &inc);
    
    compute_derivatives(x, result->y_smooth, n, result->y_deriv);
    compute_functional(x, y, result->y_smooth, n, lambda,
                      &result->data_term, &result->regularization_term, 
                      &result->total_functional);
    
    free(AB);
    free(b);
    return result;
}

/* Improved GCV with trace(H) penalty to avoid over-fitting */
static double compute_gcv_score_robust(double *x, double *y, int n, double lambda, int verbose)
{
    TikhonovResult *result;
    double rss = 0.0;
    double trace_H;
    double gcv_score;
    
    result = tikhonov_smooth(x, y, n, lambda, NULL);
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
    double scale = lambda / (h_avg * h_avg);
    
    if (n <= 5000) {
        /* Analytical trace - NOTE: approximate for non-uniform grids */
        trace_H = 2.0;
        for (int k = 1; k <= n-2; k++) {
            double theta = M_PI * k / n;
            double eigenval = 4.0 * pow(sin(theta/2.0), 2) / (h_avg * h_avg);
            trace_H += 1.0 / (1.0 + lambda * eigenval);
        }
        
        if (ratio > 2.0 && verbose) {
            printf("# WARNING: trace(H) approximation may be inaccurate for non-uniform grid (ratio=%.1f)\n", ratio);
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
        /* Exponential penalty when trace(H)/n > 0.7 */
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
static double find_lambda_lcurve(double *x, double *y, int n, double *lambda_range, int n_lambda)
{
    double *rss_vals = (double *)malloc(n_lambda * sizeof(double));
    double *reg_vals = (double *)malloc(n_lambda * sizeof(double));
    double *curv_vals = (double *)malloc(n_lambda * sizeof(double));
    
    if (!rss_vals || !reg_vals || !curv_vals) {
        free(rss_vals); free(reg_vals); free(curv_vals);
        return 0.01;
    }
    
    /* Compute L-curve points */
    for (int i = 0; i < n_lambda; i++) {
        TikhonovResult *result = tikhonov_smooth(x, y, n, lambda_range[i], NULL);
        if (result) {
            rss_vals[i] = log(result->data_term);
            reg_vals[i] = log(result->regularization_term);
            free_tikhonov_result(result);
        } else {
            rss_vals[i] = reg_vals[i] = 0.0;
        }
    }
    
    /* Compute curvature at each point */
    double max_curv = -1e20;
    int best_idx = n_lambda / 2;
    
    for (int i = 1; i < n_lambda - 1; i++) {
        /* Numerical curvature: κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2) */
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
    
    double best_lambda = lambda_range[best_idx];
    
    free(rss_vals);
    free(reg_vals);
    free(curv_vals);
    
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

    if (!grid_info) {
        fprintf(stderr, "Warning: Grid info not available\n");
        return best_lambda;
    }

    printf("# GCV optimization for n=%d points (CV = %.3f)\n", n, grid_info->cv);

    if (grid_info->cv > 0.2) {
        printf("# WARNING: Highly non-uniform grid detected (CV = %.3f)!\n", grid_info->cv);
        printf("# GCV trace approximation may be less accurate.\n");
        printf("# Consider using L-curve method or manual lambda selection.\n");
    }
    
    printf("# Format: λ | Functional | RSS | tr(H) (ratio) | GCV\n");
    
    if (n > 20000) {
        /* Large dataset: conservative range + robust GCV */
        printf("# Large dataset: using conservative lambda range with robust GCV\n");
        
        double lambdas[] = {1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0};
        int n_lambdas = sizeof(lambdas) / sizeof(lambdas[0]);
        
        /* Try GCV with penalty */
        for (int i = 0; i < n_lambdas; i++) {
            double gcv = compute_gcv_score_robust(x, y, n, lambdas[i], 1);
            if (gcv < best_gcv) {
                best_gcv = gcv;
                best_lambda = lambdas[i];
            }
        }
        
        /* Also try L-curve as sanity check */
        double lambda_lcurve = find_lambda_lcurve(x, y, n, lambdas, n_lambdas);
        printf("# L-curve suggests λ = %.6e\n", lambda_lcurve);
        
        /* If L-curve and GCV differ significantly, prefer the larger lambda (more conservative) */
        if (fabs(log10(lambda_lcurve) - log10(best_lambda)) > 0.5) {
            printf("# GCV and L-curve disagree - using more conservative choice\n");
            best_lambda = (lambda_lcurve > best_lambda) ? lambda_lcurve : best_lambda;
        }
        
    } else {
        /* Small/medium dataset: full range with standard GCV */
        printf("# Standard GCV search with full lambda range\n");
        
        double lambda_min = 1e-6;
        double lambda_max = 1e0;
        double log_min = log10(lambda_min);
        double log_max = log10(lambda_max);
        int n_points = 13;
        
        for (int i = 0; i < n_points; i++) {
            double log_lambda = log_min + (log_max - log_min) * i / (n_points - 1);
            double lambda_test = pow(10.0, log_lambda);
            
            double gcv = compute_gcv_score_robust(x, y, n, lambda_test, 1);
            
            if (gcv < best_gcv) {
                best_gcv = gcv;
                best_lambda = lambda_test;
            }
        }
        
        /* Refinement */
        if (n <= 5000) {
            printf("# Refinement around λ=%.6e\n", best_lambda);
            
            for (int i = 0; i < 8; i++) {
                double factor = 0.3 + 1.4 * i / 7.0;
                double lambda_test = best_lambda * factor;
                
                if (lambda_test > lambda_min && lambda_test < lambda_max) {
                    double gcv = compute_gcv_score_robust(x, y, n, lambda_test, 1);
                    
                    if (gcv < best_gcv) {
                        best_gcv = gcv;
                        best_lambda = lambda_test;
                    }
                }
            }
        }
    }
    
    printf("# Optimal λ: %.6e (GCV=%.3e)\n", best_lambda, best_gcv);

    /* Final check */
    TikhonovResult *final = tikhonov_smooth(x, y, n, best_lambda, NULL);
    if (final) {
        double reg_ratio = final->regularization_term / final->total_functional;
        printf("# Balance: Data=%.1f%%, Regularization=%.1f%%\n", 
               (1.0-reg_ratio)*100, reg_ratio*100);
        
        if (reg_ratio > 0.95) {
            printf("# WARNING: Over-smoothing detected. Consider manual λ values.\n");
        } else if (reg_ratio < 0.5) {
            printf("# WARNING: Under-smoothing detected. Consider smaller λ.\n");
        }
        
        free_tikhonov_result(final);
    }
    
    return best_lambda;
}

void free_tikhonov_result(TikhonovResult *result)
{
    if (result != NULL) {
        free(result->y_smooth);
        free(result->y_deriv);
        free(result);
    }
}
