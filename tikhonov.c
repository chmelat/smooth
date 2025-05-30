/*  Tikhonov regularization for data smoothing
 *  Production version with efficient band matrix implementation
 *  V3.3/2025-05-27/ Using separate grid_analysis module
 *  V3.2/2025-05-26/ Enhanced with non-uniformity analysis and optional adaptive weights
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

/* Forward declarations for internal functions */
static TikhonovResult* tikhonov_smooth_internal(double *x, double *y, int n, double lambda, int print_grid_info, int adaptive_weights);

/* Build band matrix A = I + lambda*D^T*D using correct LAPACK column-major indexing */
static void build_band_matrix(double *x, int n, double lambda, double *AB, int ldab, int kd, int adaptive_weights)
{
    double h_avg;
    int uniform = is_uniform_grid(x, n, &h_avg, 1e-10);
    
    /* Clear the band matrix */
    memset(AB, 0, ldab * n * sizeof(double));
    
    /* Set main diagonal to identity matrix */
    /* LAPACK column-major: A[j,j] stored at AB[kd + j*ldab] */
    for (int j = 0; j < n; j++) {
        AB[kd + j*ldab] = 1.0;
    }
    
    if (lambda > 0.0 && n >= 3) {
        if (uniform) {
            /* Uniform grid - use constant coefficient */
            double h2 = h_avg * h_avg;
            double c = lambda / h2;
            
            /* Interior points get full second-order penalty */
            for (int j = 1; j < n-1; j++) {
                AB[kd + j*ldab] += c * 2.0;
            }
            
            /* Boundary points get reduced penalty */
            AB[kd + 0*ldab] += c * 1.0;
            AB[kd + (n-1)*ldab] += c * 1.0;
            
            /* Superdiagonal elements */
            for (int j = 1; j < n; j++) {
                AB[0 + j*ldab] += c * (-1.0);
            }
        } else {
            /* Non-uniform grid - check if adaptive weights are requested */
            GridAnalysis *grid_info = analyze_grid(x, n, 0);
            
            if (adaptive_weights) {
                /* User explicitly requested adaptive weights */
                printf("# Using adaptive local weights (user requested)\n");
                printf("# Grid non-uniformity: h_max/h_min = %.2f, CV = %.3f\n", 
                       grid_info->ratio_max_min, grid_info->cv);
                
                /* Use local weights for each point */
                double *interval_weights = (double *)calloc(n-1, sizeof(double));
                if (interval_weights == NULL) {
                    fprintf(stderr, "Memory allocation failed for interval weights\n");
                    /* Fall back to average method */
                    adaptive_weights = 0;
                } else {
                    /* Compute weight for each interval */
                    for (int i = 0; i < n-1; i++) {
                        double h = x[i+1] - x[i];
                        interval_weights[i] = lambda / (h * h);
                    }
                    
                    /* Build matrix maintaining symmetry */
                    /* First point (j=0) */
                    AB[kd + 0*ldab] += interval_weights[0];
                    AB[0 + 1*ldab] += -interval_weights[0];
                    
                    /* Interior points */
                    for (int j = 1; j < n-1; j++) {
                        /* Contribution from left interval */
                        AB[kd + j*ldab] += interval_weights[j-1];
                        /* Contribution from right interval */
                        AB[kd + j*ldab] += interval_weights[j];
                        /* Off-diagonal for next point */
                        AB[0 + (j+1)*ldab] += -interval_weights[j];
                    }
                    
                    /* Last point (j=n-1) */
                    AB[kd + (n-1)*ldab] += interval_weights[n-2];
                    
                    free(interval_weights);
                }
            }
            
            if (!adaptive_weights) {
                /* Use average coefficient method (default for non-uniform grids) */
                printf("# Using average coefficient method for non-uniform grid\n");
                printf("# Grid non-uniformity: h_max/h_min = %.2f, CV = %.3f\n", 
                       grid_info->ratio_max_min, grid_info->cv);
                
                double h_avg_squared = 0.0;
                for (int i = 1; i < n; i++) {
                    double h_local = x[i] - x[i-1];
                    h_avg_squared += 1.0 / (h_local * h_local);
                }
                double c = lambda * h_avg_squared / (n-1);
                
                /* Same pattern as uniform case but with average c */
                for (int j = 1; j < n-1; j++) {
                    AB[kd + j*ldab] += c * 2.0;
                }
                AB[kd + 0*ldab] += c * 1.0;
                AB[kd + (n-1)*ldab] += c * 1.0;
                for (int j = 1; j < n; j++) {
                    AB[0 + j*ldab] += c * (-1.0);
                }
            }
            
            /* Print warnings only if reliability issues exist */
            if (grid_info->reliability_warning) {
                printf("# WARNING: %s\n", grid_info->warning_msg);
                if (!adaptive_weights && grid_info->ratio_max_min > 10.0) {
                    printf("# Consider using adaptive weights (-a flag) for better results.\n");
                }
            }
            
            free_grid_analysis(grid_info);
        }
    }
}

/* Compute derivatives using central differences with boundary handling */
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
    
    /* Forward difference at left boundary */
    y_deriv[0] = (y_smooth[1] - y_smooth[0]) / (x[1] - x[0]);
    
    /* Central differences for interior points */
    for (int i = 1; i < n-1; i++) {
        y_deriv[i] = (y_smooth[i+1] - y_smooth[i-1]) / (x[i+1] - x[i-1]);
    }
    
    /* Backward difference at right boundary */
    y_deriv[n-1] = (y_smooth[n-1] - y_smooth[n-2]) / (x[n-1] - x[n-2]);
}

/* Compute Tikhonov functional J = ||y - u||² + λ||D²u||² */
static void compute_functional(double *x, double *y, double *y_smooth, int n, double lambda,
                              double *data_term, double *reg_term, double *total_functional)
{
    double h_avg;
    int uniform = is_uniform_grid(x, n, &h_avg, 1e-10);
    
    /* Data fidelity term: ||y - u||² */
    *data_term = 0.0;
    for (int i = 0; i < n; i++) {
        double residual = y[i] - y_smooth[i];
        *data_term += residual * residual;
    }
    
    /* Regularization term: λ||D²u||² */
    *reg_term = 0.0;
    
    if (lambda > 0.0 && n >= 3) {
        if (uniform) {
            /* Uniform grid - use standard second differences */
            double h2 = h_avg * h_avg;
            
            /* Interior points: D²u ≈ (u[i-1] - 2u[i] + u[i+1])/h² */
            for (int i = 1; i < n-1; i++) {
                double d2u = (y_smooth[i-1] - 2.0*y_smooth[i] + y_smooth[i+1]) / h2;
                *reg_term += d2u * d2u;
            }
            
            /* Boundary points: use first-order differences for natural BC */
            /* Left boundary: approximate d²u ≈ (u[1] - u[0])/h² */
            double d2u_left = (y_smooth[1] - y_smooth[0]) / h2;
            *reg_term += 0.5 * d2u_left * d2u_left;  /* Reduced weight */
            
            /* Right boundary: approximate d²u ≈ (u[n-1] - u[n-2])/h² */
            double d2u_right = (y_smooth[n-1] - y_smooth[n-2]) / h2;
            *reg_term += 0.5 * d2u_right * d2u_right;  /* Reduced weight */
            
        } else {
            /* Non-uniform grid - use local spacing */
            for (int i = 1; i < n-1; i++) {
                double h1 = x[i] - x[i-1];
                double h2 = x[i+1] - x[i];
                double h_avg_local = (h1 + h2) / 2.0;
                
                double d2u = (y_smooth[i-1] - 2.0*y_smooth[i] + y_smooth[i+1]) / (h_avg_local * h_avg_local);
                *reg_term += d2u * d2u;
            }
            
            /* Boundary contributions for non-uniform grid */
            if (n >= 2) {
                double h_start = x[1] - x[0];
                double d2u_left = (y_smooth[1] - y_smooth[0]) / (h_start * h_start);
                *reg_term += 0.5 * d2u_left * d2u_left;
                
                double h_end = x[n-1] - x[n-2];
                double d2u_right = (y_smooth[n-1] - y_smooth[n-2]) / (h_end * h_end);
                *reg_term += 0.5 * d2u_right * d2u_right;
            }
        }
        
        *reg_term *= lambda;
    }
    
    /* Total functional */
    *total_functional = *data_term + *reg_term;
}

/* Main Tikhonov smoothing function with adaptive weights control - public interface */
TikhonovResult* tikhonov_smooth_adaptive(double *x, double *y, int n, double lambda, int adaptive_weights)
{
    return tikhonov_smooth_internal(x, y, n, lambda, 1, adaptive_weights);
}

/* Original main Tikhonov smoothing function - public interface */
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda)
{
    return tikhonov_smooth_internal(x, y, n, lambda, 1, 0);  /* adaptive_weights = 0 */
}

/* Internal smoothing function with options to suppress grid analysis output and control adaptive weights */
static TikhonovResult* tikhonov_smooth_internal(double *x, double *y, int n, double lambda, int print_grid_info, int adaptive_weights)
{
    TikhonovResult *result;
    double *AB;  /* Band matrix storage */
    double *b;   /* Right-hand side (copy of y) */
    int kd = 1;    /* Number of superdiagonals for tridiagonal matrix */
    int ldab = kd + 1;  /* Leading dimension for band storage */
    int nrhs = 1;  /* Number of right-hand sides */
    int info;      /* LAPACK info */
    int inc = 1;   /* Increment for BLAS */
    char uplo = 'U';  /* Upper triangular storage */
    
    /* Input validation */
    if (x == NULL || y == NULL || n < 1 || lambda < 0) {
        fprintf(stderr, "Error: Invalid input parameters to tikhonov_smooth\n");
        return NULL;
    }
    
    /* Check that x is monotonic */
    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "Error: x array must be strictly increasing\n");
            return NULL;
        }
    }
    
    /* Analyze grid uniformity using the new module */
    GridAnalysis *grid_info = analyze_grid(x, n, 0);
    if (grid_info == NULL) {
        fprintf(stderr, "Error: Grid analysis failed\n");
        return NULL;
    }
    
    /* Print grid analysis header only if requested */
    if (print_grid_info) {
        print_grid_analysis(grid_info, 1, "# ");
        printf("# Adaptive weights: %s\n", adaptive_weights ? "ENABLED" : "DISABLED");
    }
    
    /* Allocate result structure */
    result = (TikhonovResult *)malloc(sizeof(TikhonovResult));
    if (result == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for TikhonovResult\n");
        free_grid_analysis(grid_info);
        return NULL;
    }
    
    result->n = n;
    result->lambda = lambda;
    result->y_smooth = (double *)malloc(n * sizeof(double));
    result->y_deriv = (double *)malloc(n * sizeof(double));
    
    if (result->y_smooth == NULL || result->y_deriv == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for result arrays\n");
        free_tikhonov_result(result);
        free_grid_analysis(grid_info);
        return NULL;
    }
    
    /* Allocate working arrays */
    AB = (double *)calloc(ldab * n, sizeof(double));  /* Band matrix */
    b = (double *)malloc(n * sizeof(double));         /* RHS vector */
    
    if (AB == NULL || b == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for working arrays\n");
        /* Report memory usage for large datasets */
        if (n > 10000) {
            double mb_required = (double)(ldab * n + n) * sizeof(double) / (1024.0 * 1024.0);
            fprintf(stderr, "Required memory: %.2f MB for n=%d\n", mb_required, n);
        }
        free(AB);
        free(b);
        free_tikhonov_result(result);
        free_grid_analysis(grid_info);
        return NULL;
    }
    
    /* Copy y to b using BLAS */
    dcopy_(&n, y, &inc, b, &inc);
    
    /* Build band matrix A = I + lambda*D^T*D */
    build_band_matrix(x, n, lambda, AB, ldab, kd, adaptive_weights);
    
    /* Solve the system using LAPACK band solver */
    dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
    
    if (info != 0) {
        fprintf(stderr, "Error: LAPACK dpbsv failed with info = %d\n", info);
        if (info > 0) {
            fprintf(stderr, "Matrix is not positive definite (leading minor %d)\n", info);
            fprintf(stderr, "Try increasing lambda (current: %g)\n", lambda);
        } else {
            fprintf(stderr, "Invalid arguments to LAPACK\n");
        }
        free(AB);
        free(b);
        free_tikhonov_result(result);
        free_grid_analysis(grid_info);
        return NULL;
    }
    
    /* Copy solution to result */
    dcopy_(&n, b, &inc, result->y_smooth, &inc);
    
    /* Compute derivatives */
    compute_derivatives(x, result->y_smooth, n, result->y_deriv);
    
    /* Compute Tikhonov functional for user information */
    double data_term, reg_term, total_functional;
    compute_functional(x, y, result->y_smooth, n, lambda, &data_term, &reg_term, &total_functional);
    
    /* Store functional components in result structure */
    result->data_term = data_term;
    result->regularization_term = reg_term;
    result->total_functional = total_functional;
    
    /* Clean up */
    free(AB);
    free(b);
    free_grid_analysis(grid_info);
    
    return result;
}

/* Compute GCV score for given lambda with optional debugging */
static double compute_gcv_score_debug(double *x, double *y, int n, double lambda, int verbose)
{
    TikhonovResult *result;
    double rss = 0.0;
    double trace_H;
    double gcv_score;
    
    /* Apply smoothing - suppress grid info output during GCV optimization */
    result = tikhonov_smooth_internal(x, y, n, lambda, 0, 0);  /* Use default method for GCV */
    if (result == NULL) {
        return 1e20;  /* Return large value on error */
    }
    
    /* Compute residual sum of squares */
    for (int i = 0; i < n; i++) {
        double residual = y[i] - result->y_smooth[i];
        rss += residual * residual;
    }
    
    /* Better trace estimation with corrected eigenvalue formula */
    double h_avg = (x[n-1] - x[0]) / (n-1);
    double scale = lambda / (h_avg * h_avg);
    
    if (n <= 5000) {
        /* More accurate eigenvalue-based approximation */
        trace_H = 2.0;  /* Boundary contributions */
        for (int k = 1; k <= n-2; k++) {
            double theta = M_PI * k / (n+1);  /* Corrected denominator */
            double eigenval = 4.0 * pow(sin(theta/2.0), 2) / (h_avg * h_avg);
            trace_H += 1.0 / (1.0 + lambda * eigenval);
        }
    } else {
        /* Fast approximation for large problems with better scaling */
        trace_H = (double)n / (1.0 + sqrt(scale));
        if (trace_H < 2.0) trace_H = 2.0;
        if (trace_H > n-2) trace_H = n-2;
    }
    
    /* Compute GCV score */
    double denom = 1.0 - trace_H / n;
    if (denom > 1e-8) {
        gcv_score = (rss / n) / (denom * denom);
    } else {
        gcv_score = 1e20;
    }
    
    /* Debug output */
    if (verbose) {
        printf("# λ=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f, GCV=%9.3e, Data=%9.3e, Reg=%9.3e\n",
               lambda, result->total_functional, rss, trace_H, gcv_score, 
               result->data_term, result->regularization_term);
    }
    
    free_tikhonov_result(result);
    return gcv_score;
}

/* Find optimal lambda using GCV with enhanced diagnostics */
double find_optimal_lambda_gcv(double *x, double *y, int n)
{
    double best_lambda = 0.01;  /* Conservative default */
    double best_gcv = 1e20;
    
    /* Input validation */
    if (n < 3) {
        fprintf(stderr, "Warning: Too few points for GCV optimization (n=%d)\n", n);
        return best_lambda;
    }
    
    /* Analyze grid uniformity first */
    GridAnalysis *grid_info = analyze_grid(x, n, 0);
    if (grid_info == NULL) {
        fprintf(stderr, "Warning: Grid analysis failed in GCV optimization\n");
        return best_lambda;
    }
    
    /* Print data statistics */
    double y_min = y[0], y_max = y[0], y_mean = 0.0;
    for (int i = 0; i < n; i++) {
        if (y[i] < y_min) y_min = y[i];
        if (y[i] > y_max) y_max = y[i];
        y_mean += y[i];
    }
    y_mean /= n;
    
    double y_var = 0.0;
    for (int i = 0; i < n; i++) {
        double dev = y[i] - y_mean;
        y_var += dev * dev;
    }
    y_var /= (n-1);
    
    printf("# GCV optimization: n=%d, y_range=[%.3f,%.3f], y_mean=%.3f, y_std=%.3f\n", 
           n, y_min, y_max, y_mean, sqrt(y_var));
    
    /* Adjust lambda search range based on grid uniformity */
    double lambda_min = 1e-6;
    double lambda_max = 1e0;
    
    if (grid_info->ratio_max_min > 10.0) {
        /* For highly non-uniform grids, suggest smaller lambda range */
        lambda_min = 1e-8;
        lambda_max = 1e-2;
        printf("# Adjusted lambda range for non-uniform grid: [%.1e, %.1e]\n", 
               lambda_min, lambda_max);
    }
    
    printf("# Format: λ | Functional | RSS | tr(H) | GCV | Data_term | Reg_term\n");
    
    if (n > 20000) {
        /* For very large datasets, use simple grid search with fewer points */
        printf("# Large dataset: using coarse grid search\n");
        
        double lambdas[] = {1e-4, 1e-3, 1e-2, 1e-1, 1e0};
        int n_lambdas = sizeof(lambdas) / sizeof(lambdas[0]);
        
        for (int i = 0; i < n_lambdas; i++) {
            double gcv = compute_gcv_score_debug(x, y, n, lambdas[i], 1);
            if (gcv < best_gcv) {
                best_gcv = gcv;
                best_lambda = lambdas[i];
            }
        }
        
    } else {
        /* Enhanced grid search for smaller problems */
        printf("# Phase 1: Coarse grid search\n");
        
        /* Generate logarithmically spaced lambda values */
        double log_min = log10(lambda_min);
        double log_max = log10(lambda_max);
        int n_points = 13;
        
        double best_functional = 1e20;
        double best_functional_lambda = 0.01;
        
        for (int i = 0; i < n_points; i++) {
            double log_lambda = log_min + (log_max - log_min) * i / (n_points - 1);
            double lambda_test = pow(10.0, log_lambda);
            
            double gcv = compute_gcv_score_debug(x, y, n, lambda_test, 1);
            
            /* Track which lambda gives best functional (for comparison) */
            TikhonovResult *test_result = tikhonov_smooth_internal(x, y, n, lambda_test, 0, 0);
            if (test_result) {
                if (test_result->total_functional < best_functional) {
                    best_functional = test_result->total_functional;
                    best_functional_lambda = lambda_test;
                }
                free_tikhonov_result(test_result);
            }
            
            if (gcv < best_gcv) {
                best_gcv = gcv;
                best_lambda = lambda_test;
            }
        }
        
        printf("# Best λ by GCV: %.6e, Best λ by functional: %.6e\n", 
               best_lambda, best_functional_lambda);
        
        /* Phase 2: Refinement around best lambda (only for smaller problems) */
        if (n <= 5000) {
            printf("# Phase 2: Refinement around λ=%.6e\n", best_lambda);
            
            for (int i = 0; i < 8; i++) {
                double factor = 0.3 + 1.4 * i / 7.0;  /* 0.3x to 1.7x around best */
                double lambda_test = best_lambda * factor;
                
                if (lambda_test > lambda_min && lambda_test < lambda_max) {
                    double gcv = compute_gcv_score_debug(x, y, n, lambda_test, 1);
                    
                    if (gcv < best_gcv) {
                        best_gcv = gcv;
                        best_lambda = lambda_test;
                    }
                }
            }
        }
    }
    
    printf("# Final optimal λ: %.6e (GCV=%.3e)\n", best_lambda, best_gcv);
    
    /* Sanity check and warning */
    TikhonovResult *final_result = tikhonov_smooth_internal(x, y, n, best_lambda, 0, 0);
    if (final_result) {
        double reg_ratio = final_result->regularization_term / final_result->total_functional;
        double data_ratio = final_result->data_term / final_result->total_functional;
        
        printf("# Final balance: Data=%.1f%%, Regularization=%.1f%%\n", 
               data_ratio*100, reg_ratio*100);
        
        if (reg_ratio > 0.98) {
            printf("# WARNING: Regularization strongly dominates (%.1f%%)\n", reg_ratio*100);
            printf("# This suggests over-smoothing. Consider manual λ values:\n");
            printf("# Try: λ=0.1, 0.01, 0.001 and compare results\n");
        } else if (data_ratio > 0.98) {
            printf("# WARNING: Data term dominates - possibly under-smoothed\n");
            printf("# Consider smaller λ for more smoothing\n");
        }
        
        /* Additional warning for non-uniform grids */
        if (grid_info->reliability_warning) {
            printf("# ADDITIONAL WARNING: Due to grid non-uniformity issues,\n");
            printf("# automatic lambda selection may be suboptimal.\n");
            printf("# Manual tuning is recommended for best results.\n");
        }
        
        free_tikhonov_result(final_result);
    }
    
    free_grid_analysis(grid_info);
    
    return best_lambda;
}

/* Free allocated memory */
void free_tikhonov_result(TikhonovResult *result)
{
    if (result != NULL) {
        free(result->y_smooth);
        free(result->y_deriv);
        free(result);
    }
}

