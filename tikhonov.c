/*  Tikhonov regularization for data smoothing
 *  Simplified version with correct D² operator for non-uniform grids
 *  V4.0/2025-10-04/ Simplified implementation with correct discretization
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

/* Build band matrix A = I + lambda*D^T*D with correct D² operator */
static void build_band_matrix(double *x, int n, double lambda, double *AB, int ldab, int kd)
{
    int j;
    
    /* Clear the band matrix */
    memset(AB, 0, ldab * n * sizeof(double));
    
    /* Set main diagonal to identity */
    for (j = 0; j < n; j++) {
        AB[kd + j*ldab] = 1.0;
    }
    
    if (lambda <= 0.0 || n < 3) return;
    
    /* Interior points: use local spacing for D² operator
     * For point j with left spacing h_left and right spacing h_right:
     * D²u ≈ 2/(h_left+h_right) * [u[j-1]/h_left - u[j]*(1/h_left+1/h_right) + u[j+1]/h_right]
     */
    for (j = 1; j < n-1; j++) {
        double h_left = x[j] - x[j-1];
        double h_right = x[j+1] - x[j];
        double h_sum = h_left + h_right;
        
        /* Weights for second difference */
        double w = 2.0 * lambda / h_sum;
        
        AB[kd + j*ldab] += w * (1.0/h_left + 1.0/h_right);  /* diagonal */
        AB[0 + (j+1)*ldab] += -w / h_right;                  /* superdiagonal */
    }
    
    /* Boundary points: use first-order differences for natural BC */
    double h_first = x[1] - x[0];
    AB[kd + 0*ldab] += lambda / (h_first * h_first);
    
    double h_last = x[n-1] - x[n-2];
    AB[kd + (n-1)*ldab] += lambda / (h_last * h_last);
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
    /* Data fidelity term: ||y - u||² */
    *data_term = 0.0;
    for (int i = 0; i < n; i++) {
        double residual = y[i] - y_smooth[i];
        *data_term += residual * residual;
    }
    
    /* Regularization term: λ||D²u||² using local spacings */
    *reg_term = 0.0;
    
    if (lambda > 0.0 && n >= 3) {
        /* Interior points */
        for (int i = 1; i < n-1; i++) {
            double h1 = x[i] - x[i-1];
            double h2 = x[i+1] - x[i];
            double h_avg = (h1 + h2) / 2.0;
            
            double d2u = (y_smooth[i-1] - 2.0*y_smooth[i] + y_smooth[i+1]) / (h_avg * h_avg);
            *reg_term += d2u * d2u;
        }
        
        /* Boundary contributions */
        if (n >= 2) {
            double h_start = x[1] - x[0];
            double d2u_left = (y_smooth[1] - y_smooth[0]) / (h_start * h_start);
            *reg_term += 0.5 * d2u_left * d2u_left;
            
            double h_end = x[n-1] - x[n-2];
            double d2u_right = (y_smooth[n-1] - y_smooth[n-2]) / (h_end * h_end);
            *reg_term += 0.5 * d2u_right * d2u_right;
        }
        
        *reg_term *= lambda;
    }
    
    /* Total functional */
    *total_functional = *data_term + *reg_term;
}

/* Main Tikhonov smoothing function */
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda)
{
    TikhonovResult *result;
    double *AB, *b;
    int kd = 1;
    int ldab = kd + 1;
    int nrhs = 1;
    int info;
    int inc = 1;
    char uplo = 'U';
    
    /* Input validation */
    if (x == NULL || y == NULL || n < 1 || lambda < 0) {
        fprintf(stderr, "Error: Invalid input parameters\n");
        return NULL;
    }
    
    /* Check monotonicity */
    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "Error: x array must be strictly increasing\n");
            return NULL;
        }
    }
    
    /* Grid analysis - show warnings if significant non-uniformity */
    GridAnalysis *grid_info = analyze_grid(x, n, 0);
    if (grid_info) {
        if (grid_info->reliability_warning) {
            printf("# Grid analysis:\n");
            print_grid_analysis(grid_info, 0, "# ");
        }
        free_grid_analysis(grid_info);
    }
    
    /* Allocate result */
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
    
    /* Allocate working arrays */
    AB = (double *)calloc(ldab * n, sizeof(double));
    b = (double *)malloc(n * sizeof(double));
    
    if (!AB || !b) {
        free(AB); 
        free(b);
        free_tikhonov_result(result);
        return NULL;
    }
    
    /* Copy y to b */
    dcopy_(&n, y, &inc, b, &inc);
    
    /* Build matrix with correct D² operator */
    build_band_matrix(x, n, lambda, AB, ldab, kd);
    
    /* Solve the system */
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
    
    /* Copy solution */
    dcopy_(&n, b, &inc, result->y_smooth, &inc);
    
    /* Compute derivatives and functional */
    compute_derivatives(x, result->y_smooth, n, result->y_deriv);
    compute_functional(x, y, result->y_smooth, n, lambda,
                      &result->data_term, &result->regularization_term, 
                      &result->total_functional);
    
    free(AB);
    free(b);
    return result;
}

/* Compute GCV score for given lambda */
static double compute_gcv_score(double *x, double *y, int n, double lambda, int verbose)
{
    TikhonovResult *result;
    double rss = 0.0;
    double trace_H;
    double gcv_score;
    
    /* Apply smoothing */
    result = tikhonov_smooth(x, y, n, lambda);
    if (result == NULL) {
        return 1e20;
    }
    
    /* Compute residual sum of squares */
    for (int i = 0; i < n; i++) {
        double residual = y[i] - result->y_smooth[i];
        rss += residual * residual;
    }
    
    /* Trace estimation using eigenvalues for natural BC */
    double h_avg = (x[n-1] - x[0]) / (n-1);
    double scale = lambda / (h_avg * h_avg);
    
    if (n <= 5000) {
        /* Eigenvalue-based approximation for natural BC */
        trace_H = 2.0;  /* Boundary contributions */
        for (int k = 1; k <= n-2; k++) {
            double theta = M_PI * k / n;  /* Natural BC eigenvalues */
            double eigenval = 4.0 * pow(sin(theta/2.0), 2) / (h_avg * h_avg);
            trace_H += 1.0 / (1.0 + lambda * eigenval);
        }
    } else {
        /* Fast approximation for large problems */
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
        printf("# λ=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f, GCV=%9.3e\n",
               lambda, result->total_functional, rss, trace_H, gcv_score);
    }
    
    free_tikhonov_result(result);
    return gcv_score;
}

/* Find optimal lambda using GCV */
double find_optimal_lambda_gcv(double *x, double *y, int n)
{
    double best_lambda = 0.01;
    double best_gcv = 1e20;
    
    if (n < 3) {
        fprintf(stderr, "Warning: Too few points for GCV (n=%d)\n", n);
        return best_lambda;
    }
    
    printf("# GCV optimization for n=%d points\n", n);
    printf("# Format: λ | Functional | RSS | tr(H) | GCV\n");
    
    /* Logarithmic grid search */
    double lambda_min = 1e-6;
    double lambda_max = 1e0;
    double log_min = log10(lambda_min);
    double log_max = log10(lambda_max);
    int n_points = 13;
    
    for (int i = 0; i < n_points; i++) {
        double log_lambda = log_min + (log_max - log_min) * i / (n_points - 1);
        double lambda_test = pow(10.0, log_lambda);
        
        double gcv = compute_gcv_score(x, y, n, lambda_test, 1);
        
        if (gcv < best_gcv) {
            best_gcv = gcv;
            best_lambda = lambda_test;
        }
    }
    
    /* Refinement for smaller problems */
    if (n <= 5000) {
        printf("# Refinement around λ=%.6e\n", best_lambda);
        
        for (int i = 0; i < 8; i++) {
            double factor = 0.3 + 1.4 * i / 7.0;
            double lambda_test = best_lambda * factor;
            
            if (lambda_test > lambda_min && lambda_test < lambda_max) {
                double gcv = compute_gcv_score(x, y, n, lambda_test, 1);
                
                if (gcv < best_gcv) {
                    best_gcv = gcv;
                    best_lambda = lambda_test;
                }
            }
        }
    }
    
    printf("# Optimal λ: %.6e (GCV=%.3e)\n", best_lambda, best_gcv);
    
    /* Check balance */
    TikhonovResult *final = tikhonov_smooth(x, y, n, best_lambda);
    if (final) {
        double reg_ratio = final->regularization_term / final->total_functional;
        printf("# Balance: Data=%.1f%%, Regularization=%.1f%%\n", 
               (1.0-reg_ratio)*100, reg_ratio*100);
        
        if (reg_ratio > 0.95) {
            printf("# WARNING: Over-smoothing detected. Try manual λ values.\n");
        }
        free_tikhonov_result(final);
    }
    
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
