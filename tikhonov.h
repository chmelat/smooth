/*  Tikhonov regularization for data smoothing
 *  Header file with functional support
 *  V3.1/2025-05-26/ Production version with optional adaptive weights
 */

#ifndef TIKHONOV_H
#define TIKHONOV_H

/* Structure for Tikhonov results */
typedef struct {
    double *y_smooth;            /* Smoothed values */
    double *y_deriv;             /* First derivatives */
    double lambda;               /* Used regularization parameter */
    int n;                       /* Number of points */
    double data_term;            /* Data fidelity term: ||y - u||² */
    double regularization_term;  /* Regularization term: λ||D²u||² */
    double total_functional;     /* Total functional: data_term + regularization_term */
} TikhonovResult;

/* Main smoothing function
 * 
 * Performs Tikhonov regularization using second-order penalty term.
 * Solves: min ||y - u||² + λ||D²u||²
 * where D² is the second-order difference operator.
 *
 * Parameters:
 *   x      - Array of x-coordinates (must be strictly monotonic increasing)
 *   y      - Array of y-values to be smoothed
 *   n      - Number of data points (must be >= 1)
 *   lambda - Regularization parameter (>= 0)
 *            Small lambda: less smoothing, follows data closely
 *            Large lambda: more smoothing, smoother result
 * 
 * Returns:
 *   Pointer to TikhonovResult structure containing:
 *   - y_smooth: smoothed values
 *   - y_deriv: first derivatives
 *   - lambda: used regularization parameter
 *   - n: number of points
 *   - data_term: ||y - u||² (data fidelity)
 *   - regularization_term: λ||D²u||² (smoothness penalty)
 *   - total_functional: sum of both terms
 *   Returns NULL on error.
 * 
 * Notes:
 *   - Uses natural boundary conditions (second derivative = 0 at ends)
 *   - Efficient band matrix implementation for O(n) memory complexity
 *   - Uses average coefficient method for non-uniform grids (default)
 *   - Memory must be freed using free_tikhonov_result()
 *   - For large datasets (n > 50000), uses optimized algorithms
 */
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda);

/* Main smoothing function with adaptive weights control
 * 
 * Performs Tikhonov regularization with optional local weighting for non-uniform grids.
 *
 * Parameters:
 *   x              - Array of x-coordinates (must be strictly monotonic increasing)
 *   y              - Array of y-values to be smoothed
 *   n              - Number of data points (must be >= 1)
 *   lambda         - Regularization parameter (>= 0)
 *   adaptive_weights - 0: use average coefficient method (default, faster)
 *                     1: use local weights for non-uniform grids (more accurate)
 * 
 * Returns:
 *   Same as tikhonov_smooth()
 * 
 * Notes:
 *   - adaptive_weights=0: Uses harmonic mean of 1/h² for all intervals (faster)
 *   - adaptive_weights=1: Uses individual 1/h² weight for each interval (slower but more accurate)
 *   - For uniform grids, both methods are equivalent
 *   - Local weighting is recommended for highly non-uniform grids (h_max/h_min > 10)
 */
TikhonovResult* tikhonov_smooth_adaptive(double *x, double *y, int n, double lambda, int adaptive_weights);

/* Automatic lambda selection using Generalized Cross Validation (GCV)
 * 
 * Finds the optimal regularization parameter by minimizing the GCV criterion:
 * GCV(λ) = (RSS(λ)/n) / ((1 - tr(H(λ))/n)²)
 * where RSS is residual sum of squares and tr(H) is trace of influence matrix.
 *
 * Parameters:
 *   x - Array of x-coordinates (must be strictly monotonic increasing)
 *   y - Array of y-values
 *   n - Number of data points (should be >= 5 for reliable results)
 * 
 * Returns:
 *   Optimal lambda value minimizing GCV criterion
 * 
 * Notes:
 *   - Uses enhanced grid search with diagnostics
 *   - Prints detailed optimization information to stdout
 *   - Search range: 1e-8 to 1e3
 *   - For small datasets (n < 5), returns conservative default
 *   - For large datasets (n > 20000), uses coarse grid search
 *   - Warns if regularization term dominates excessively
 *   - Uses average coefficient method (not adaptive weights)
 *   - Computational cost: O(n × log(search_range) × iterations)
 * 
 * Example diagnostic output:
 *   # GCV optimization: n=1000, y_range=[1.2,5.8], y_mean=3.4, y_std=1.1
 *   # λ=1.000e-03: J=1.234e+02, RSS=1.200e+01, tr(H)=85.2, GCV=2.345e+01
 *   # Final optimal λ: 1.000e-02 (GCV=1.876e+01)
 *   # Final balance: Data=75.2%, Regularization=24.8%
 */
double find_optimal_lambda_gcv(double *x, double *y, int n);

/* Free allocated memory for TikhonovResult structure
 * 
 * Parameters:
 *   result - Pointer to TikhonovResult structure to be freed
 * 
 * Notes:
 *   - Safe to call with NULL pointer
 *   - Frees all internal arrays and the structure itself
 *   - Always call this function to avoid memory leaks
 */
void free_tikhonov_result(TikhonovResult *result);

/* Compilation requirements:
 * 
 * Requires LAPACK and BLAS libraries. Compile with:
 *   gcc -o program program.c tikhonov.c -llapack -lblas -lm
 * 
 * For optimized performance with Intel MKL:
 *   gcc -o program program.c tikhonov.c -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm
 * 
 * Example usage:
 * 
 * #include "tikhonov.h"
 * 
 * int main() {
 *     double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
 *     double y[] = {1.0, 2.1, 3.9, 6.2, 8.1};  // Noisy data
 *     int n = 5;
 *     
 *     // Method 1: Manual lambda selection with default method
 *     TikhonovResult *result = tikhonov_smooth(x, y, n, 0.1);
 *     
 *     // Method 2: Manual lambda with adaptive weights for non-uniform grids
 *     TikhonovResult *result = tikhonov_smooth_adaptive(x, y, n, 0.1, 1);
 *     
 *     // Method 3: Automatic lambda selection with diagnostics
 *     double optimal_lambda = find_optimal_lambda_gcv(x, y, n);
 *     TikhonovResult *result = tikhonov_smooth(x, y, n, optimal_lambda);
 *     
 *     if (result != NULL) {
 *         printf("# Tikhonov Results (lambda = %.2e)\n", result->lambda);
 *         printf("# Functional J = %.3e (Data: %.3e + Reg: %.3e)\n", 
 *                result->total_functional, result->data_term, result->regularization_term);
 *         printf("# Data/Total = %.1f%%, Reg/Total = %.1f%%\n",
 *                100.0 * result->data_term / result->total_functional,
 *                100.0 * result->regularization_term / result->total_functional);
 *         
 *         for (int i = 0; i < n; i++) {
 *             printf("%8.4f %8.4f %8.4f\n", 
 *                    x[i], result->y_smooth[i], result->y_deriv[i]);
 *         }
 *         
 *         free_tikhonov_result(result);
 *     }
 *     
 *     return 0;
 * }
 * 
 * Performance characteristics:
 *   - Time complexity: O(n) for smoothing, O(n·log(search_range)) for GCV
 *   - Space complexity: O(n) - uses efficient band matrix storage
 *   - Suitable for datasets from n=3 to n=10⁶⁺
 *   - Numerical stability: Good for λ ∈ [1e-8, 1e3]
 *   - Memory usage: ~16n bytes (vs ~8n² for full matrix methods)
 * 
 * Algorithm details:
 *   - Uses LAPACK dpbsv for efficient band matrix solving
 *   - Second-order finite differences for regularization operator
 *   - Natural boundary conditions (free ends, d²u/dx²|endpoints = 0)
 *   - Enhanced GCV with eigenvalue-based trace estimation
 *   - Adaptive algorithms based on problem size
 *   - Comprehensive error checking and diagnostics
 * 
 * Functional interpretation:
 *   - data_term: measures fit to original data (smaller = better fit)
 *   - regularization_term: measures smoothness (larger λ = more smooth)
 *   - total_functional: objective function being minimized
 *   - Typical good balance: data_term 50-90% of total
 *   - Warning signs: regularization_term > 95% (over-smoothing)
 * 
 * Grid handling strategies:
 *   - Uniform grids: Exact second-order finite differences with spacing h
 *   - Non-uniform grids (default): Average coefficient method - robust and fast
 *   - Non-uniform grids (adaptive): Local weighting - more accurate for highly non-uniform data
 *   - Recommendation: Use adaptive weights when h_max/h_min > 10 or CV > 0.5
 */

#endif /* TIKHONOV_H */

