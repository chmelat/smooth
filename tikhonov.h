/*  Tikhonov regularization for data smoothing
 *  Second derivative penalty via (D²)ᵀ W D² Gram matrix
 *  V5.0/2026-02-07/ Corrected to true 2nd order penalty (pentadiagonal)
 */

#ifndef TIKHONOV_H
#define TIKHONOV_H

#include "grid_analysis.h"

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
 * where D² is the second-order difference operator with correct
 * discretization for non-uniform grids.
 *
 * Parameters:
 *   x         - Array of x-coordinates (must be strictly monotonic increasing)
 *   y         - Array of y-values to be smoothed
 *   n         - Number of data points (must be >= 1)
 *   lambda    - Regularization parameter (>= 0)
 *               Small lambda: less smoothing, follows data closely
 *               Large lambda: more smoothing, smoother result
 *   grid_info - Grid analysis results (required, non-NULL); the module uses a
 *               single integral-measure discretization and does not select
 *               between schemes based on grid statistics
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
 *   - Efficient pentadiagonal band matrix (bandwidth 2) for O(n) memory
 *   - Correct discretization for both uniform and non-uniform grids
 *   - Memory must be freed using free_tikhonov_result()
 */
TikhonovResult* tikhonov_smooth(const double *x, const double *y, int n, double lambda,
                                const GridAnalysis *grid_info);

/* Automatic lambda selection using Generalized Cross Validation (GCV)
 * 
 * Finds the optimal regularization parameter by minimizing:
 * GCV(λ) = (RSS(λ)/n) / ((1 - tr(H(λ))/n)²)
 *
 * Parameters:
 *   x         - Array of x-coordinates (must be strictly monotonic increasing)
 *   y         - Array of y-values
 *   n         - Number of data points (should be >= 5 for reliable results)
 *   grid_info - Grid analysis results (used for GCV optimization warnings)
 * 
 * Returns:
 *   Optimal lambda value minimizing GCV criterion
 * 
 * Notes:
 *   - Uses grid search with refinement (refinement only runs for n <= 5000;
 *     see tikhonov.c for the exact thresholds)
 *   - Prints detailed optimization information to stdout
 *   - Search range: 1e-8 to 1e0 (tikhonov.c is the source of truth)
 *   - For small datasets (n < 3), returns conservative default
 *   - Warns if the chosen lambda lands on the edge of the search range
 */
double find_optimal_lambda_gcv(const double *x, const double *y, int n, const GridAnalysis *grid_info);

/* Free allocated memory for TikhonovResult structure
 * 
 * Parameters:
 *   result - Pointer to TikhonovResult structure to be freed
 * 
 * Notes:
 *   - Safe to call with NULL pointer
 */
void free_tikhonov_result(TikhonovResult *result);

/* Example usage:
 *
 * #include <stdio.h>
 * #include "tikhonov.h"
 * #include "grid_analysis.h"
 *
 * int main(void) {
 *     double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
 *     double y[] = {1.0, 2.1, 3.9, 6.2, 8.1};
 *     int n = 5;
 *
 *     // Grid analysis is required; every method takes it as input.
 *     GridAnalysis *grid = analyze_grid(x, n);
 *     if (grid == NULL) return 1;
 *
 *     // Manual lambda:   tikhonov_smooth(x, y, n, 0.1, grid)
 *     // Automatic lambda (GCV):
 *     double lambda = find_optimal_lambda_gcv(x, y, n, grid);
 *     TikhonovResult *result = tikhonov_smooth(x, y, n, lambda, grid);
 *
 *     if (result != NULL) {
 *         printf("# Functional J = %.3e\n", result->total_functional);
 *
 *         for (int i = 0; i < n; i++) {
 *             printf("%8.4f %8.4f %8.4f\n",
 *                    x[i], result->y_smooth[i], result->y_deriv[i]);
 *         }
 *
 *         free_tikhonov_result(result);
 *     }
 *
 *     free_grid_analysis(grid);
 *     return 0;
 * }
 *
 * Compilation:
 *   gcc -o program program.c tikhonov.c grid_analysis.c -llapack -lblas -lm
 */

#endif /* TIKHONOV_H */
