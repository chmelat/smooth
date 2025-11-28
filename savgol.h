/*  Savitzky-Golay filter for data smoothing
 *  Header file
 *  V2.4/2025-11-28/ FIXED: Allow deriv_order > poly_degree (returns zero coefficients)
 *  V2.3/2025-11-28/ savgol_coefficients now returns error code
 *  V2.1/2025-10-04/ Added uniform grid requirement documentation
 */

#ifndef SAVGOL_H
#define SAVGOL_H

#include "grid_analysis.h"

/* Structure for Savitzky-Golay results */
typedef struct {
    double *y_smooth;     /* Smoothed values */
    double *y_deriv;      /* First derivatives */
    int n;                /* Number of points */
    int poly_degree;      /* Polynomial degree used */
    int window_size;      /* Window size used */
} SavgolResult;

/* Main Savitzky-Golay smoothing function
 * 
 * Performs Savitzky-Golay filtering for data smoothing and derivative estimation.
 * This method provides optimal polynomial smoothing in the least-squares sense.
 *
 * IMPORTANT: This method REQUIRES uniformly spaced data points!
 * The function will check grid uniformity and REJECT non-uniform grids.
 *
 * Parameters:
 *   x            - Array of x-coordinates (must be strictly monotonic increasing)
 *   y            - Array of y-values to be smoothed
 *   n            - Number of data points
 *   window_size  - Size of sliding window (must be odd, >= 3)
 *   poly_degree  - Degree of approximating polynomial (0 to 6)
 *   grid_info    - Grid analysis results (used for uniformity check)
 * 
 * Returns:
 *   Pointer to SavgolResult structure containing smoothed values and derivatives
 *   Returns NULL on error or if grid is non-uniform.
 * 
 * Grid Uniformity Requirements:
 *   - The method checks coefficient of variation (CV) of spacing
 *   - CV = std_dev(spacing) / avg(spacing)
 *   - Grid is rejected if CV > 0.05 (5% variation)
 *   - For non-uniform grids, use Tikhonov or Polyfit methods instead
 * 
 * Notes:
 *   - Memory must be freed using free_savgol_result()
 *   - Uses asymmetric windows at boundaries
 *   - Provides better derivative estimates than simple finite differences
 *   - For non-uniform grids, error message will recommend alternatives
 * 
 * Algorithm:
 *   - Computes convolution coefficients for polynomial fitting
 *   - Assumes uniform spacing (integer indices)
 *   - Applies filter to each point with appropriate window
 *   - Derivative scaling assumes uniform grid
 * 
 * Error handling:
 *   If grid is non-uniform, function returns NULL and prints:
 *   - Grid statistics (CV, h_min, h_max, ratio)
 *   - Recommended alternative methods
 *   - Suggestions for data preprocessing
 * 
 * Example usage:
 *   // This will work (uniform grid):
 *   double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
 *   double y[] = {1.0, 2.1, 2.9, 4.2, 5.1, 6.0};
 *   SavgolResult *result = savgol_smooth(x, y, 6, 5, 2);
 * 
 *   // This will FAIL (non-uniform grid):
 *   double x[] = {0.0, 1.0, 1.5, 4.0, 5.0, 10.0};
 *   double y[] = {1.0, 2.1, 2.9, 4.2, 5.1, 6.0};
 *   SavgolResult *result = savgol_smooth(x, y, 6, 5, 2);
 *   // Returns NULL with error message suggesting Tikhonov method
 */
SavgolResult* savgol_smooth(double *x, double *y, int n, int window_size, int poly_degree,
                            GridAnalysis *grid_info);

/* Calculate Savitzky-Golay coefficients
 *
 * Computes the convolution coefficients for a given window configuration.
 * This is a low-level function used internally by savgol_smooth().
 *
 * Parameters:
 *   nl           - Number of points to the left of center
 *   nr           - Number of points to the right of center
 *   poly_degree  - Degree of polynomial
 *   deriv_order  - Order of derivative (0 for smoothing, 1 for first derivative)
 *   c            - Output array for coefficients (must be pre-allocated with size nl+nr+1)
 *
 * Returns:
 *   0 on success, -1 on error.
 *   On error, output array 'c' is zeroed for safety.
 *
 * Special cases:
 *   - If deriv_order > poly_degree, returns 0 with zero coefficients
 *     (mathematically correct: derivative of lower-degree polynomial is zero)
 *
 * Notes:
 *   - For symmetric windows: nl = nr
 *   - For boundary handling: nl != nr
 *   - Total window size = nl + nr + 1
 *   - Assumes integer spacing (normalized coordinates)
 */
int savgol_coefficients(int nl, int nr, int poly_degree, int deriv_order, double *c);

/* Free allocated memory for SavgolResult structure
 * 
 * Parameters:
 *   result - Pointer to SavgolResult structure to be freed
 */
void free_savgol_result(SavgolResult *result);

/* Recommendations for non-uniform grids:
 * 
 * If savgol_smooth() rejects your data due to non-uniform spacing, consider:
 * 
 * 1. TIKHONOV METHOD (best for non-uniform grids):
 *    TikhonovResult *result = tikhonov_smooth(x, y, n, lambda);
 *    - Works correctly with any spacing
 *    - Global optimization approach
 *    - Use automatic lambda: find_optimal_lambda_gcv()
 * 
 * 2. POLYFIT METHOD (local fitting):
 *    PolyfitResult *result = polyfit_smooth(x, y, n, window_size, poly_degree);
 *    - Less sensitive to spacing variations
 *    - Local polynomial fitting
 *    - Similar to SG but more robust
 * 
 * 3. RESAMPLE TO UNIFORM GRID:
 *    - Interpolate data to uniform spacing
 *    - Then apply Savitzky-Golay
 *    - Be aware this adds interpolation error
 * 
 * Grid uniformity check details:
 *   CV < 0.01  : Perfectly uniform, SG optimal
 *   CV < 0.05  : Nearly uniform, SG acceptable (warning shown)
 *   CV >= 0.05 : Non-uniform, SG rejected (error shown)
 */

#endif /* SAVGOL_H */
