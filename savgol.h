/*  Savitzky-Golay filter for data smoothing
 *  Header file
 *  V1.0/2025-05-27/ Extracted from smooth.c
 */

#ifndef SAVGOL_H
#define SAVGOL_H

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
 * Parameters:
 *   x            - Array of x-coordinates (must be strictly monotonic increasing)
 *   y            - Array of y-values to be smoothed
 *   n            - Number of data points
 *   window_size  - Size of sliding window (must be odd, >= 3)
 *   poly_degree  - Degree of approximating polynomial (0 to 6)
 * 
 * Returns:
 *   Pointer to SavgolResult structure containing smoothed values and derivatives
 *   Returns NULL on error.
 * 
 * Notes:
 *   - Memory must be freed using free_savgol_result()
 *   - Handles non-uniform grids with adaptive spacing
 *   - Uses asymmetric windows at boundaries
 *   - Provides better derivative estimates than simple finite differences
 * 
 * Algorithm:
 *   - Computes convolution coefficients for polynomial fitting
 *   - Applies filter to each point with appropriate window
 *   - Accounts for non-uniform spacing in derivative calculations
 */
SavgolResult* savgol_smooth(double *x, double *y, int n, int window_size, int poly_degree);

/* Calculate Savitzky-Golay coefficients
 * 
 * Computes the convolution coefficients for a given window configuration.
 *
 * Parameters:
 *   nl           - Number of points to the left of center
 *   nr           - Number of points to the right of center
 *   poly_degree  - Degree of polynomial
 *   deriv_order  - Order of derivative (0 for smoothing, 1 for first derivative)
 *   c            - Output array for coefficients (must be pre-allocated with size nl+nr+1)
 * 
 * Notes:
 *   - For symmetric windows: nl = nr
 *   - For boundary handling: nl != nr
 *   - Total window size = nl + nr + 1
 */
void savgol_coefficients(int nl, int nr, int poly_degree, int deriv_order, double *c);

/* Free allocated memory for SavgolResult structure
 * 
 * Parameters:
 *   result - Pointer to SavgolResult structure to be freed
 */
void free_savgol_result(SavgolResult *result);

#endif /* SAVGOL_H */
