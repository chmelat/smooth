/*  Polynomial fitting for data smoothing
 *  Header file for least squares polynomial approximation
 *  V1.0/2025-05-27/ Extracted from smooth.c
 */

#ifndef POLYFIT_H
#define POLYFIT_H

/* Structure for polynomial fitting results */
typedef struct {
    double *y_smooth;     /* Smoothed values */
    double *y_deriv;      /* First derivatives */
    double *coeffs;       /* Polynomial coefficients (if needed) */
    int n;                /* Number of points */
    int poly_degree;      /* Polynomial degree used */
    int window_size;      /* Window size used */
} PolyfitResult;

/* Main polynomial fitting function
 * 
 * Performs local polynomial fitting using least squares method
 * in a sliding window.
 *
 * Parameters:
 *   x            - Array of x-coordinates (must be strictly monotonic increasing)
 *   y            - Array of y-values to be smoothed
 *   n            - Number of data points
 *   window_size  - Size of sliding window (must be odd, >= 3)
 *   poly_degree  - Degree of approximating polynomial (0 to 6)
 * 
 * Returns:
 *   Pointer to PolyfitResult structure containing smoothed values and derivatives
 *   Returns NULL on error.
 * 
 * Notes:
 *   - Memory must be freed using free_polyfit_result()
 *   - Boundary points use asymmetric windows
 *   - Uses Gaussian elimination for solving normal equations
 */
PolyfitResult* polyfit_smooth(double *x, double *y, int n, int window_size, int poly_degree);

/* Free allocated memory for PolyfitResult structure
 * 
 * Parameters:
 *   result - Pointer to PolyfitResult structure to be freed
 */
void free_polyfit_result(PolyfitResult *result);

#endif /* POLYFIT_H */
