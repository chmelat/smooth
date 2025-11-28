/*  Polynomial fitting for data smoothing
 *  Header file for least squares polynomial approximation
 *  V2.1/2025-11-28/ Updated documentation for V3.1 SVD implementation
 *  V2.0/2025-11-28/ Updated documentation for V3.0 implementation
 *  V1.1/2025-11-28/ Removed unused coeffs field
 *  V1.0/2025-05-27/ Extracted from smooth.c
 */

#ifndef POLYFIT_H
#define POLYFIT_H

/* Structure for polynomial fitting results */
typedef struct {
    double *y_smooth;     /* Smoothed values */
    double *y_deriv;      /* First derivatives */
    int n;                /* Number of points */
    int poly_degree;      /* Polynomial degree requested */
    int window_size;      /* Window size used */
} PolyfitResult;

/* Main polynomial fitting function
 * 
 * Performs local polynomial fitting using least squares method
 * in a sliding window. Uses SVD decomposition for maximum numerical stability.
 *
 * Parameters:
 *   x            - Array of x-coordinates (must be strictly monotonic increasing)
 *   y            - Array of y-values to be smoothed
 *   n            - Number of data points
 *   window_size  - Size of sliding window (must be odd, >= 3)
 *   poly_degree  - Degree of approximating polynomial (0 to DPMAX, default 12)
 * 
 * Returns:
 *   Pointer to PolyfitResult structure containing smoothed values and derivatives
 *   Returns NULL on error.
 * 
 * Notes:
 *   - Memory must be freed using free_polyfit_result()
 *   - Boundary points use polynomial extrapolation from nearest interior fit
 *   - Uses LAPACK dgelss (SVD decomposition) for solving least squares
 *   - Singular values below rcond * s_max are treated as zero (regularization)
 *   - Reports condition number and effective rank for diagnostics
 * 
 * Algorithm:
 *   For each interior point i, fits polynomial p(x) = sum_{k=0}^{d} c_k (x - x_i)^k
 *   by minimizing sum_{j in window} (y_j - p(x_j))^2 using SVD decomposition
 *   of the Vandermonde matrix. The smoothed value is c_0, derivative is c_1.
 *   
 *   SVD provides implicit regularization: if the Vandermonde matrix is 
 *   ill-conditioned, small singular values are truncated, yielding a stable
 *   solution with minimal norm.
 */
PolyfitResult* polyfit_smooth(double *x, double *y, int n, int window_size, int poly_degree);

/* Free allocated memory for PolyfitResult structure
 * 
 * Parameters:
 *   result - Pointer to PolyfitResult structure to be freed
 * 
 * Notes:
 *   - Safe to call with NULL pointer
 */
void free_polyfit_result(PolyfitResult *result);

#endif /* POLYFIT_H */
