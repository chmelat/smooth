/*  Butterworth filter for data smoothing
 *  Implementation of 4th-order digital Butterworth low-pass filter with filtfilt
 *  V1.1/2025-11-03/ Updated grid uniformity check to use CV (consistent with savgol)
 *  V1.0/2025-11-03/ Initial implementation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "butterworth.h"
#include "grid_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Threshold for grid uniformity (same as Savitzky-Golay) */
#define UNIFORMITY_CV_THRESHOLD 0.05

/* Internal function prototypes */
static void butterworth_coefficients(double fc, double *b, double *a);
static void apply_iir_filter(double *b, double *a, double *x, double *y, int n,
                             double *zi, double *zf);
static void compute_initial_conditions(double *b, double *a, double *zi);
static double* pad_signal(double *y, int n, int pad_len, int *total_len);

/* Design 4th-order Butterworth low-pass filter
 *
 * Uses bilinear transformation to convert analog Butterworth prototype
 * to digital IIR filter.
 *
 * Butterworth filter has maximally flat magnitude response:
 * |H(ω)|² = 1 / (1 + (ω/ωc)^(2N))
 * where N = order (4 in our case), ωc = cutoff frequency
 *
 * Parameters:
 *   fc - Normalized cutoff frequency (0 < fc < 0.5)
 *   b  - Output numerator coefficients (5 elements)
 *   a  - Output denominator coefficients (5 elements, a[0]=1)
 */
static void butterworth_coefficients(double fc, double *b, double *a)
{
    /* Design 4th-order Butterworth filter using complex pole calculation
     * Following scipy.signal.butter algorithm
     */

    int N = 4;  /* Filter order */

    /* Step 1: Calculate analog prototype poles on unit circle
     * Butterworth poles: exp(j*θ_k) where θ_k = π/2 + π*(2k+1)/(2N)
     * k = 0, 1, 2, 3 for N=4
     */
    double complex s_poles[4];
    for (int k = 0; k < N; k++) {
        double theta = M_PI/2.0 + M_PI*(2*k + 1)/(2.0*N);
        s_poles[k] = cexp(I * theta);  /* exp(j*θ) = cos(θ) + j*sin(θ) */
    }

    /* Step 2: Prewarp the cutoff frequency for bilinear transform */
    double wc = tan(M_PI * fc);

    /* Step 3: Scale poles by cutoff frequency */
    for (int k = 0; k < N; k++) {
        s_poles[k] *= wc;
    }

    /* Step 4: Apply bilinear transform to get z-domain poles
     * Bilinear transform: s = 2*(z-1)/(z+1)
     * Solving for z: z = (2+s)/(2-s)
     */
    double complex z_poles[4];
    for (int k = 0; k < N; k++) {
        z_poles[k] = (2.0 + s_poles[k]) / (2.0 - s_poles[k]);
    }

    /* Step 5: Form biquad sections from conjugate pole pairs
     * Pair 0: poles k=0 and k=3 (conjugates)
     * Pair 1: poles k=1 and k=2 (conjugates)
     */

    /* Biquad 1: poles z_poles[0] and z_poles[3]
     * Denominator: (z - z_pole[0])*(z - z_pole[3]) = z² - (p0+p3)*z + p0*p3
     * In standard form: 1 + a1*z^-1 + a2*z^-2
     */
    double complex p_sum1 = z_poles[0] + z_poles[3];
    double complex p_prod1 = z_poles[0] * z_poles[3];
    double a1_1 = -creal(p_sum1);   /* Real part */
    double a1_2 = creal(p_prod1);    /* Real part */

    /* Biquad 2: poles z_poles[1] and z_poles[2] */
    double complex p_sum2 = z_poles[1] + z_poles[2];
    double complex p_prod2 = z_poles[1] * z_poles[2];
    double a2_1 = -creal(p_sum2);
    double a2_2 = creal(p_prod2);

    /* Step 6: Calculate gain to ensure unity DC gain
     * For lowpass, numerator is (1 + z^-1)^2 = 1 + 2*z^-1 + z^-2
     * Gain is chosen so that H(z=1) = 1
     */
    double gain1 = (1.0 + a1_1 + a1_2);  /* Denominator at z=1 */
    double gain2 = (1.0 + a2_1 + a2_2);

    /* Numerator coefficients for each biquad (normalized) */
    double b1_0 = gain1 / 4.0;  /* (1+z^-1)^2 has coefficients [1,2,1], sum=4 */
    double b1_1 = 2.0 * gain1 / 4.0;
    double b1_2 = gain1 / 4.0;

    double b2_0 = gain2 / 4.0;
    double b2_1 = 2.0 * gain2 / 4.0;
    double b2_2 = gain2 / 4.0;

    /* Step 7: Cascade the two biquad sections
     * H(z) = H1(z) * H2(z)
     */

    /* Numerator: multiply polynomials */
    b[0] = b1_0 * b2_0;
    b[1] = b1_0*b2_1 + b1_1*b2_0;
    b[2] = b1_0*b2_2 + b1_1*b2_1 + b1_2*b2_0;
    b[3] = b1_1*b2_2 + b1_2*b2_1;
    b[4] = b1_2 * b2_2;

    /* Denominator: multiply polynomials */
    a[0] = 1.0;
    a[1] = a1_1 + a2_1;
    a[2] = a1_2 + a1_1*a2_1 + a2_2;
    a[3] = a1_2*a2_1 + a1_1*a2_2;
    a[4] = a1_2 * a2_2;
}

/* Compute initial conditions for IIR filter
 *
 * Implements scipy.signal.lfilter_zi algorithm.
 * Computes initial state for step response so that y[0] = x[0] for constant input.
 *
 * For Direct Form II Transposed, we solve:
 * (I - A) * zi = B
 * where A and B are derived from filter coefficients.
 *
 * Parameters:
 *   b  - Numerator coefficients (5 elements)
 *   a  - Denominator coefficients (5 elements, a[0]=1)
 *   zi - Output: initial conditions (4 elements)
 */
static void compute_initial_conditions(double *b, double *a, double *zi)
{
    /* For a 4th order filter (5 coefficients), we have 4 state variables.
     * We need to solve a 4x4 linear system.
     *
     * The system is: (I - A) * zi = B
     * where:
     * A[i,j] comes from the filter recursion
     * B[i] = b[i+1] - a[i+1]*b[0]
     *
     * For transposed Direct Form II:
     * zi[0] = b[1] - a[1]*b[0] + sum(...)
     * zi[1] = b[2] - a[2]*b[0] + sum(...)
     * etc.
     *
     * Simplified solution for our specific case (order 4):
     * We use the fact that for constant input, zi should give constant output.
     */

    /* Build the (I - A) matrix and B vector */
    /* For order 4, the matrix is 4x4 */
    double IminusA[4][4];
    double B[4];

    /* Initialize to zero */
    for (int i = 0; i < 4; i++) {
        B[i] = 0.0;
        for (int j = 0; j < 4; j++) {
            IminusA[i][j] = 0.0;
        }
    }

    /* Diagonal: I */
    for (int i = 0; i < 4; i++) {
        IminusA[i][i] = 1.0;
    }

    /* Build the system based on filter structure */
    /* This follows scipy's lfilter_zi implementation exactly */

    /* B vector: B = b[1:] - a[1:] * b[0] */
    for (int i = 0; i < 4; i++) {
        B[i] = b[i+1] - a[i+1] * b[0];
    }

    /* IminusA = I - companion(a).T
     * Companion matrix for [1, a1, a2, a3, a4]:
     *   [a1, -a2, a3, -a4]
     *   [1,   0,   0,   0 ]
     *   [0,   1,   0,   0 ]
     *   [0,   0,   1,   0 ]
     *
     * Transposed:
     *   [a1,  1,  0,  0]
     *   [-a2, 0,  1,  0]
     *   [a3,  0,  0,  1]
     *   [-a4, 0,  0,  0]
     *
     * I - companion(a).T:
     *   [1-a1,  -1,   0,   0]
     *   [a2,     1,  -1,   0]
     *   [-a3,    0,   1,  -1]
     *   [a4,     0,   0,   1]
     */

    IminusA[0][0] = 1.0 - a[1];
    IminusA[0][1] = -1.0;
    IminusA[0][2] = 0.0;
    IminusA[0][3] = 0.0;

    IminusA[1][0] = a[2];
    IminusA[1][1] = 1.0;
    IminusA[1][2] = -1.0;
    IminusA[1][3] = 0.0;

    IminusA[2][0] = -a[3];
    IminusA[2][1] = 0.0;
    IminusA[2][2] = 1.0;
    IminusA[2][3] = -1.0;

    IminusA[3][0] = a[4];
    IminusA[3][1] = 0.0;
    IminusA[3][2] = 0.0;
    IminusA[3][3] = 1.0;

    /* Solve the system using Gaussian elimination (simple for 4x4) */
    /* Forward elimination */
    for (int k = 0; k < 4; k++) {
        /* Find pivot */
        double pivot = IminusA[k][k];
        if (fabs(pivot) < 1e-10) {
            /* Singular matrix - use zero initial conditions */
            for (int i = 0; i < 4; i++) zi[i] = 0.0;
            return;
        }

        /* Eliminate column k */
        for (int i = k+1; i < 4; i++) {
            double factor = IminusA[i][k] / pivot;
            for (int j = k; j < 4; j++) {
                IminusA[i][j] -= factor * IminusA[k][j];
            }
            B[i] -= factor * B[k];
        }
    }

    /* Back substitution */
    for (int i = 3; i >= 0; i--) {
        zi[i] = B[i];
        for (int j = i+1; j < 4; j++) {
            zi[i] -= IminusA[i][j] * zi[j];
        }
        zi[i] /= IminusA[i][i];
    }
}

/* Apply IIR filter using transposed Direct Form II
 *
 * This is the core filtering function implementing the difference equation:
 * y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[4]*x[n-4]
 *                  - a[1]*y[n-1] - ... - a[4]*y[n-4]
 *
 * Parameters:
 *   b  - Numerator coefficients (5 elements)
 *   a  - Denominator coefficients (5 elements, a[0]=1)
 *   x  - Input signal
 *   y  - Output signal (pre-allocated)
 *   n  - Signal length
 *   zi - Initial conditions (4 elements), NULL for zero
 *   zf - Final conditions (4 elements), NULL if not needed
 */
static void apply_iir_filter(double *b, double *a, double *x, double *y, int n,
                             double *zi, double *zf)
{
    double z[4] = {0.0, 0.0, 0.0, 0.0};  /* State variables */

    /* Initialize state if provided */
    if (zi != NULL) {
        z[0] = zi[0];
        z[1] = zi[1];
        z[2] = zi[2];
        z[3] = zi[3];
    }

    /* Apply filter using transposed Direct Form II
     * This form is numerically more stable than Direct Form I
     */
    for (int i = 0; i < n; i++) {
        /* Output = b[0]*input + state[0] */
        y[i] = b[0] * x[i] + z[0];

        /* Update states */
        z[0] = b[1]*x[i] - a[1]*y[i] + z[1];
        z[1] = b[2]*x[i] - a[2]*y[i] + z[2];
        z[2] = b[3]*x[i] - a[3]*y[i] + z[3];
        z[3] = b[4]*x[i] - a[4]*y[i];
    }

    /* Return final state if requested */
    if (zf != NULL) {
        zf[0] = z[0];
        zf[1] = z[1];
        zf[2] = z[2];
        zf[3] = z[3];
    }
}

/* Pad signal using reflection for edge handling
 *
 * Creates padded signal: [y_reflected | y | y_reflected]
 * This minimizes edge effects in filtfilt.
 *
 * Parameters:
 *   y        - Input signal
 *   n        - Signal length
 *   pad_len  - Length of padding on each side
 *   total_len- Output: total length of padded signal
 *
 * Returns:
 *   Pointer to padded signal (must be freed by caller)
 */
static double* pad_signal(double *y, int n, int pad_len, int *total_len)
{
    *total_len = n + 2 * pad_len;
    double *padded = (double*)malloc(*total_len * sizeof(double));

    if (padded == NULL) {
        return NULL;
    }

    /* Left padding: reflect around y[0] */
    for (int i = 0; i < pad_len; i++) {
        int idx = pad_len - 1 - i;  /* index into original array */
        if (idx >= n) idx = n - 1;
        padded[i] = 2.0 * y[0] - y[idx];
    }

    /* Copy original signal */
    memcpy(padded + pad_len, y, n * sizeof(double));

    /* Right padding: reflect around y[n-1] */
    for (int i = 0; i < pad_len; i++) {
        int idx = n - 2 - i;  /* index into original array */
        if (idx < 0) idx = 0;
        padded[pad_len + n + i] = 2.0 * y[n-1] - y[idx];
    }

    return padded;
}

/* Automatic cutoff frequency estimation
 *
 * Simple heuristic: use 10% of Nyquist frequency as default
 * TODO: For better estimation, could analyze power spectrum
 */
double estimate_cutoff_frequency(double *x, double *y, int n)
{
    /* Avoid unused parameter warnings */
    (void)x;
    (void)y;
    (void)n;

    /* Use conservative default: 10% of Nyquist frequency */
    return 0.1;
}

/* Main filtfilt function */
ButterworthResult* butterworth_filtfilt(double *x, double *y, int n,
                                        double cutoff_freq, int auto_cutoff)
{
    /* Input validation */
    if (x == NULL || y == NULL || n < 20) {
        fprintf(stderr, "ERROR: Invalid input (need at least 20 points)\n");
        return NULL;
    }

    /* Automatic cutoff selection if requested */
    double fc = cutoff_freq;
    if (auto_cutoff > 0) {
        fc = estimate_cutoff_frequency(x, y, n);
        printf("# Auto-selected cutoff frequency: fc = %.4f\n", fc);
    }

    /* Validate cutoff frequency */
    if (fc <= 0.0 || fc >= 0.5) {
        fprintf(stderr, "ERROR: Cutoff frequency must be in range (0, 0.5)\n");
        fprintf(stderr, "       Got fc = %.4f\n", fc);
        return NULL;
    }

    /* CRITICAL: Check grid uniformity - Butterworth filter requires uniform grid */
    GridAnalysis *grid = analyze_grid(x, n, 0);
    if (grid == NULL) {
        fprintf(stderr, "Error: Grid analysis failed\n");
        return NULL;
    }

    /* Check if grid is sufficiently uniform for Butterworth filter */
    if (grid->cv > UNIFORMITY_CV_THRESHOLD) {
        fprintf(stderr, "\n");
        fprintf(stderr, "========================================\n");
        fprintf(stderr, "ERROR: Butterworth filter not suitable for non-uniform grid!\n");
        fprintf(stderr, "========================================\n");
        fprintf(stderr, "Grid analysis:\n");
        fprintf(stderr, "  Coefficient of variation (CV) = %.4f\n", grid->cv);
        fprintf(stderr, "  Threshold for uniformity = %.4f\n", UNIFORMITY_CV_THRESHOLD);
        fprintf(stderr, "  h_min = %.6e, h_max = %.6e, h_avg = %.6e\n",
                grid->h_min, grid->h_max, grid->h_avg);
        fprintf(stderr, "  Ratio h_max/h_min = %.2f\n", grid->ratio_max_min);
        fprintf(stderr, "\n");
        fprintf(stderr, "The Butterworth filter assumes uniformly spaced data points.\n");
        fprintf(stderr, "Your data has significant spacing variation.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "RECOMMENDED ALTERNATIVES:\n");
        fprintf(stderr, "  1. Use Tikhonov method: -m 2 -l auto\n");
        fprintf(stderr, "     (Works correctly with non-uniform grids)\n");
        fprintf(stderr, "  2. Resample your data to uniform grid before filtering\n");
        fprintf(stderr, "\n");

        free_grid_analysis(grid);
        return NULL;
    }

    /* Grid is sufficiently uniform - proceed with filtering */
    if (grid->cv > 0.01) {
        /* Grid is uniform enough but not perfectly uniform - just warn */
        printf("# Butterworth: Grid is nearly uniform (CV=%.4f)\n", grid->cv);
    }

    double sample_rate = 1.0 / grid->h_avg;
    free_grid_analysis(grid);

    /* Allocate result structure */
    ButterworthResult *result = (ButterworthResult*)malloc(sizeof(ButterworthResult));
    if (result == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        return NULL;
    }

    result->y_smooth = (double*)malloc(n * sizeof(double));
    if (result->y_smooth == NULL) {
        free(result);
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        return NULL;
    }

    result->n = n;
    result->order = 4;
    result->cutoff_freq = fc;
    result->sample_rate = sample_rate;

    /* Design Butterworth filter */
    double b[5], a[5];
    butterworth_coefficients(fc, b, a);

    /* Padding length: 3 times the filter order is typical */
    int pad_len = 3 * result->order;
    if (pad_len > n / 2) {
        pad_len = n / 2;
    }

    /* Pad the signal */
    int padded_len;
    double *y_padded = pad_signal(y, n, pad_len, &padded_len);
    if (y_padded == NULL) {
        free(result->y_smooth);
        free(result);
        fprintf(stderr, "ERROR: Signal padding failed\n");
        return NULL;
    }

    /* Allocate temporary arrays for filtering */
    double *y_fwd = (double*)malloc(padded_len * sizeof(double));
    double *y_bwd = (double*)malloc(padded_len * sizeof(double));
    double *y_temp = (double*)malloc(padded_len * sizeof(double));

    if (y_fwd == NULL || y_bwd == NULL || y_temp == NULL) {
        free(y_padded);
        free(y_fwd);
        free(y_bwd);
        free(y_temp);
        free(result->y_smooth);
        free(result);
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        return NULL;
    }

    /* Compute initial conditions once using scipy algorithm */
    double zi_base[4];
    double zi[4];
    compute_initial_conditions(b, a, zi_base);

    /* Forward filtering: scale initial conditions by first value */
    for (int i = 0; i < 4; i++) {
        zi[i] = zi_base[i] * y_padded[0];
    }
    apply_iir_filter(b, a, y_padded, y_fwd, padded_len, zi, NULL);

    /* Reverse the signal for backward pass */
    for (int i = 0; i < padded_len; i++) {
        y_temp[i] = y_fwd[padded_len - 1 - i];
    }

    /* Backward filtering: scale initial conditions by first value of reversed signal */
    for (int i = 0; i < 4; i++) {
        zi[i] = zi_base[i] * y_temp[0];
    }
    apply_iir_filter(b, a, y_temp, y_bwd, padded_len, zi, NULL);

    /* Reverse back and extract original signal (remove padding) */
    for (int i = 0; i < n; i++) {
        result->y_smooth[i] = y_bwd[padded_len - 1 - pad_len - i];
    }

    /* Cleanup */
    free(y_padded);
    free(y_fwd);
    free(y_bwd);
    free(y_temp);

    return result;
}

/* Free result structure */
void free_butterworth_result(ButterworthResult *result)
{
    if (result != NULL) {
        free(result->y_smooth);
        free(result);
    }
}
