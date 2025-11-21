/* Butterworth filter for data smoothing
 * Implementation of 4th-order digital Butterworth low-pass filter with filtfilt
 * V1.2/2025-11-21/ Refactored error handling (goto pattern) for memory safety
 * V1.1/2025-11-03/ Updated grid uniformity check to use CV
 * V1.0/2025-11-03/ Initial implementation
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

/* Design 4th-order Butterworth low-pass filter */
static void butterworth_coefficients(double fc, double *b, double *a)
{
    int N = 4;  /* Filter order */

    /* Step 1: Calculate analog prototype poles on unit circle */
    double complex s_poles[4];
    for (int k = 0; k < N; k++) {
        double theta = M_PI/2.0 + M_PI*(2*k + 1)/(2.0*N);
        s_poles[k] = cexp(I * theta);
    }

    /* Step 2: Prewarp the cutoff frequency for bilinear transform */
    double wc = tan(M_PI / 2.0 * fc);

    /* Step 3: Scale poles by cutoff frequency */
    for (int k = 0; k < N; k++) {
        s_poles[k] *= wc;
    }

    /* Step 4: Apply bilinear transform to get z-domain poles */
    double complex z_poles[4];
    for (int k = 0; k < N; k++) {
        z_poles[k] = (2.0 + s_poles[k]) / (2.0 - s_poles[k]);
    }

    /* Step 5: Form biquad sections from conjugate pole pairs */
    /* Biquad 1: poles z_poles[0] and z_poles[3] */
    double complex p_sum1 = z_poles[0] + z_poles[3];
    double complex p_prod1 = z_poles[0] * z_poles[3];
    double a1_1 = -creal(p_sum1);
    double a1_2 = creal(p_prod1);

    /* Biquad 2: poles z_poles[1] and z_poles[2] */
    double complex p_sum2 = z_poles[1] + z_poles[2];
    double complex p_prod2 = z_poles[1] * z_poles[2];
    double a2_1 = -creal(p_sum2);
    double a2_2 = creal(p_prod2);

    /* Step 6: Calculate gain to ensure unity DC gain */
    double gain1 = (1.0 + a1_1 + a1_2);
    double gain2 = (1.0 + a2_1 + a2_2);

    /* Numerator coefficients for each biquad (normalized) */
    double b1_0 = gain1 / 4.0;
    double b1_1 = 2.0 * gain1 / 4.0;
    double b1_2 = gain1 / 4.0;

    double b2_0 = gain2 / 4.0;
    double b2_1 = 2.0 * gain2 / 4.0;
    double b2_2 = gain2 / 4.0;

    /* Step 7: Cascade the two biquad sections */
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

/* Compute initial conditions for IIR filter using LAPACK */
static void compute_initial_conditions(double *b, double *a, double *zi)
{
    /* Build the (I - A) matrix and B vector */
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

    /* B vector: B = b[1:] - a[1:] * b[0] */
    for (int i = 0; i < 4; i++) {
        B[i] = b[i+1] - a[i+1] * b[0];
    }

    /* IminusA = I - companion(a).T */
    IminusA[0][0] = 1.0 + a[1];
    IminusA[0][1] = -1.0;
    IminusA[0][2] = 0.0;
    IminusA[0][3] = 0.0;

    IminusA[1][0] = a[2];
    IminusA[1][1] = 1.0;
    IminusA[1][2] = -1.0;
    IminusA[1][3] = 0.0;

    IminusA[2][0] = a[3];
    IminusA[2][1] = 0.0;
    IminusA[2][2] = 1.0;
    IminusA[2][3] = -1.0;

    IminusA[3][0] = a[4];
    IminusA[3][1] = 0.0;
    IminusA[3][2] = 0.0;
    IminusA[3][3] = 1.0;

    /* Solve the system using LAPACK dgesv */
    double A_colmajor[16];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            A_colmajor[j*4 + i] = IminusA[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        zi[i] = B[i];
    }

    int n = 4;
    int nrhs = 1;
    int lda = 4;
    int ipiv[4];
    int ldb = 4;
    int info;

    extern void dgesv_(int *n, int *nrhs, double *A, int *lda,
                       int *ipiv, double *B, int *ldb, int *info);

    dgesv_(&n, &nrhs, A_colmajor, &lda, ipiv, zi, &ldb, &info);

    if (info != 0) {
        fprintf(stderr, "Warning: LAPACK dgesv failed (info=%d), using zero initial conditions\n", info);
        for (int i = 0; i < 4; i++) zi[i] = 0.0;
    }
}

/* Apply IIR filter using transposed Direct Form II */
static void apply_iir_filter(double *b, double *a, double *x, double *y, int n,
                             double *zi, double *zf)
{
    double z[4] = {0.0, 0.0, 0.0, 0.0};

    if (zi != NULL) {
        z[0] = zi[0];
        z[1] = zi[1];
        z[2] = zi[2];
        z[3] = zi[3];
    }

    for (int i = 0; i < n; i++) {
        y[i] = b[0] * x[i] + z[0];

        z[0] = b[1]*x[i] - a[1]*y[i] + z[1];
        z[1] = b[2]*x[i] - a[2]*y[i] + z[2];
        z[2] = b[3]*x[i] - a[3]*y[i] + z[3];
        z[3] = b[4]*x[i] - a[4]*y[i];
    }

    if (zf != NULL) {
        zf[0] = z[0];
        zf[1] = z[1];
        zf[2] = z[2];
        zf[3] = z[3];
    }
}

/* Pad signal using odd reflection */
static double* pad_signal(double *y, int n, int pad_len, int *total_len)
{
    *total_len = n + 2 * pad_len;
    double *padded = (double*)malloc(*total_len * sizeof(double));

    if (padded == NULL) {
        return NULL;
    }

    memcpy(padded + pad_len, y, n * sizeof(double));

    double left_end = y[0];
    for (int i = 0; i < pad_len; i++) {
        int src_idx = pad_len - i;
        if (src_idx >= n) src_idx = n - 1;
        padded[i] = 2.0 * left_end - y[src_idx];
    }

    double right_end = y[n-1];
    for (int i = 0; i < pad_len; i++) {
        int src_idx = n - 2 - i;
        if (src_idx < 0) src_idx = 0;
        padded[pad_len + n + i] = 2.0 * right_end - y[src_idx];
    }

    return padded;
}

/* Automatic cutoff frequency estimation */
double estimate_cutoff_frequency(double *x, double *y, int n)
{
    (void)x; (void)y; (void)n;
    /* TODO: Implement actual estimation (e.g., zero-crossing rate) */
    return 0.2;
}

/* Main filtfilt function with improved Error Handling */
ButterworthResult* butterworth_filtfilt(double *x, double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        GridAnalysis *grid_info)
{
    /* Initialize pointers to NULL for safe cleanup */
    ButterworthResult *result = NULL;
    double *y_padded = NULL;
    double *y_fwd = NULL;
    double *y_bwd = NULL;
    double *y_temp = NULL;

    /* Input validation */
    if (x == NULL || y == NULL || n < 20) {
        fprintf(stderr, "ERROR: Invalid input (need at least 20 points)\n");
        return NULL;
    }

    if (grid_info == NULL) {
        fprintf(stderr, "Error: Grid info not available\n");
        return NULL;
    }

    /* Automatic cutoff selection */
    double fc = cutoff_freq;
    if (auto_cutoff > 0) {
        fc = estimate_cutoff_frequency(x, y, n);
        printf("# Auto-selected cutoff frequency: fc = %.4f\n", fc);
    }

    /* Validate cutoff frequency */
    if (fc <= 0.0 || fc >= 1.0) {
        fprintf(stderr, "ERROR: Cutoff frequency must be in range (0, 1)\n");
        return NULL;
    }

    /* Check grid uniformity */
    if (grid_info->cv > UNIFORMITY_CV_THRESHOLD) {
        fprintf(stderr, "ERROR: Butterworth filter not suitable for non-uniform grid (CV=%.4f)!\n", grid_info->cv);
        return NULL;
    }

    if (grid_info->cv > 0.01) {
        printf("# Butterworth: Grid is nearly uniform (CV=%.4f)\n", grid_info->cv);
    }

    double sample_rate = 1.0 / grid_info->h_avg;

    /* --- ALLOCATIONS START --- */

    /* 1. Allocate result structure */
    result = (ButterworthResult*)malloc(sizeof(ButterworthResult));
    if (result == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (result structure)\n");
        goto error;
    }
    /* Important: Initialize internal pointer to NULL immediately */
    result->y_smooth = NULL;

    /* 2. Allocate result data array */
    result->y_smooth = (double*)malloc(n * sizeof(double));
    if (result->y_smooth == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (result buffer)\n");
        goto error;
    }

    /* Fill metadata */
    result->n = n;
    result->order = 4;
    result->cutoff_freq = fc;
    result->sample_rate = sample_rate;

    /* Design filter */
    double b[5], a[5];
    butterworth_coefficients(fc, b, a);

    /* Determine padding */
    int ntaps = 5;
    int pad_len = 3 * ntaps;
    if (pad_len >= n) pad_len = n - 1;

    /* 3. Allocate padded signal */
    int padded_len;
    y_padded = pad_signal(y, n, pad_len, &padded_len);
    if (y_padded == NULL) {
        fprintf(stderr, "ERROR: Signal padding failed\n");
        goto error;
    }

    /* 4. Allocate temporary buffers */
    y_fwd = (double*)malloc(padded_len * sizeof(double));
    y_bwd = (double*)malloc(padded_len * sizeof(double));
    y_temp = (double*)malloc(padded_len * sizeof(double));

    if (y_fwd == NULL || y_bwd == NULL || y_temp == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (temp buffers)\n");
        goto error;
    }

    /* --- PROCESSING --- */

    /* Initial conditions */
    double zi_base[4];
    double zi[4];
    compute_initial_conditions(b, a, zi_base);

    /* Forward filtering */
    for (int i = 0; i < 4; i++) zi[i] = zi_base[i] * y_padded[0];
    apply_iir_filter(b, a, y_padded, y_fwd, padded_len, zi, NULL);

    /* Reverse */
    for (int i = 0; i < padded_len; i++) {
        y_temp[i] = y_fwd[padded_len - 1 - i];
    }

    /* Backward filtering */
    for (int i = 0; i < 4; i++) zi[i] = zi_base[i] * y_temp[0];
    apply_iir_filter(b, a, y_temp, y_bwd, padded_len, zi, NULL);

    /* Final Reverse */
    for (int i = 0; i < padded_len; i++) {
        y_temp[i] = y_bwd[padded_len - 1 - i];
    }

    /* Extract result */
    for (int i = 0; i < n; i++) {
        result->y_smooth[i] = y_temp[pad_len + i];
    }

    /* --- CLEANUP (Success) --- */
    free(y_padded);
    free(y_fwd);
    free(y_bwd);
    free(y_temp);

    return result;

    /* --- ERROR HANDLER --- */
error:
    /* Safe cleanup of all resources */
    if (y_padded) free(y_padded);
    if (y_fwd) free(y_fwd);
    if (y_bwd) free(y_bwd);
    if (y_temp) free(y_temp);
    
    /* free_butterworth_result checks for NULL internally, 
     * and handles result->y_smooth correctly if we initialized it to NULL */
    free_butterworth_result(result);
    
    return NULL;
}

/* Free result structure */
void free_butterworth_result(ButterworthResult *result)
{
    if (result != NULL) {
        if (result->y_smooth != NULL) {
            free(result->y_smooth);
        }
        free(result);
    }
}
