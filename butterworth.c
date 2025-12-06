/* Butterworth filter for data smoothing
 * Implementation of 4th-order digital Butterworth low-pass filter with filtfilt
 * V1.3/2025-12-06/ Memory optimization (2 buffers instead of 4), symbolic constants,
 *                  improved validation, const-correctness
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
#define NEARLY_UNIFORM_CV_THRESHOLD 0.01

/* Cutoff frequency limits */
#define CUTOFF_FREQ_MIN 0.0
#define CUTOFF_FREQ_MAX 1.0
#define CUTOFF_FREQ_STABILITY_WARN 0.95

/* Internal function prototypes */
static void butterworth_coefficients(double fc, double *b, double *a);
static void apply_iir_filter(const double *b, const double *a,
                             const double *x, double *y, size_t n,
                             const double *zi, double *zf);
static void compute_initial_conditions(const double *b, const double *a, double *zi);
static double* pad_signal(const double *y, int n, int pad_len, size_t *total_len);
static void reverse_array_inplace(double *arr, size_t n);

/* Calculate padding length (SciPy convention: 3 * max(len(a), len(b)) - 1) */
static inline int calculate_pad_length(int n)
{
    int pad_len = 3 * BUTTERWORTH_NUM_COEFFS - 1;
    if (pad_len >= n) {
        pad_len = n - 1;
    }
    return pad_len;
}

/* Estimate memory usage in bytes */
static inline size_t estimate_memory_usage(int n, int pad_len)
{
    size_t padded_len = (size_t)n + 2 * (size_t)pad_len;
    /* 2 temp buffers + 1 result buffer */
    return (2 * padded_len + (size_t)n) * sizeof(double) + sizeof(ButterworthResult);
}

/* Design 4th-order Butterworth low-pass filter */
static void butterworth_coefficients(double fc, double *b, double *a)
{
    /* Step 1: Calculate analog prototype poles on unit circle */
    double complex s_poles[BUTTERWORTH_ORDER];
    for (int k = 0; k < BUTTERWORTH_ORDER; k++) {
        double theta = M_PI/2.0 + M_PI*(2*k + 1)/(2.0*BUTTERWORTH_ORDER);
        s_poles[k] = cexp(I * theta);
    }

    /* Step 2: Prewarp the cutoff frequency for bilinear transform */
    double wc = tan(M_PI / 2.0 * fc);

    /* Step 3: Scale poles by cutoff frequency */
    for (int k = 0; k < BUTTERWORTH_ORDER; k++) {
        s_poles[k] *= wc;
    }

    /* Step 4: Apply bilinear transform to get z-domain poles */
    double complex z_poles[BUTTERWORTH_ORDER];
    for (int k = 0; k < BUTTERWORTH_ORDER; k++) {
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

    /* Step 6: Calculate normalization factors for unity DC gain
     * For a biquad with numerator [1, 2, 1] and denominator [1, a1, a2]:
     * H(z=1) = (1+2+1)/(1+a1+a2) = 4/(1+a1+a2)
     * To normalize: multiply numerator by (1+a1+a2)/4
     */
    double norm1 = (1.0 + a1_1 + a1_2) / 4.0;
    double norm2 = (1.0 + a2_1 + a2_2) / 4.0;

    /* Numerator coefficients for each biquad (normalized) */
    double b1_0 = norm1;
    double b1_1 = 2.0 * norm1;
    double b1_2 = norm1;

    double b2_0 = norm2;
    double b2_1 = 2.0 * norm2;
    double b2_2 = norm2;

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
static void compute_initial_conditions(const double *b, const double *a, double *zi)
{
    /* Build the (I - A) matrix and B vector */
    double IminusA[BUTTERWORTH_ORDER][BUTTERWORTH_ORDER];
    double B[BUTTERWORTH_ORDER];

    /* Initialize to zero */
    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
        B[i] = 0.0;
        for (int j = 0; j < BUTTERWORTH_ORDER; j++) {
            IminusA[i][j] = 0.0;
        }
    }

    /* Diagonal: I */
    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
        IminusA[i][i] = 1.0;
    }

    /* B vector: B = b[1:] - a[1:] * b[0] */
    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
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
    double A_colmajor[BUTTERWORTH_ORDER * BUTTERWORTH_ORDER];
    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
        for (int j = 0; j < BUTTERWORTH_ORDER; j++) {
            A_colmajor[j*BUTTERWORTH_ORDER + i] = IminusA[i][j];
        }
    }

    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
        zi[i] = B[i];
    }

    int n_sys = BUTTERWORTH_ORDER;
    int nrhs = 1;
    int lda = BUTTERWORTH_ORDER;
    int ipiv[BUTTERWORTH_ORDER];
    int ldb = BUTTERWORTH_ORDER;
    int info;

    extern void dgesv_(int *n, int *nrhs, double *A, int *lda,
                       int *ipiv, double *B, int *ldb, int *info);

    dgesv_(&n_sys, &nrhs, A_colmajor, &lda, ipiv, zi, &ldb, &info);

    if (info != 0) {
        fprintf(stderr, "Warning: LAPACK dgesv failed (info=%d), using zero initial conditions\n", info);
        for (int i = 0; i < BUTTERWORTH_ORDER; i++) zi[i] = 0.0;
    }
}

/* Apply IIR filter using transposed Direct Form II */
static void apply_iir_filter(const double *b, const double *a,
                             const double *x, double *y, size_t n,
                             const double *zi, double *zf)
{
    double z[BUTTERWORTH_ORDER] = {0.0, 0.0, 0.0, 0.0};

    if (zi != NULL) {
        for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
            z[i] = zi[i];
        }
    }

    for (size_t i = 0; i < n; i++) {
        y[i] = b[0] * x[i] + z[0];

        z[0] = b[1]*x[i] - a[1]*y[i] + z[1];
        z[1] = b[2]*x[i] - a[2]*y[i] + z[2];
        z[2] = b[3]*x[i] - a[3]*y[i] + z[3];
        z[3] = b[4]*x[i] - a[4]*y[i];
    }

    if (zf != NULL) {
        for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
            zf[i] = z[i];
        }
    }
}

/* Reverse array in-place */
static void reverse_array_inplace(double *arr, size_t n)
{
    size_t i = 0;
    size_t j = n - 1;
    while (i < j) {
        double tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
        i++;
        j--;
    }
}

/* Pad signal using odd reflection */
static double* pad_signal(const double *y, int n, int pad_len, size_t *total_len)
{
    *total_len = (size_t)n + 2 * (size_t)pad_len;
    double *padded = (double*)malloc(*total_len * sizeof(double));

    if (padded == NULL) {
        return NULL;
    }

    memcpy(padded + pad_len, y, (size_t)n * sizeof(double));

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
double estimate_cutoff_frequency(const double *x, const double *y, int n)
{
    (void)x; (void)y; (void)n;
    /* TODO: Implement actual estimation (e.g., zero-crossing rate, spectral analysis) */
    return 0.2;
}

/* Main filtfilt function with optimized memory usage */
ButterworthResult* butterworth_filtfilt(const double *x, const double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        const GridAnalysis *grid_info)
{
    /* Initialize pointers to NULL for safe cleanup */
    ButterworthResult *result = NULL;
    double *y_padded = NULL;
    double *y_work = NULL;

    /* Input validation */
    if (x == NULL || y == NULL) {
        fprintf(stderr, "ERROR: NULL input pointer\n");
        return NULL;
    }

    if (n < BUTTERWORTH_MIN_POINTS) {
        fprintf(stderr, "ERROR: Need at least %d points (got %d)\n",
                BUTTERWORTH_MIN_POINTS, n);
        return NULL;
    }

    if (n > BUTTERWORTH_MAX_POINTS_LIMIT) {
        fprintf(stderr, "ERROR: Dataset too large (%d points, limit is %d)\n",
                n, BUTTERWORTH_MAX_POINTS_LIMIT);
        return NULL;
    }

    if (grid_info == NULL) {
        fprintf(stderr, "ERROR: Grid info not available\n");
        return NULL;
    }

    /* Calculate padding and check memory requirements */
    int pad_len = calculate_pad_length(n);
    size_t mem_estimate = estimate_memory_usage(n, pad_len);

    if (n > BUTTERWORTH_MAX_POINTS_WARNING) {
        fprintf(stderr, "Warning: Large dataset (%d points) requires ~%.1f GB RAM\n",
                n, (double)mem_estimate / (1024.0 * 1024.0 * 1024.0));
    }

    /* Automatic cutoff selection */
    double fc = cutoff_freq;
    if (auto_cutoff > 0) {
        fc = estimate_cutoff_frequency(x, y, n);
        printf("# Auto-selected cutoff frequency: fc = %.4f\n", fc);
    }

    /* Validate cutoff frequency */
    if (fc <= CUTOFF_FREQ_MIN || fc >= CUTOFF_FREQ_MAX) {
        fprintf(stderr, "ERROR: Cutoff frequency must be in range (%.1f, %.1f), got %.4f\n",
                CUTOFF_FREQ_MIN, CUTOFF_FREQ_MAX, fc);
        return NULL;
    }

    if (fc > CUTOFF_FREQ_STABILITY_WARN) {
        fprintf(stderr, "Warning: fc = %.4f is close to Nyquist limit, "
                "may cause numerical instability\n", fc);
    }

    /* Check grid uniformity */
    if (grid_info->cv > UNIFORMITY_CV_THRESHOLD) {
        fprintf(stderr, "ERROR: Butterworth filter not suitable for "
                "non-uniform grid (CV=%.4f, threshold=%.4f)!\n",
                grid_info->cv, UNIFORMITY_CV_THRESHOLD);
        return NULL;
    }

    if (grid_info->cv > NEARLY_UNIFORM_CV_THRESHOLD) {
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
    result->y_smooth = NULL;

    /* 2. Allocate result data array */
    result->y_smooth = (double*)malloc((size_t)n * sizeof(double));
    if (result->y_smooth == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (result buffer)\n");
        goto error;
    }

    /* Fill metadata */
    result->n = n;
    result->order = BUTTERWORTH_ORDER;
    result->cutoff_freq = fc;
    result->sample_rate = sample_rate;

    /* Design filter */
    double b[BUTTERWORTH_NUM_COEFFS], a[BUTTERWORTH_NUM_COEFFS];
    butterworth_coefficients(fc, b, a);

    /* 3. Allocate padded signal */
    size_t padded_len;
    y_padded = pad_signal(y, n, pad_len, &padded_len);
    if (y_padded == NULL) {
        fprintf(stderr, "ERROR: Signal padding failed\n");
        goto error;
    }

    /* 4. Allocate single work buffer (optimization: reduced from 3 to 1) */
    y_work = (double*)malloc(padded_len * sizeof(double));
    if (y_work == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (work buffer)\n");
        goto error;
    }

    /* --- PROCESSING --- */

    /* Initial conditions */
    double zi_base[BUTTERWORTH_ORDER];
    double zi[BUTTERWORTH_ORDER];
    compute_initial_conditions(b, a, zi_base);

    /* Forward filtering: y_padded -> y_work */
    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
        zi[i] = zi_base[i] * y_padded[0];
    }
    apply_iir_filter(b, a, y_padded, y_work, padded_len, zi, NULL);

    /* Reverse y_work in-place */
    reverse_array_inplace(y_work, padded_len);

    /* Backward filtering: y_work -> y_padded (reuse as output buffer) */
    for (int i = 0; i < BUTTERWORTH_ORDER; i++) {
        zi[i] = zi_base[i] * y_work[0];
    }
    apply_iir_filter(b, a, y_work, y_padded, padded_len, zi, NULL);

    /* Reverse y_padded in-place */
    reverse_array_inplace(y_padded, padded_len);

    /* Extract result (from reversed padded buffer) */
    for (int i = 0; i < n; i++) {
        result->y_smooth[i] = y_padded[pad_len + i];
    }

    /* --- CLEANUP (Success) --- */
    free(y_padded);
    free(y_work);

    return result;

    /* --- ERROR HANDLER --- */
error:
    if (y_padded) free(y_padded);
    if (y_work) free(y_work);
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
