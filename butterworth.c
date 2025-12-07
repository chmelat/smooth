/* Butterworth filter for data smoothing
 * V1.4/2025-12-07/ Biquad cascade implementation with proper initial conditions
 *                  Combines numerical stability of biquads with robust IC computation
 * V1.3/2025-12-06/ Memory optimization, symbolic constants, improved validation
 * V1.2/2025-11-21/ Refactored error handling (goto pattern) for memory safety
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "butterworth.h"
#include "grid_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Thresholds */
#define UNIFORMITY_CV_THRESHOLD 0.05
#define NEARLY_UNIFORM_CV_THRESHOLD 0.01
#define CUTOFF_FREQ_MIN 0.0
#define CUTOFF_FREQ_MAX 1.0
#define CUTOFF_FREQ_STABILITY_WARN 0.95

/* Internal function prototypes */
static void design_biquad_sections(double fc, BiquadSection *sections);
static void compute_biquad_ic(const BiquadSection *bq, double *zi_base);
static void apply_biquad(const BiquadSection *bq, const double *x, double *y, 
                         size_t n, double z[2]);
static double* pad_signal(const double *y, int n, int pad_len, size_t *total_len);
static void reverse_array_inplace(double *arr, size_t n);

/* Calculate padding length (3 * filter_order for biquad cascade) */
static inline int calculate_pad_length(int n)
{
    int pad_len = 3 * (BUTTERWORTH_ORDER + 1) - 1;
    if (pad_len >= n) {
        pad_len = n - 1;
    }
    return pad_len;
}

/* Estimate memory usage in bytes */
static inline size_t estimate_memory_usage(int n, int pad_len)
{
    size_t padded_len = (size_t)n + 2 * (size_t)pad_len;
    return (2 * padded_len + (size_t)n) * sizeof(double) + sizeof(ButterworthResult);
}

/* Design 4th-order Butterworth as 2 cascaded biquad sections
 * Uses bilinear transform with prewarping
 */
static void design_biquad_sections(double fc, BiquadSection *sections)
{
    /* Prewarp cutoff frequency */
    double Wc = tan(M_PI * fc / 2.0);
    double Wc_sq = Wc * Wc;

    /* 4th-order Butterworth pole angles (2 conjugate pairs) */
    /* Poles at angles: π/8, 3π/8, 5π/8, 7π/8 from positive real axis */
    /* We use the two unique angles for the two biquad sections */
    double theta[NUM_BIQUADS] = {
        M_PI * (1.0 / 8.0),   /* First conjugate pair */
        M_PI * (3.0 / 8.0)    /* Second conjugate pair */
    };
    
    for (int i = 0; i < NUM_BIQUADS; i++) {
        /* Analog prototype: H(s) = 1 / (s² + 2·cos(θ)·s + 1) */
        /* Scaled by Wc: H(s) = Wc² / (s² + 2·Wc·cos(θ)·s + Wc²) */
        double cos_theta = cos(theta[i]);
        double alpha = 2.0 * Wc * cos_theta;  /* 2·Wc·cos(θ) */

        /* Bilinear transform: s = 2·(z-1)/(z+1)
         * Denominator coefficients after transform:
         * A(z) = (4 + 2α + Wc²) + (2Wc² - 8)z⁻¹ + (4 - 2α + Wc²)z⁻²
         */
        double A0 = 4.0 + 2.0 * alpha + Wc_sq;  /* Normalization factor */
        
        /* Normalized denominator: a[0] = 1 */
        sections[i].a[0] = 1.0;
        sections[i].a[1] = (2.0 * Wc_sq - 8.0) / A0;
        sections[i].a[2] = (4.0 - 2.0 * alpha + Wc_sq) / A0;

        /* Numerator: B(z) = Wc²·(1 + 2z⁻¹ + z⁻²) for lowpass
         * Normalized for unity DC gain: H(z=1) = 1
         */
        double gain = Wc_sq / A0;
        sections[i].b[0] = gain;
        sections[i].b[1] = 2.0 * gain;
        sections[i].b[2] = gain;
    }
}

/* Compute initial conditions for one biquad section
 * Solves (I - A)·zi = B for steady-state response to step input
 * where A is the companion matrix of the denominator
 *
 * For TDF-II biquad with a[0]=1:
 *   (I - A) = [[1+a1, -1], [a2, 1]]
 *   B = [b1 - a1·b0, b2 - a2·b0]
 */
static void compute_biquad_ic(const BiquadSection *bq, double *zi_base)
{
    double b0 = bq->b[0], b1 = bq->b[1], b2 = bq->b[2];
    double a1 = bq->a[1], a2 = bq->a[2];
    
    /* Determinant of (I - A): det = (1+a1)·1 - (-1)·a2 = 1 + a1 + a2 */
    double det = 1.0 + a1 + a2;
    
    /* B vector */
    double B0 = b1 - a1 * b0;
    double B1 = b2 - a2 * b0;
    
    /* Solve using Cramer's rule / direct inverse of 2x2 matrix
     * (I-A)^(-1) = (1/det) · [[1, 1], [-a2, 1+a1]]
     */
    if (fabs(det) > 1e-10) {
        zi_base[0] = (B0 + B1) / det;
        zi_base[1] = (-a2 * B0 + (1.0 + a1) * B1) / det;
    } else {
        /* Fallback for degenerate case (should not happen with valid fc) */
        zi_base[0] = 0.0;
        zi_base[1] = 0.0;
    }
}

/* Apply one biquad section using Transposed Direct Form II
 * This form is numerically stable and efficient
 *
 * y[n] = b0·x[n] + z[0]
 * z[0] = b1·x[n] - a1·y[n] + z[1]
 * z[1] = b2·x[n] - a2·y[n]
 */
static void apply_biquad(const BiquadSection *bq, const double *x, double *y, 
                         size_t n, double z[2])
{
    double b0 = bq->b[0], b1 = bq->b[1], b2 = bq->b[2];
    double a1 = bq->a[1], a2 = bq->a[2];
    
    for (size_t i = 0; i < n; i++) {
        double xi = x[i];
        double yi = b0 * xi + z[0];
        
        z[0] = b1 * xi - a1 * yi + z[1];
        z[1] = b2 * xi - a2 * yi;
        
        y[i] = yi;
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

/* Pad signal using odd reflection (anti-symmetric extension)
 * This minimizes edge discontinuities for filtfilt
 */
static double* pad_signal(const double *y, int n, int pad_len, size_t *total_len)
{
    *total_len = (size_t)n + 2 * (size_t)pad_len;
    double *padded = (double*)malloc(*total_len * sizeof(double));

    if (padded == NULL) {
        return NULL;
    }

    /* Copy original data to center */
    memcpy(padded + pad_len, y, (size_t)n * sizeof(double));

    /* Left padding: odd reflection around y[0]
     * padded[pad_len-1-i] = 2·y[0] - y[1+i]
     */
    double left_val = y[0];
    for (int i = 0; i < pad_len; i++) {
        int src_idx = pad_len - i;
        if (src_idx >= n) src_idx = n - 1;
        padded[i] = 2.0 * left_val - y[src_idx];
    }

    /* Right padding: odd reflection around y[n-1]
     * padded[pad_len+n+i] = 2·y[n-1] - y[n-2-i]
     */
    double right_val = y[n-1];
    for (int i = 0; i < pad_len; i++) {
        int src_idx = n - 2 - i;
        if (src_idx < 0) src_idx = 0;
        padded[pad_len + n + i] = 2.0 * right_val - y[src_idx];
    }

    return padded;
}

/* Automatic cutoff frequency estimation (placeholder) */
double estimate_cutoff_frequency(const double *x, const double *y, int n)
{
    (void)x; (void)y; (void)n;
    return 0.2;
}

/* Main filtfilt function using biquad cascade with proper IC */
ButterworthResult* butterworth_filtfilt(const double *x, const double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        const GridAnalysis *grid_info)
{
    /* Initialize pointers for safe cleanup */
    ButterworthResult *result = NULL;
    double *y_padded = NULL;
    double *y_work = NULL;

    /* --- Input validation --- */
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

    /* Calculate padding and memory estimate */
    int pad_len = calculate_pad_length(n);
    size_t mem_estimate = estimate_memory_usage(n, pad_len);

    if (n > BUTTERWORTH_MAX_POINTS_WARNING) {
        fprintf(stderr, "Warning: Large dataset (%d points) requires ~%.1f GB RAM\n",
                n, (double)mem_estimate / (1024.0 * 1024.0 * 1024.0));
    }

    /* Auto cutoff selection */
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
        fprintf(stderr, "Warning: fc = %.4f is close to Nyquist limit\n", fc);
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

    /* --- ALLOCATIONS --- */
    result = (ButterworthResult*)malloc(sizeof(ButterworthResult));
    if (result == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (result structure)\n");
        goto error;
    }
    result->y_smooth = NULL;

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

    /* Design filter (2 biquad sections) */
    BiquadSection sections[NUM_BIQUADS];
    design_biquad_sections(fc, sections);

    /* Compute initial conditions for each biquad */
    double zi_base[NUM_BIQUADS][2];
    for (int i = 0; i < NUM_BIQUADS; i++) {
        compute_biquad_ic(&sections[i], zi_base[i]);
    }

    /* Pad signal */
    size_t padded_len;
    y_padded = pad_signal(y, n, pad_len, &padded_len);
    if (y_padded == NULL) {
        fprintf(stderr, "ERROR: Signal padding failed\n");
        goto error;
    }

    /* Allocate work buffer */
    y_work = (double*)malloc(padded_len * sizeof(double));
    if (y_work == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (work buffer)\n");
        goto error;
    }

    /* --- FORWARD FILTERING --- */
    /* Copy input to work buffer */
    memcpy(y_work, y_padded, padded_len * sizeof(double));
    
    /* Apply biquads in cascade with scaled IC */
    double zi[2];
    double first_val = y_work[0];
    
    for (int s = 0; s < NUM_BIQUADS; s++) {
        /* Scale IC by first sample value */
        zi[0] = zi_base[s][0] * first_val;
        zi[1] = zi_base[s][1] * first_val;
        
        /* Apply biquad (in-place: y_work -> y_work) */
        apply_biquad(&sections[s], y_work, y_work, padded_len, zi);
        
        /* Update first_val for next section (output of this section) */
        first_val = y_work[0];
    }

    /* --- BACKWARD FILTERING --- */
    /* Reverse the forward-filtered signal */
    reverse_array_inplace(y_work, padded_len);
    
    /* Apply biquads again with IC based on reversed signal's first value */
    first_val = y_work[0];
    
    for (int s = 0; s < NUM_BIQUADS; s++) {
        zi[0] = zi_base[s][0] * first_val;
        zi[1] = zi_base[s][1] * first_val;
        
        apply_biquad(&sections[s], y_work, y_work, padded_len, zi);
        
        first_val = y_work[0];
    }

    /* Reverse back to original order */
    reverse_array_inplace(y_work, padded_len);

    /* --- EXTRACT RESULT --- */
    for (int i = 0; i < n; i++) {
        result->y_smooth[i] = y_work[pad_len + i];
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
