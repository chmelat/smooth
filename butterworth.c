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
#define UNIFORMITY_CV_THRESHOLD 0.15
#define UNIFORMITY_CV_WARNING   0.05
#define NEARLY_UNIFORM_CV_THRESHOLD 0.01
#define CUTOFF_FREQ_MIN 0.0
#define CUTOFF_FREQ_MAX 1.0
#define FC_MIN_PRACTICAL 1e-4  /* Poles approach unit circle below this */
/* Above this threshold the filter is numerically stable (pole radius < 0.75
 * at fc=0.95, far from unit circle) but passes almost the entire spectrum
 * unattenuated. See check_pole_stability() for the actual numerical limit. */
#define CUTOFF_FREQ_INEFFECTIVE_WARN 0.95
#define POLE_RADIUS_WARN   0.99   /* warn when poles approach unit circle */
#define POLE_RADIUS_ERROR  1.0    /* filter marginally stable / unstable */

/* Auto cutoff selection (Morozov's discrepancy principle) */
#define NOISE_MAD_NORMALIZATION  1.6528553  /* sqrt(6) * 0.6745 */
#define DISCREPANCY_TOLERANCE    1.1
#define AUTO_CUTOFF_FALLBACK     0.2
#define N_AUTO_CANDIDATES        6

/* Internal function prototypes */
static void design_biquad_sections(double fc, BiquadSection *sections);
static int  check_pole_stability(const BiquadSection *sections);
static void compute_biquad_ic(const BiquadSection *bq, double *zi_base);
static void apply_biquad(const BiquadSection *bq, const double *x, double *y,
                         size_t n, double z[2]);
static void apply_cascade(const BiquadSection *sections,
                          const double zi_base[][2],
                          double *buf, size_t padded_len);
static double* pad_signal(const double *y, int n, int pad_len, size_t *total_len);
static void reverse_array_inplace(double *arr, size_t n);
static void compute_derivatives_5pt(const double *y_smooth, int n,
                                    double h, double *y_deriv);
static int    compare_double(const void *a, const void *b);
static double estimate_noise_sigma(const double *y, int n);
static double residual_std(const double *y, const double *y_smooth, int n);
static int    run_filtfilt_trial(const double *y, double *out, int n, double fc);

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
    return (padded_len + (size_t)n) * sizeof(double) + sizeof(ButterworthResult);
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

/* Check that filter poles lie safely inside the unit circle.
 * Returns 0 on success, -1 if any pole is on or outside the unit circle.
 * Prints a warning to stderr if any pole exceeds POLE_RADIUS_WARN.
 *
 * For biquad denominator 1 + a1·z^-1 + a2·z^-2, poles are roots of
 * z^2 + a1·z + a2 = 0.
 */
static int check_pole_stability(const BiquadSection *sections)
{
    double max_radius = 0.0;
    int worst_section = 0;

    for (int i = 0; i < NUM_BIQUADS; i++) {
        double a1 = sections[i].a[1];
        double a2 = sections[i].a[2];
        double disc = a1 * a1 - 4.0 * a2;
        double radius;

        if (disc < 0.0) {
            /* Complex conjugate poles: |z| = sqrt(a2) */
            radius = sqrt(fabs(a2));
        } else {
            /* Real poles (shouldn't happen for Butterworth, handle defensively) */
            double sq = sqrt(disc);
            double r1 = fabs((-a1 + sq) * 0.5);
            double r2 = fabs((-a1 - sq) * 0.5);
            radius = (r1 > r2) ? r1 : r2;
        }

        if (radius > max_radius) {
            max_radius = radius;
            worst_section = i;
        }
    }

    if (max_radius >= POLE_RADIUS_ERROR) {
        fprintf(stderr, "ERROR: Filter is unstable: pole radius %.6f >= 1.0 "
                "in biquad section %d\n", max_radius, worst_section);
        return -1;
    }

    if (max_radius > POLE_RADIUS_WARN) {
        fprintf(stderr, "Warning: Filter poles very close to unit circle "
                "(max radius %.6f in biquad %d). Numerical precision may suffer. "
                "Consider a less extreme cutoff frequency.\n",
                max_radius, worst_section);
    }

    return 0;
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

/* Apply biquad cascade in-place with IC scaled by the signal's first value.
 * Used by both forward and backward passes of filtfilt.
 * Each biquad's IC is scaled to match the step amplitude seen at its input
 * (output of previous biquad), relying on unity DC gain of each section.
 */
static void apply_cascade(const BiquadSection *sections,
                          const double zi_base[][2],
                          double *buf, size_t padded_len)
{
    double zi[2];
    double first_val = buf[0];

    for (int s = 0; s < NUM_BIQUADS; s++) {
        zi[0] = zi_base[s][0] * first_val;
        zi[1] = zi_base[s][1] * first_val;
        apply_biquad(&sections[s], buf, buf, padded_len, zi);
        first_val = buf[0];
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

/* Compute first derivatives of filtered signal using 5-point stencils,
 * O(h^4) accuracy with uniform spacing assumption (h = grid->h_avg).
 * Butterworth already requires CV <= 0.15, so using h_avg adds only
 * O(CV*h^2) additional error, well within filter assumptions.
 *
 * Requires n >= 5 (guaranteed by BUTTERWORTH_MIN_POINTS = 20).
 */
static void compute_derivatives_5pt(const double *y_smooth, int n,
                                    double h, double *y_deriv)
{
    const double inv12h = 1.0 / (12.0 * h);

    /* Left boundary (forward 5-point stencils) */
    y_deriv[0] = (-25.0*y_smooth[0] + 48.0*y_smooth[1]
                  - 36.0*y_smooth[2] + 16.0*y_smooth[3]
                  -  3.0*y_smooth[4]) * inv12h;
    y_deriv[1] = (-3.0*y_smooth[0] - 10.0*y_smooth[1]
                  + 18.0*y_smooth[2] -  6.0*y_smooth[3]
                  +       y_smooth[4]) * inv12h;

    /* Interior (central 5-point stencil) */
    for (int i = 2; i < n - 2; i++) {
        y_deriv[i] = (-y_smooth[i+2] + 8.0*y_smooth[i+1]
                      - 8.0*y_smooth[i-1] + y_smooth[i-2]) * inv12h;
    }

    /* Right boundary (backward 5-point stencils) */
    y_deriv[n-2] = (       -y_smooth[n-5] +  6.0*y_smooth[n-4]
                    - 18.0*y_smooth[n-3] + 10.0*y_smooth[n-2]
                    +  3.0*y_smooth[n-1]) * inv12h;
    y_deriv[n-1] = ( 3.0*y_smooth[n-5] - 16.0*y_smooth[n-4]
                    + 36.0*y_smooth[n-3] - 48.0*y_smooth[n-2]
                    + 25.0*y_smooth[n-1]) * inv12h;
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

/* Ascending double comparator for qsort */
static int compare_double(const void *a, const void *b)
{
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

/* Estimate noise standard deviation from second differences of y.
 * Uses MAD of Delta^2 y with MAD-to-sigma normalization for Gaussian noise.
 * For white noise: Var(Delta^2 y) = 6 * sigma^2, MAD = 0.6745 * sigma.
 * Returns -1.0 on failure (n too small, allocation, or zero MAD).
 */
static double estimate_noise_sigma(const double *y, int n)
{
    if (n < 3) return -1.0;

    int m = n - 2;
    double *d2 = (double*)malloc((size_t)m * sizeof(double));
    if (d2 == NULL) return -1.0;

    for (int i = 0; i < m; i++) {
        d2[i] = y[i+2] - 2.0 * y[i+1] + y[i];
    }

    qsort(d2, (size_t)m, sizeof(double), compare_double);
    double median = (m % 2) ? d2[m/2] : 0.5 * (d2[m/2 - 1] + d2[m/2]);

    for (int i = 0; i < m; i++) {
        d2[i] = fabs(d2[i] - median);
    }
    qsort(d2, (size_t)m, sizeof(double), compare_double);
    double mad = (m % 2) ? d2[m/2] : 0.5 * (d2[m/2 - 1] + d2[m/2]);

    free(d2);

    if (mad <= 0.0) return -1.0;
    return mad / NOISE_MAD_NORMALIZATION;
}

/* Compute sample standard deviation of residuals r = y - y_smooth */
static double residual_std(const double *y, const double *y_smooth, int n)
{
    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        double r = y[i] - y_smooth[i];
        sum += r;
        sum_sq += r * r;
    }
    double mean = sum / (double)n;
    double var = sum_sq / (double)n - mean * mean;
    return (var > 0.0) ? sqrt(var) : 0.0;
}

/* Lightweight filtfilt for a single cutoff candidate — no validation, no prints.
 * Writes n smoothed samples to `out`. Returns 0 on success, -1 on alloc failure.
 */
static int run_filtfilt_trial(const double *y, double *out, int n, double fc)
{
    BiquadSection sections[NUM_BIQUADS];
    double zi_base[NUM_BIQUADS][2];

    design_biquad_sections(fc, sections);
    for (int i = 0; i < NUM_BIQUADS; i++) {
        compute_biquad_ic(&sections[i], zi_base[i]);
    }

    int pad_len = calculate_pad_length(n);
    size_t padded_len;
    double *buf = pad_signal(y, n, pad_len, &padded_len);
    if (buf == NULL) return -1;

    apply_cascade(sections, zi_base, buf, padded_len);
    reverse_array_inplace(buf, padded_len);
    apply_cascade(sections, zi_base, buf, padded_len);
    reverse_array_inplace(buf, padded_len);

    memcpy(out, buf + pad_len, (size_t)n * sizeof(double));
    free(buf);
    return 0;
}

/* Automatic cutoff frequency selection via Morozov's discrepancy principle.
 * 1. Estimate noise sigma from MAD of second differences.
 * 2. For increasing fc candidates, find smallest fc where residual std
 *    does not exceed DISCREPANCY_TOLERANCE * sigma_hat (signal preserved).
 * 3. On failure, return AUTO_CUTOFF_FALLBACK.
 */
double estimate_cutoff_frequency(const double *y, int n)
{
    static const double fc_candidates[N_AUTO_CANDIDATES] = {
        0.02, 0.05, 0.1, 0.2, 0.35, 0.5
    };

    double sigma_hat = estimate_noise_sigma(y, n);
    if (sigma_hat <= 0.0) {
        printf("# Auto cutoff: noise estimation failed, using fallback fc=%.2f\n",
               AUTO_CUTOFF_FALLBACK);
        return AUTO_CUTOFF_FALLBACK;
    }

    double *trial = (double*)malloc((size_t)n * sizeof(double));
    if (trial == NULL) {
        printf("# Auto cutoff: allocation failed, using fallback fc=%.2f\n",
               AUTO_CUTOFF_FALLBACK);
        return AUTO_CUTOFF_FALLBACK;
    }

    double selected = AUTO_CUTOFF_FALLBACK;
    double selected_res = -1.0;

    for (int k = 0; k < N_AUTO_CANDIDATES; k++) {
        double fc = fc_candidates[k];
        if (run_filtfilt_trial(y, trial, n, fc) != 0) continue;
        double res = residual_std(y, trial, n);
        if (res <= DISCREPANCY_TOLERANCE * sigma_hat) {
            selected = fc;
            selected_res = res;
            break;
        }
    }

    free(trial);

    printf("# Auto cutoff: noise sigma estimate = %.4e\n", sigma_hat);
    if (selected_res >= 0.0) {
        printf("# Auto cutoff: selected fc = %.4f (residual std = %.4e)\n",
               selected, selected_res);
    } else {
        printf("# Auto cutoff: selected fc = %.4f (fallback, no candidate satisfied discrepancy)\n",
               selected);
    }

    return selected;
}

/* Main filtfilt function using biquad cascade with proper IC */
ButterworthResult* butterworth_filtfilt(const double *x, const double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        const GridAnalysis *grid_info)
{
    /* Initialize pointers for safe cleanup */
    ButterworthResult *result = NULL;
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
        double mb = (double)mem_estimate / (1024.0 * 1024.0);
        if (mb >= 1024.0) {
            fprintf(stderr, "Warning: Large dataset (%d points) requires ~%.1f GB RAM\n",
                    n, mb / 1024.0);
        } else {
            fprintf(stderr, "Warning: Large dataset (%d points) requires ~%.0f MB RAM\n",
                    n, mb);
        }
    }

    /* Auto cutoff selection */
    double fc = cutoff_freq;
    if (auto_cutoff > 0) {
        fc = estimate_cutoff_frequency(y, n);
    }

    /* Validate cutoff frequency */
    if (fc <= CUTOFF_FREQ_MIN || fc >= CUTOFF_FREQ_MAX) {
        fprintf(stderr, "ERROR: Cutoff frequency must be in range (%.1f, %.1f), got %.4f\n",
                CUTOFF_FREQ_MIN, CUTOFF_FREQ_MAX, fc);
        return NULL;
    }

    if (fc < FC_MIN_PRACTICAL) {
        fprintf(stderr, "ERROR: Cutoff frequency too small (fc=%.4e < %.4e). "
                "Filter would be numerically ill-conditioned "
                "(poles approach unit circle). "
                "Use a larger fc or a different smoothing method.\n",
                fc, FC_MIN_PRACTICAL);
        return NULL;
    }

    if (fc > CUTOFF_FREQ_INEFFECTIVE_WARN) {
        fprintf(stderr, "Warning: fc = %.4f is close to Nyquist limit. "
                "Filter passes nearly the entire spectrum unattenuated; "
                "consider a smaller fc for meaningful smoothing.\n", fc);
    }

    /* Check grid uniformity */
    if (grid_info->cv > UNIFORMITY_CV_THRESHOLD) {
        fprintf(stderr, "ERROR: Butterworth filter not suitable for "
                "non-uniform grid (CV=%.4f, threshold=%.4f)!\n",
                grid_info->cv, UNIFORMITY_CV_THRESHOLD);
        return NULL;
    }

    if (grid_info->cv > UNIFORMITY_CV_WARNING) {
        printf("# WARNING: Grid is moderately non-uniform (CV=%.4f)\n", grid_info->cv);
        printf("#   Butterworth frequency response may be distorted.\n");
        printf("#   Consider using Tikhonov method (-m 2 -l auto) for better results.\n");
    } else if (grid_info->cv > NEARLY_UNIFORM_CV_THRESHOLD) {
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
    result->y_deriv = NULL;

    result->y_smooth = (double*)malloc((size_t)n * sizeof(double));
    if (result->y_smooth == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (result buffer)\n");
        goto error;
    }

    result->y_deriv = (double*)malloc((size_t)n * sizeof(double));
    if (result->y_deriv == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed (derivative buffer)\n");
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

    /* Verify pole stability (warn if close to unit circle, error if outside) */
    if (check_pole_stability(sections) != 0) {
        goto error;
    }

    /* Compute initial conditions for each biquad */
    double zi_base[NUM_BIQUADS][2];
    for (int i = 0; i < NUM_BIQUADS; i++) {
        compute_biquad_ic(&sections[i], zi_base[i]);
    }

    /* Pad signal — use directly as work buffer */
    size_t padded_len;
    y_work = pad_signal(y, n, pad_len, &padded_len);
    if (y_work == NULL) {
        fprintf(stderr, "ERROR: Signal padding failed\n");
        goto error;
    }

    /* --- FORWARD + BACKWARD FILTERING (filtfilt) --- */
    apply_cascade(sections, zi_base, y_work, padded_len);
    reverse_array_inplace(y_work, padded_len);
    apply_cascade(sections, zi_base, y_work, padded_len);
    reverse_array_inplace(y_work, padded_len);

    /* --- EXTRACT RESULT --- */
    for (int i = 0; i < n; i++) {
        result->y_smooth[i] = y_work[pad_len + i];
    }

    /* Compute first derivatives via 5-point stencils, O(h^4) accuracy.
     * Using h_avg is valid because filter requires CV <= 0.15. */
    compute_derivatives_5pt(result->y_smooth, n, grid_info->h_avg,
                            result->y_deriv);

    /* --- CLEANUP (Success) --- */
    free(y_work);

    return result;

    /* --- ERROR HANDLER --- */
error:
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
        if (result->y_deriv != NULL) {
            free(result->y_deriv);
        }
        free(result);
    }
}
