/* Butterworth filter for data smoothing
 * Digital low-pass filter with filtfilt (zero-phase forward-backward filtering)
 * V1.4/2025-12-07/ Biquad cascade with proper initial conditions (no LAPACK)
 * V1.3/2025-12-06/ Memory optimization, symbolic constants, improved validation
 */

#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include <stddef.h>
#include "grid_analysis.h"

/* Filter parameters */
#define BUTTERWORTH_ORDER      4
#define NUM_BIQUADS            2   /* 4th order = 2 cascaded 2nd-order sections */
#define BUTTERWORTH_NUM_COEFFS 3   /* Each biquad has 3 coefficients (b0,b1,b2 / a0,a1,a2) */

/* Memory limits */
#define BUTTERWORTH_MAX_POINTS_WARNING  (50 * 1000 * 1000)
#define BUTTERWORTH_MAX_POINTS_LIMIT    (500 * 1000 * 1000)
#define BUTTERWORTH_MIN_POINTS          20

/* Structure for one biquad section (2nd-order IIR filter) */
typedef struct {
    double b[3];  /* Numerator coefficients: [b0, b1, b2] */
    double a[3];  /* Denominator coefficients: [a0=1, a1, a2] */
} BiquadSection;

/* Structure for complete filter (2 cascaded biquads for 4th order) */
typedef struct {
    BiquadSection sections[NUM_BIQUADS];
} ButterworthCoeffs;

/* Structure for Butterworth filter results */
typedef struct {
    double *y_smooth;     /* Smoothed values */
    int n;                /* Number of points */
    int order;            /* Filter order (BUTTERWORTH_ORDER) */
    double cutoff_freq;   /* Normalized cutoff frequency (0 < fc < 1) */
    double sample_rate;   /* Effective sample rate from data spacing */
} ButterworthResult;

/* Main Butterworth filtfilt function
 *
 * Zero-phase filtering using 4th-order Butterworth low-pass filter
 * implemented as cascade of 2 biquad sections for numerical stability.
 *
 * Parameters:
 *   x           - Array of x-coordinates (strictly monotonic increasing)
 *   y           - Array of y-values to be smoothed
 *   n           - Number of data points (>= BUTTERWORTH_MIN_POINTS)
 *   cutoff_freq - Normalized cutoff frequency (0 < fc < 1)
 *                 fc = f_cutoff / f_Nyquist
 *   auto_cutoff - If > 0, automatically determine cutoff (overrides cutoff_freq)
 *   grid_info   - Grid analysis results (uniformity check, sample rate)
 *
 * Returns:
 *   Pointer to ButterworthResult, or NULL on error.
 *   Free with free_butterworth_result().
 *
 * Implementation notes:
 *   - Uses biquad cascade (2 × 2nd-order) instead of monolithic 4th-order
 *   - Proper initial conditions computed analytically (no LAPACK dependency)
 *   - Odd-reflection padding to minimize edge effects
 *   - Effective filter order after filtfilt: 2 × 4 = 8
 */
ButterworthResult* butterworth_filtfilt(const double *x, const double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        const GridAnalysis *grid_info);

/* Automatic cutoff frequency estimation */
double estimate_cutoff_frequency(const double *x, const double *y, int n);

/* Free allocated memory */
void free_butterworth_result(ButterworthResult *result);

#endif /* BUTTERWORTH_H */
