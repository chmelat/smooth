/*  Butterworth filter for data smoothing
 *  Digital low-pass filter with filtfilt (zero-phase forward-backward filtering)
 *  V1.0/2025-11-03/ Initial implementation
 */

#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

/* Structure for Butterworth filter results */
typedef struct {
    double *y_smooth;     /* Smoothed values */
    int n;                /* Number of points */
    int order;            /* Filter order used (fixed at 4) */
    double cutoff_freq;   /* Normalized cutoff frequency used (0 < fc < 0.5) */
    double sample_rate;   /* Effective sample rate from data spacing */
} ButterworthResult;

/* Main Butterworth filtfilt function
 *
 * Performs zero-phase filtering using 4th-order Butterworth low-pass filter.
 * Uses forward-backward filtering (filtfilt) to eliminate phase distortion.
 *
 * IMPORTANT: This method works best with uniformly or nearly-uniformly spaced data!
 * For highly non-uniform grids, consider resampling or using Tikhonov method.
 *
 * Parameters:
 *   x           - Array of x-coordinates (must be strictly monotonic increasing)
 *   y           - Array of y-values to be smoothed
 *   n           - Number of data points (must be >= 20 for reliable filtering)
 *   cutoff_freq - Normalized cutoff frequency (0 < fc < 0.5)
 *                 fc = f_cutoff / f_sample
 *                 where f_sample = 1 / h_avg (h_avg = average spacing)
 *                 Smaller fc = more smoothing
 *                 Typical range: 0.01 - 0.3
 *   auto_cutoff - If > 0, automatically determine cutoff from data
 *                 (overrides cutoff_freq parameter)
 *
 * Returns:
 *   Pointer to ButterworthResult structure containing:
 *   - y_smooth: filtered values
 *   - n: number of points
 *   - order: filter order (always 4)
 *   - cutoff_freq: normalized cutoff frequency used
 *   - sample_rate: effective sample rate (1/h_avg)
 *   Returns NULL on error.
 *
 * Normalized Cutoff Frequency (fc):
 *   fc = f_cutoff / f_sample, where:
 *   - f_cutoff = desired cutoff frequency in your physical units
 *   - f_sample = sampling rate = 1/h_avg
 *   - Must satisfy: 0 < fc < 0.5 (Nyquist criterion)
 *
 *   Example: If data spacing h_avg = 0.1 sec, then f_sample = 10 Hz
 *            To filter out frequencies above 1 Hz: fc = 1/10 = 0.1
 *
 * Cutoff Frequency Guidelines:
 *   fc = 0.01 - 0.05: Heavy smoothing (removes high-frequency details)
 *   fc = 0.05 - 0.15: Moderate smoothing (typical for noisy data)
 *   fc = 0.15 - 0.30: Light smoothing (preserves most details)
 *   fc > 0.30:        Minimal smoothing (may not remove enough noise)
 *
 * Algorithm:
 *   1. Design 4th-order Butterworth filter (bilinear transform)
 *   2. Pad signal at both ends (reflect padding)
 *   3. Apply filter forward (left to right)
 *   4. Apply filter backward (right to left)
 *   5. Result has zero phase distortion and 8th-order effective filtering
 *
 * Notes:
 *   - Memory must be freed using free_butterworth_result()
 *   - Effective filter order after filtfilt: 2 × 4 = 8
 *   - Zero phase lag (no signal delay)
 *   - Edge effects minimized by padding
 *   - For uniform grids, very efficient and accurate
 *
 * Error handling:
 *   Returns NULL if:
 *   - n < 20 (too few points for reliable filtering)
 *   - cutoff_freq <= 0 or >= 0.5 (invalid cutoff)
 *   - Memory allocation fails
 *   - Grid is excessively non-uniform (ratio > 20)
 */
ButterworthResult* butterworth_filtfilt(double *x, double *y, int n,
                                        double cutoff_freq, int auto_cutoff);

/* Automatic cutoff frequency selection
 *
 * Estimates optimal cutoff frequency based on signal characteristics.
 *
 * Parameters:
 *   x - Array of x-coordinates
 *   y - Array of y-values
 *   n - Number of data points
 *
 * Returns:
 *   Recommended normalized cutoff frequency (0 < fc < 0.5)
 *
 * Notes:
 *   - Conservative approach: prefers over-smoothing to under-smoothing
 *   - Works best with n > 50
 *   - For small datasets, returns default fc = 0.1
 */
double estimate_cutoff_frequency(double *x, double *y, int n);

/* Free allocated memory for ButterworthResult structure
 *
 * Parameters:
 *   result - Pointer to ButterworthResult structure to be freed
 *
 * Notes:
 *   - Safe to call with NULL pointer
 */
void free_butterworth_result(ButterworthResult *result);

#endif /* BUTTERWORTH_H */
