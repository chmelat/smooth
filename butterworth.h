/*  Butterworth filter for data smoothing
 *  Digital low-pass filter with filtfilt (zero-phase forward-backward filtering)
 *  V1.0/2025-11-03/ Initial implementation
 */

#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include "grid_analysis.h"

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
 *   cutoff_freq - Normalized cutoff frequency (0 < fc < 1)
 *                 fc = f_cutoff / f_Nyquist = f_cutoff / (f_sample / 2)
 *                 where f_sample = 1 / h_avg (h_avg = average spacing)
 *                 fc = 1 corresponds to Nyquist frequency (maximum representable)
 *                 Smaller fc = more smoothing
 *                 Typical range: 0.02 - 0.6
 *   auto_cutoff - If > 0, automatically determine cutoff from data
 *                 (overrides cutoff_freq parameter)
 *   grid_info   - Grid analysis results (used for uniformity check and sample rate)
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
 *   fc = f_cutoff / f_Nyquist = f_cutoff / (f_sample / 2), where:
 *   - f_cutoff = desired cutoff frequency in your physical units
 *   - f_sample = sampling rate = 1/h_avg
 *   - f_Nyquist = f_sample / 2 (Nyquist frequency - maximum representable frequency)
 *   - Must satisfy: 0 < fc < 1 (where fc = 1 corresponds to Nyquist frequency)
 *
 *   IMPORTANT: fc is normalized with respect to Nyquist frequency (fs/2), not fs!
 *              This is the standard convention in digital signal processing.
 *              fc = 1 represents Nyquist frequency - the theoretical limit.
 *              Frequencies above Nyquist cannot be represented in sampled data
 *              and cause aliasing.
 *
 *   Example: If data spacing h_avg = 0.1 sec, then:
 *            f_sample = 10 Hz, f_Nyquist = 5 Hz
 *            To filter out frequencies above 1 Hz: fc = 1/5 = 0.2
 *            To filter at Nyquist frequency: fc = 5/5 = 1.0
 *
 * Cutoff Frequency Guidelines:
 *   fc = 0.02 - 0.10: Heavy smoothing (removes high-frequency details)
 *   fc = 0.10 - 0.30: Moderate smoothing (typical for noisy data)
 *   fc = 0.30 - 0.60: Light smoothing (preserves most details)
 *   fc = 0.60 - 0.90: Minimal smoothing (may not remove enough noise)
 *   fc > 0.90:        Very minimal filtering (close to Nyquist limit)
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
 *   - cutoff_freq <= 0 or >= 1.0 (invalid cutoff)
 *   - Memory allocation fails
 *   - Grid is excessively non-uniform (ratio > 20)
 */
ButterworthResult* butterworth_filtfilt(double *x, double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        GridAnalysis *grid_info);

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
 *   Recommended normalized cutoff frequency (0 < fc < 1)
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
