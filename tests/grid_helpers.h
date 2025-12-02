#ifndef GRID_HELPERS_H
#define GRID_HELPERS_H

#include "../grid_analysis.h"

/**
 * Grid helper functions for smooth test suite
 * Provides functions to create various types of grids for testing
 */

/**
 * Create perfectly uniform grid with spacing h
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param h Spacing between points
 */
void create_uniform_grid(double *x, int n, double x_start, double h);

/**
 * Create nearly uniform grid with small random perturbations
 * Generates CV < 0.01 (perfectly uniform)
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param h_avg Average spacing
 * @param seed Random seed for reproducibility
 */
void create_nearly_uniform_grid(double *x, int n, double x_start, double h_avg, unsigned int seed);

/**
 * Create moderately non-uniform grid
 * Generates 0.05 < CV < 0.15 (suitable for Tikhonov)
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param h_avg Average spacing
 * @param seed Random seed for reproducibility
 */
void create_moderately_nonuniform_grid(double *x, int n, double x_start, double h_avg, unsigned int seed);

/**
 * Create highly non-uniform grid
 * Generates CV > 0.15 (unsuitable for Savgol)
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param h_avg Average spacing
 * @param seed Random seed for reproducibility
 */
void create_highly_nonuniform_grid(double *x, int n, double x_start, double h_avg, unsigned int seed);

/**
 * Create grid and analyze it in one call
 * Convenience function that creates grid and fills GridAnalysis structure
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param h Average spacing
 * @param grid_type 0=uniform, 1=nearly uniform, 2=moderately non-uniform, 3=highly non-uniform
 * @param seed Random seed (ignored for uniform grid)
 * @param analysis Output GridAnalysis structure
 */
void create_and_analyze_grid(double *x, int n, double x_start, double h,
                              int grid_type, unsigned int seed,
                              GridAnalysis *analysis);

/**
 * Create grid with specific CV value
 * Uses trial-and-error to generate grid with target CV ± tolerance
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param h_avg Average spacing
 * @param target_cv Target coefficient of variation
 * @param seed Random seed
 * @return Actual CV achieved (may differ slightly from target)
 */
double create_grid_with_cv(double *x, int n, double x_start, double h_avg,
                           double target_cv, unsigned int seed);

/**
 * Create non-uniform grid using sin² pattern
 * Used in tikhonov and polyfit tests
 * Spacing: h[i] = base_spacing * sin²(i) + min_spacing
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param base_spacing Base spacing multiplier (typically 0.1)
 * @param min_spacing Minimum spacing to add (typically 0.05 or 0.0)
 */
void create_sin_squared_grid(double *x, int n, double x_start,
                             double base_spacing, double min_spacing);

/**
 * Create non-uniform grid with alternating spacing
 * Used in savgol tests - alternates between two spacing values
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param spacing_dense Spacing for even indices (dense regions)
 * @param spacing_sparse Spacing for odd indices (sparse regions)
 */
void create_alternating_grid(double *x, int n, double x_start,
                              double spacing_dense, double spacing_sparse);

/**
 * Create non-uniform grid with sinusoidal variation
 * Used in tikhonov tests for moderate non-uniformity
 * Spacing: h[i] = base_spacing + amplitude * sin(frequency * i)
 * @param x Output array (must be pre-allocated with size n)
 * @param n Number of points
 * @param x_start Starting x value
 * @param base_spacing Average spacing
 * @param amplitude Amplitude of variation
 * @param frequency Frequency of sinusoidal variation
 */
void create_sinusoidal_grid(double *x, int n, double x_start,
                             double base_spacing, double amplitude, double frequency);

#endif /* GRID_HELPERS_H */
