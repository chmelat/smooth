#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

/**
 * Test helper functions for smooth test suite
 * Provides statistical and signal analysis utilities
 */

/**
 * Calculate variance of an array
 * @param data Array of values
 * @param n Number of elements
 * @return Sample variance
 */
double calculate_variance(const double *data, int n);

/**
 * Calculate root mean square error between two arrays
 * @param actual Actual values
 * @param predicted Predicted values
 * @param n Number of elements
 * @return RMSE value
 */
double calculate_rmse(const double *actual, const double *predicted, int n);

/**
 * Find index of first local maximum in array
 * A local maximum is where data[i] > data[i-1] and data[i] > data[i+1]
 * @param data Array of values
 * @param n Number of elements
 * @return Index of first maximum, or -1 if not found
 */
int find_first_maximum(const double *data, int n);

/**
 * Find index of first local minimum in array
 * @param data Array of values
 * @param n Number of elements
 * @return Index of first minimum, or -1 if not found
 */
int find_first_minimum(const double *data, int n);

/**
 * Calculate mean of an array
 * @param data Array of values
 * @param n Number of elements
 * @return Mean value
 */
double calculate_mean(const double *data, int n);

/**
 * Calculate standard deviation of an array
 * @param data Array of values
 * @param n Number of elements
 * @return Standard deviation
 */
double calculate_std(const double *data, int n);

/**
 * Add uniform noise to data (pseudo-Gaussian approximation)
 * Generates noise in range [-amplitude, +amplitude]
 * @param y Original clean data
 * @param y_noisy Output array for noisy data (must be pre-allocated)
 * @param n Number of elements
 * @param amplitude Noise amplitude (peak-to-peak is 2*amplitude)
 * @param seed Random seed for reproducibility
 */
void add_noise(const double *y, double *y_noisy, int n, double amplitude, unsigned int seed);

#endif /* TEST_HELPERS_H */
