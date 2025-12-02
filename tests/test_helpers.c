#include "test_helpers.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

double calculate_mean(const double *data, int n) {
    if (data == NULL || n <= 0) {
        return 0.0;
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / n;
}

double calculate_variance(const double *data, int n) {
    if (data == NULL || n <= 1) {
        return 0.0;
    }

    double mean = calculate_mean(data, n);
    double sum_sq_diff = 0.0;

    for (int i = 0; i < n; i++) {
        double diff = data[i] - mean;
        sum_sq_diff += diff * diff;
    }

    return sum_sq_diff / (n - 1);  // Sample variance
}

double calculate_std(const double *data, int n) {
    return sqrt(calculate_variance(data, n));
}

double calculate_rmse(const double *actual, const double *predicted, int n) {
    if (actual == NULL || predicted == NULL || n <= 0) {
        return 0.0;
    }

    double sum_sq_error = 0.0;
    for (int i = 0; i < n; i++) {
        double error = actual[i] - predicted[i];
        sum_sq_error += error * error;
    }

    return sqrt(sum_sq_error / n);
}

int find_first_maximum(const double *data, int n) {
    if (data == NULL || n < 3) {
        return -1;  // Need at least 3 points for local maximum
    }

    for (int i = 1; i < n - 1; i++) {
        if (data[i] > data[i - 1] && data[i] > data[i + 1]) {
            return i;
        }
    }

    return -1;  // No maximum found
}

int find_first_minimum(const double *data, int n) {
    if (data == NULL || n < 3) {
        return -1;  // Need at least 3 points for local minimum
    }

    for (int i = 1; i < n - 1; i++) {
        if (data[i] < data[i - 1] && data[i] < data[i + 1]) {
            return i;
        }
    }

    return -1;  // No minimum found
}

void add_noise(const double *y, double *y_noisy, int n, double amplitude, unsigned int seed) {
    if (y == NULL || y_noisy == NULL || n <= 0) {
        return;
    }

    srand(seed);
    for (int i = 0; i < n; i++) {
        // Uniform noise in range [-amplitude, +amplitude]
        double noise = amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);
        y_noisy[i] = y[i] + noise;
    }
}
