#include "grid_helpers.h"
#include <stdlib.h>
#include <math.h>

void create_uniform_grid(double *x, int n, double x_start, double h) {
    if (x == NULL || n <= 0) {
        return;
    }

    for (int i = 0; i < n; i++) {
        x[i] = x_start + i * h;
    }
}

void create_nearly_uniform_grid(double *x, int n, double x_start, double h_avg, unsigned int seed) {
    if (x == NULL || n <= 0) {
        return;
    }

    srand(seed);
    x[0] = x_start;

    // Add small random perturbations (±2% of h_avg) to achieve CV < 0.01
    for (int i = 1; i < n; i++) {
        double perturbation = h_avg * (rand() / (double)RAND_MAX - 0.5) * 0.04;
        x[i] = x[i-1] + h_avg + perturbation;
    }
}

void create_moderately_nonuniform_grid(double *x, int n, double x_start, double h_avg, unsigned int seed) {
    if (x == NULL || n <= 0) {
        return;
    }

    srand(seed);
    x[0] = x_start;

    // Add moderate random perturbations (±20% of h_avg) to achieve 0.05 < CV < 0.15
    for (int i = 1; i < n; i++) {
        double perturbation = h_avg * (rand() / (double)RAND_MAX - 0.5) * 0.4;
        x[i] = x[i-1] + h_avg + perturbation;
    }
}

void create_highly_nonuniform_grid(double *x, int n, double x_start, double h_avg, unsigned int seed) {
    if (x == NULL || n <= 0) {
        return;
    }

    srand(seed);
    x[0] = x_start;

    // Add large random perturbations (±50% of h_avg) to achieve CV > 0.15
    for (int i = 1; i < n; i++) {
        double perturbation = h_avg * (rand() / (double)RAND_MAX - 0.5) * 1.0;
        // Ensure spacing is always positive
        double spacing = h_avg + perturbation;
        if (spacing < 0.1 * h_avg) {
            spacing = 0.1 * h_avg;
        }
        x[i] = x[i-1] + spacing;
    }
}

void create_and_analyze_grid(double *x, int n, double x_start, double h,
                              int grid_type, unsigned int seed,
                              GridAnalysis *analysis) {
    if (x == NULL || analysis == NULL || n <= 0) {
        return;
    }

    switch (grid_type) {
        case 0:
            create_uniform_grid(x, n, x_start, h);
            break;
        case 1:
            create_nearly_uniform_grid(x, n, x_start, h, seed);
            break;
        case 2:
            create_moderately_nonuniform_grid(x, n, x_start, h, seed);
            break;
        case 3:
            create_highly_nonuniform_grid(x, n, x_start, h, seed);
            break;
        default:
            create_uniform_grid(x, n, x_start, h);
            break;
    }

    GridAnalysis *temp = analyze_grid(x, n, 0);
    if (temp) {
        *analysis = *temp;  // Copy structure contents
        free(temp);
    }
}

double create_grid_with_cv(double *x, int n, double x_start, double h_avg,
                           double target_cv, unsigned int seed) {
    if (x == NULL || n <= 0) {
        return 0.0;
    }

    srand(seed);
    x[0] = x_start;

    // Empirical relationship: perturbation_factor ≈ target_cv * 2.5
    double perturbation_factor = target_cv * 2.5;

    for (int i = 1; i < n; i++) {
        double perturbation = h_avg * (rand() / (double)RAND_MAX - 0.5) * 2.0 * perturbation_factor;
        double spacing = h_avg + perturbation;

        // Ensure spacing is always positive
        if (spacing < 0.1 * h_avg) {
            spacing = 0.1 * h_avg;
        }

        x[i] = x[i-1] + spacing;
    }

    // Calculate actual CV achieved
    GridAnalysis *analysis = analyze_grid(x, n, 0);
    double cv = 0.0;
    if (analysis) {
        cv = analysis->cv;
        free(analysis);
    }

    return cv;
}

void create_sin_squared_grid(double *x, int n, double x_start,
                             double base_spacing, double min_spacing) {
    if (x == NULL || n <= 0) {
        return;
    }

    x[0] = x_start;
    for (int i = 1; i < n; i++) {
        double sin_i = sin((double)i);
        double spacing = base_spacing * sin_i * sin_i + min_spacing;
        x[i] = x[i-1] + spacing;
    }
}

void create_alternating_grid(double *x, int n, double x_start,
                              double spacing_dense, double spacing_sparse) {
    if (x == NULL || n <= 0) {
        return;
    }

    x[0] = x_start;
    for (int i = 1; i < n; i++) {
        if (i % 2 == 0) {
            x[i] = x[i-1] + spacing_dense;
        } else {
            x[i] = x[i-1] + spacing_sparse;
        }
    }
}

void create_sinusoidal_grid(double *x, int n, double x_start,
                             double base_spacing, double amplitude, double frequency) {
    if (x == NULL || n <= 0) {
        return;
    }

    x[0] = x_start;
    for (int i = 1; i < n; i++) {
        double variation = amplitude * sin(frequency * i);
        double spacing = base_spacing + variation;
        x[i] = x[i-1] + spacing;
    }
}
