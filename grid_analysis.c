/* Grid uniformity analysis
 * Implementation of grid analysis utilities
 * V2.0/2025-11-23/ Single-pass optimization, safe strings, macros
 * V1.0/2025-05-27/ Extracted from tikhonov.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grid_analysis.h"

/* Configuration Thresholds */
#define THRESH_CV_UNIFORM       1e-10 /* Tolerance for perfect uniformity */
#define THRESH_CV_NOTICEABLE    0.2   /* CV level where non-uniformity is noted */
#define THRESH_CV_HIGH          0.5   /* CV level suggesting adaptive methods */
#define THRESH_CV_SEVERE        1.0   /* CV level indicating unreliable data */

#define CLUSTER_RATIO_SMALL     0.1   /* Factor of avg for "small" gap */
#define CLUSTER_RATIO_LARGE     10.0  /* Factor of avg for "large" gap */
#define UNIFORMITY_DECAY_FACTOR 2.0   /* For score calculation: exp(-CV * factor) */

/* Helper for safe string appending */
static void append_warning(char *buffer, size_t size, const char *msg) {
    size_t len = strlen(buffer);
    if (len < size - 1) {
        snprintf(buffer + len, size - len, "%s", msg);
    }
}

GridAnalysis* analyze_grid(const double *x, int n, int store_spacings)
{
    GridAnalysis *analysis;
    int i;
    double h_prev = 0.0; /* Stores spacing from previous iteration for cluster detection */
    
    /* Input validation */
    if (x == NULL || n < 1) {
        return NULL;
    }
    
    /* Allocate analysis structure */
    analysis = (GridAnalysis *)calloc(1, sizeof(GridAnalysis));
    if (analysis == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in analyze_grid\n");
        return NULL;
    }
    
    analysis->n_points = n;
    analysis->n_intervals = n - 1;
    
    /* Handle trivial cases */
    if (n < 2) {
        analysis->is_uniform = 1;
        analysis->reliability_warning = 0;
        analysis->uniformity_score = 1.0;
        return analysis;
    }
    
    /* Allocate spacings array if requested */
    if (store_spacings) {
        analysis->spacings = (double *)malloc((n-1) * sizeof(double));
        if (analysis->spacings == NULL) {
            fprintf(stderr, "Error: Memory allocation failed for spacings array\n");
            free(analysis);
            return NULL;
        }
    }
    
    /* * OPTIMIZATION: Calculate h_avg first (O(1)). 
     * This allows us to calculate STD and everything else in a single pass.
     */
    analysis->h_avg = (x[n-1] - x[0]) / (n-1);
    
    /* Initialize accumulators */
    analysis->h_min = INFINITY;
    analysis->h_max = 0.0;
    double sum_sq_diff = 0.0;
    
    /* * MAIN LOOP (Single Pass) 
     * Calculates: Min, Max, STD, Spacings, and Clusters
     */
    for (i = 0; i < n-1; i++) {
        double h_curr = x[i+1] - x[i];

        /* Check for non-monotonic data */
        if (h_curr <= 0) {
            fprintf(stderr, "Error: Non-monotonic x data at index %d (dx=%g)\n", i, h_curr);
            free_grid_analysis(analysis);
            return NULL;
        }

        /* 1. Store spacing if requested */
        if (store_spacings) {
            analysis->spacings[i] = h_curr;
        }

        /* 2. Min/Max updates */
        if (h_curr < analysis->h_min) analysis->h_min = h_curr;
        if (h_curr > analysis->h_max) analysis->h_max = h_curr;

        /* 3. Standard Deviation accumulation */
        double dev = h_curr - analysis->h_avg;
        sum_sq_diff += dev * dev;

        /* 4. Cluster Detection */
        /* Logic: A cluster boundary is a small gap followed by a large gap */
        if (i > 0) {
            if (h_prev < CLUSTER_RATIO_SMALL * analysis->h_avg && 
                h_curr > CLUSTER_RATIO_LARGE * analysis->h_avg) {
                analysis->n_clusters++;
            }
        }
        
        h_prev = h_curr; /* Save for next iteration */
    }

    /* Finalize statistics */
    analysis->h_std = sqrt(sum_sq_diff / (n-1));
    analysis->ratio_max_min = analysis->h_max / analysis->h_min;
    analysis->cv = analysis->h_std / analysis->h_avg;
    
    /* Determine uniformity */
    analysis->is_uniform = (analysis->cv < THRESH_CV_UNIFORM) ? 1 : 0;
    
    if (analysis->is_uniform) {
        analysis->uniformity_score = 1.0;
    } else {
        analysis->uniformity_score = exp(-analysis->cv * UNIFORMITY_DECAY_FACTOR);
    }
    
    /* Assess reliability and generate warnings using SAFE string handling */
    analysis->reliability_warning = 0;
    analysis->warning_msg[0] = '\0';
    
    char temp_msg[256]; /* Temp buffer for formatted parts */

    if (analysis->cv > THRESH_CV_SEVERE) {
        analysis->reliability_warning = 1;
        snprintf(temp_msg, sizeof(temp_msg),
                "SEVERE grid non-uniformity detected: CV = %.2f\n"
                "This may lead to unreliable smoothing results!\n"
                "Consider resampling data to a more uniform grid.",
                analysis->cv);
        append_warning(analysis->warning_msg, sizeof(analysis->warning_msg), temp_msg);
    }
    else if (analysis->cv > THRESH_CV_HIGH) {
        analysis->reliability_warning = 1;
        snprintf(temp_msg, sizeof(temp_msg),
                "HIGH grid non-uniformity: CV = %.2f\n"
                "Adaptive methods are strongly recommended.\n"
                "Consider using smaller regularization parameters.",
                analysis->cv);
        append_warning(analysis->warning_msg, sizeof(analysis->warning_msg), temp_msg);
    }
    else if (analysis->cv > THRESH_CV_NOTICEABLE) {
        analysis->reliability_warning = 1;
        snprintf(temp_msg, sizeof(temp_msg),
                "Significant spacing variation detected: CV = %.2f\n"
                "Adaptive methods may improve results.",
                analysis->cv);
        append_warning(analysis->warning_msg, sizeof(analysis->warning_msg), temp_msg);
    }
    
    /* Append cluster warning if needed */
    if (analysis->n_clusters > 0) {
        analysis->reliability_warning = 1;
        snprintf(temp_msg, sizeof(temp_msg), 
                "\nWARNING: %d abrupt spacing changes detected (possible data clustering).\n"
                "Standard methods may over-smooth clustered regions.", 
                analysis->n_clusters);
        append_warning(analysis->warning_msg, sizeof(analysis->warning_msg), temp_msg);
    }
    
    return analysis;
}

/* Get recommended method based on grid analysis */
const char* get_grid_recommendation(GridAnalysis *analysis)
{
    if (analysis == NULL) {
        return "Invalid analysis";
    }
    
    if (analysis->is_uniform) {
        return "Grid is uniform - all methods suitable";
    }
    /* Note: Using scores roughly derived from CV thresholds */
    else if (analysis->uniformity_score > 0.8) {
        return "Grid is nearly uniform - standard methods work well";
    }
    else if (analysis->uniformity_score > 0.5) {
        return "Moderate non-uniformity - consider adaptive methods";
    }
    else if (analysis->uniformity_score > 0.2) {
        return "High non-uniformity - adaptive methods recommended";
    }
    else {
        return "Severe non-uniformity - consider resampling data";
    }
}

/* Print grid analysis report */
void print_grid_analysis(GridAnalysis *analysis, int verbose, const char *prefix)
{
    if (analysis == NULL) return;
    
    if (prefix == NULL) prefix = "";
    
    /* Basic report (verbose = 0) */
    printf("%sGrid uniformity analysis:\n", prefix);
    printf("%s  n = %d points\n", prefix, analysis->n_points);
    printf("%s  h_min = %.6e, h_max = %.6e, h_avg = %.6e\n",
           prefix, analysis->h_min, analysis->h_max, analysis->h_avg);
    printf("%s  CV = %.3f\n", prefix, analysis->cv);
    printf("%s  Grid type: %s\n", prefix, analysis->is_uniform ? "UNIFORM" : "NON-UNIFORM");
    printf("%s  Uniformity score: %.2f\n", prefix, analysis->uniformity_score);
    
    if (verbose >= 1) {
        /* Detailed report */
        printf("%s  Standard deviation: %.6e\n", prefix, analysis->h_std);
        printf("%s  Detected clusters: %d\n", prefix, analysis->n_clusters);
        printf("%s  Recommendation: %s\n", prefix, get_grid_recommendation(analysis));
        
        /* Warnings */
        if (analysis->reliability_warning) {
            printf("%sWARNING: %s\n", prefix, analysis->warning_msg);
        }
    }
    
    if (verbose >= 2 && analysis->spacings != NULL) {
        /* Full report with spacing details */
        printf("%s  Spacing details:\n", prefix);
        for (int i = 0; i < analysis->n_intervals; i++) {
            printf("%s    h[%d] = %.6e (%.1f%% of average)\n", 
                   prefix, i, analysis->spacings[i], 
                   100.0 * analysis->spacings[i] / analysis->h_avg);
        }
    }
}

/* Check if adaptive methods should be used */
int should_use_adaptive(GridAnalysis *analysis)
{
    if (analysis == NULL) {
        return 0;
    }
    
    /* Recommend adaptive methods for poor uniformity */
    return (analysis->uniformity_score < 0.5 || 
            analysis->n_clusters > 0);
}

/* Free allocated memory */
void free_grid_analysis(GridAnalysis *analysis)
{
    if (analysis != NULL) {
        free(analysis->spacings);
        free(analysis);
    }
}
