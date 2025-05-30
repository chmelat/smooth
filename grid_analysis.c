/*  Grid uniformity analysis
 *  Implementation of grid analysis utilities
 *  V1.0/2025-05-27/ Extracted from tikhonov.c for general use
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grid_analysis.h"

/* Analyze grid uniformity and return detailed statistics */
GridAnalysis* analyze_grid(double *x, int n, int store_spacings)
{
    GridAnalysis *analysis;
    int i;
    
    /* Input validation */
    if (x == NULL || n < 1) {
        return NULL;
    }
    
    /* Allocate analysis structure */
    analysis = (GridAnalysis *)calloc(1, sizeof(GridAnalysis));
    if (analysis == NULL) {
        fprintf(stderr, "Memory allocation failed in analyze_grid\n");
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
            fprintf(stderr, "Memory allocation failed for spacings array\n");
            free(analysis);
            return NULL;
        }
    }
    
    /* Calculate spacings and basic statistics */
    analysis->h_min = 1e20;
    analysis->h_max = 0.0;
    analysis->h_avg = 0.0;
    
    for (i = 0; i < n-1; i++) {
        double h = x[i+1] - x[i];
        
        /* Check for non-monotonic data */
        if (h <= 0) {
            fprintf(stderr, "Error: Non-monotonic x data at index %d\n", i);
            free_grid_analysis(analysis);
            return NULL;
        }
        
        if (store_spacings) {
            analysis->spacings[i] = h;
        }
        
        if (h < analysis->h_min) analysis->h_min = h;
        if (h > analysis->h_max) analysis->h_max = h;
        analysis->h_avg += h;
    }
    analysis->h_avg /= (n-1);
    
    /* Calculate standard deviation */
    analysis->h_std = 0.0;
    for (i = 0; i < n-1; i++) {
        double h = x[i+1] - x[i];
        double dev = h - analysis->h_avg;
        analysis->h_std += dev * dev;
    }
    analysis->h_std = sqrt(analysis->h_std / (n-1));
    
    /* Calculate uniformity metrics */
    analysis->ratio_max_min = analysis->h_max / analysis->h_min;
    analysis->cv = analysis->h_std / analysis->h_avg;
    
    /* Determine if grid is uniform */
    analysis->is_uniform = (analysis->cv < 1e-10) ? 1 : 0;
    
    /* Calculate uniformity score (0-1 scale) */
    if (analysis->is_uniform) {
        analysis->uniformity_score = 1.0;
    } else {
        /* Score based on coefficient of variation and ratio */
        double cv_score = exp(-analysis->cv * 2.0);  /* Exponential decay */
        double ratio_score = 1.0 / (1.0 + log(analysis->ratio_max_min));
        analysis->uniformity_score = 0.7 * cv_score + 0.3 * ratio_score;
        
        /* Ensure score is in [0,1] */
        if (analysis->uniformity_score < 0.0) analysis->uniformity_score = 0.0;
        if (analysis->uniformity_score > 1.0) analysis->uniformity_score = 1.0;
    }
    
    /* Detect clusters */
    analysis->n_clusters = 0;
    for (i = 0; i < n-2; i++) {
        double h1 = x[i+1] - x[i];
        double h2 = x[i+2] - x[i+1];
        
        /* Detect abrupt changes in spacing */
        if (h1 < 0.1 * analysis->h_avg && h2 > 10.0 * analysis->h_avg) {
            analysis->n_clusters++;
        }
    }
    
    /* Assess reliability and generate warnings */
    analysis->reliability_warning = 0;
    analysis->warning_msg[0] = '\0';
    
    /* Severe non-uniformity warnings */
    if (analysis->ratio_max_min > 100.0) {
        analysis->reliability_warning = 1;
        sprintf(analysis->warning_msg, 
                "SEVERE grid non-uniformity detected: h_max/h_min = %.1f\n"
                "This may lead to unreliable smoothing results!\n"
                "Consider resampling data to a more uniform grid.", 
                analysis->ratio_max_min);
    }
    else if (analysis->ratio_max_min > 20.0) {
        analysis->reliability_warning = 1;
        sprintf(analysis->warning_msg, 
                "HIGH grid non-uniformity: h_max/h_min = %.1f\n"
                "Adaptive methods are strongly recommended.\n"
                "Consider using smaller regularization parameters.", 
                analysis->ratio_max_min);
    }
    else if (analysis->cv > 0.5) {
        analysis->reliability_warning = 1;
        sprintf(analysis->warning_msg, 
                "Significant spacing variation detected: CV = %.2f\n"
                "Adaptive methods may improve results.", 
                analysis->cv);
    }
    
    /* Check for clustering */
    if (analysis->n_clusters > 0) {
        analysis->reliability_warning = 1;
        char cluster_msg[256];
        sprintf(cluster_msg, 
                "\nWARNING: %d abrupt spacing changes detected (possible data clustering).\n"
                "Standard methods may over-smooth clustered regions.", 
                analysis->n_clusters);
        
        /* Append to existing warning if space permits */
        if (strlen(analysis->warning_msg) + strlen(cluster_msg) < 510) {
            strcat(analysis->warning_msg, cluster_msg);
        }
    }
    
    return analysis;
}

/* Check if grid is uniform within tolerance */
int is_uniform_grid(double *x, int n, double *h_avg, double tolerance)
{
    int i;
    double h_local, h_mean;
    
    if (n < 2) return 1;
    
    /* Calculate average spacing */
    h_mean = (x[n-1] - x[0]) / (n - 1);
    
    if (h_avg != NULL) {
        *h_avg = h_mean;
    }
    
    /* Check uniformity */
    for (i = 1; i < n; i++) {
        h_local = x[i] - x[i-1];
        if (fabs(h_local - h_mean) > tolerance * fabs(h_mean)) {
            return 0;
        }
    }
    
    return 1;
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
    printf("%s  h_max/h_min = %.2f, CV = %.3f\n", 
           prefix, analysis->ratio_max_min, analysis->cv);
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

/* Calculate effective number of points for regularization */
double effective_grid_points(GridAnalysis *analysis)
{
    if (analysis == NULL || analysis->n_points < 2) {
        return 1.0;
    }
    
    /* For uniform grids, effective points = actual points */
    if (analysis->is_uniform) {
        return (double)analysis->n_points;
    }
    
    /* For non-uniform grids, adjust based on uniformity score */
    /* This is a heuristic that reduces effective points for highly non-uniform grids */
    double base_points = (double)analysis->n_points;
    double adjustment = 0.7 + 0.3 * analysis->uniformity_score;
    
    return base_points * adjustment;
}

/* Calculate optimal window size for given grid */
int optimal_window_size(GridAnalysis *analysis, int min_window, int max_window)
{
    int window;
    
    if (analysis == NULL) {
        return min_window;
    }
    
    /* Ensure odd window size */
    if (min_window % 2 == 0) min_window++;
    if (max_window % 2 == 0) max_window--;
    
    if (analysis->is_uniform) {
        /* For uniform grids, use larger windows */
        window = (min_window + max_window) / 2;
    } else {
        /* For non-uniform grids, adjust based on uniformity */
        double factor = 0.3 + 0.7 * analysis->uniformity_score;
        window = min_window + (int)(factor * (max_window - min_window));
    }
    
    /* Ensure odd */
    if (window % 2 == 0) window++;
    
    /* Ensure within bounds */
    if (window < min_window) window = min_window;
    if (window > max_window) window = max_window;
    
    return window;
}

/* Check if adaptive methods should be used */
int should_use_adaptive(GridAnalysis *analysis)
{
    if (analysis == NULL) {
        return 0;
    }
    
    /* Recommend adaptive methods for poor uniformity */
    return (analysis->uniformity_score < 0.5 || 
            analysis->ratio_max_min > 10.0 ||
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
