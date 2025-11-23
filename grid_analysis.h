/* Grid uniformity analysis
 * Header file for grid analysis utilities
 * V2.0/2025-11-23/ Cleaned up, optimized, safe string handling
 * V1.0/2025-05-27/ Extracted from tikhonov.c for general use
 */

#ifndef GRID_ANALYSIS_H
#define GRID_ANALYSIS_H

/* Structure for grid uniformity analysis */
typedef struct {
    double h_min;           /* Minimum spacing */
    double h_max;           /* Maximum spacing */
    double h_avg;           /* Average spacing */
    double h_std;           /* Standard deviation of spacing */
    double ratio_max_min;   /* Ratio h_max/h_min */
    double cv;              /* Coefficient of variation (std/avg) */
    int is_uniform;         /* 1 if uniform, 0 otherwise */
    int reliability_warning; /* 1 if reliability concerns exist */
    char warning_msg[512];  /* Warning message buffer */
    
    /* Additional statistics */
    int n_points;           /* Number of points */
    int n_intervals;        /* Number of intervals (n_points - 1) */
    double *spacings;       /* Array of spacings (can be NULL if not requested) */
    int n_clusters;         /* Number of detected clusters */
    double uniformity_score; /* Score 0-1, where 1 is perfectly uniform */
} GridAnalysis;

/* Analyze grid uniformity and return detailed statistics
 * * Parameters:
 * x              - Array of x-coordinates (must be strictly monotonic increasing)
 * n              - Number of data points
 * store_spacings - If 1, stores spacing array in result (must be freed)
 * * Returns:
 * Pointer to GridAnalysis structure containing analysis results
 * Returns NULL on error.
 */
GridAnalysis* analyze_grid(const double *x, int n, int store_spacings);

/* Get recommended method based on grid analysis
 * * Returns:
 * String with recommendation (static, do not free)
 */
const char* get_grid_recommendation(GridAnalysis *analysis);

/* Print grid analysis report
 * * Parameters:
 * analysis     - Pointer to GridAnalysis structure
 * verbose      - Verbosity level (0=summary, 1=detailed, 2=full)
 * prefix       - String prefix for each line (e.g., "# ")
 */
void print_grid_analysis(GridAnalysis *analysis, int verbose, const char *prefix);

/* Free allocated memory for GridAnalysis structure */
void free_grid_analysis(GridAnalysis *analysis);

/* Check if adaptive methods should be used */
int should_use_adaptive(GridAnalysis *analysis);

#endif /* GRID_ANALYSIS_H */
