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
    int is_uniform;         /* 1 if CV <= 0.01 (practical uniformity), 0 otherwise */
    int reliability_warning; /* 1 if reliability concerns exist */
    char warning_msg[512];  /* Warning message buffer */
    
    /* Additional statistics */
    int n_points;           /* Number of points */
    int n_clusters;         /* Number of detected clusters */
    double uniformity_score; /* Score 0-1, where 1 is perfectly uniform */

    /* Dropout detection: a nominally regular grid with missing samples leaves
     * gaps that are integer multiples of the base period. Distinct from
     * n_clusters (abrupt spacing change) and from cv (overall spread) — a grid
     * can have any combination of the three. Advisory only: no smoothing method
     * reads these fields. All zero when has_dropouts == 0.
     *
     * Known limit: a mesh whose spacings are genuinely an exact multiple of one
     * another (e.g. mostly h with an occasional 4h by design, unbalanced enough
     * that the median lands on h) is indistinguishable from dropouts by
     * construction and will be reported as missing samples. */
    double h_base;       /* Median spacing — robust base period estimate */
    double integer_fit;  /* Fraction of spacings within DROPOUT_TOL of k*h_base */
    int    n_gaps;       /* Gaps that are an integer multiple k >= 2 of h_base */
    int    n_missing;    /* Estimated missing samples, sum of (k-1) */
    int    max_run;      /* Longest run of consecutive missing samples */
    int    has_dropouts; /* 1 if integer_fit >= DROPOUT_FIT_MIN and n_missing > 0 */
} GridAnalysis;

/* Analyze grid uniformity and return detailed statistics
 * * Parameters:
 * x              - Array of x-coordinates (must be strictly monotonic increasing)
 * n              - Number of data points
 * * Returns:
 * Pointer to GridAnalysis structure containing analysis results
 * Returns NULL on error.
 */
GridAnalysis* analyze_grid(const double *x, int n);

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

#endif /* GRID_ANALYSIS_H */
