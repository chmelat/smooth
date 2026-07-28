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

    /* Sampling regime. All local: no global mean is involved, because the
     * events these describe are local. max_jump is the largest ratio between
     * neighbouring spacings, regime_shift how much the typical spacing differs
     * between the two halves of the record. Advisory only; no method reads them. */
    double max_jump;      /* max(h_i,h_i+1)/min(h_i,h_i+1); 1.0 on a uniform grid */
    double max_jump_x;    /* x where max_jump occurs: the point BEFORE the gap */
    int    n_jumps;       /* Neighbouring pairs whose ratio exceeds JUMP_RATIO */
    double regime_shift;  /* Ratio of the median spacing of the first half of the
                           * record to that of the second half (>= 1.0). Measures
                           * directly what "the sampling regime changed" means.
                           * 1.0 for dropouts of any severity, 10.0 for a
                           * 1 Hz -> 10 Hz change. */
    int    multi_regime;  /* 1 if the record mixes two or more sampling regimes */
    double uniformity_score; /* Score 0-1, where 1 is perfectly uniform */

    /* Dropout detection: a nominally regular grid with missing samples leaves
     * gaps that are integer multiples of the base period. Distinct from
     * multi_regime (the sampling rate changed) and from cv (overall spread) — a grid
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
