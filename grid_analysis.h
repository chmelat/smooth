/*  Grid uniformity analysis
 *  Header file for grid analysis utilities
 *  V1.0/2025-05-27/ Extracted from tikhonov.c for general use
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
 * 
 * Parameters:
 *   x              - Array of x-coordinates (must be strictly monotonic increasing)
 *   n              - Number of data points
 *   store_spacings - If 1, stores spacing array in result (must be freed)
 * 
 * Returns:
 *   Pointer to GridAnalysis structure containing analysis results
 *   Returns NULL on error.
 * 
 * Notes:
 *   - Memory must be freed using free_grid_analysis()
 *   - Warning messages are generated for highly non-uniform grids
 *   - Uniformity score: 1.0 = perfect, < 0.5 = poor
 */
GridAnalysis* analyze_grid(double *x, int n, int store_spacings);

/* Check if grid is uniform within tolerance
 * 
 * Parameters:
 *   x        - Array of x-coordinates
 *   n        - Number of points
 *   h_avg    - Output: average spacing (can be NULL)
 *   tolerance - Relative tolerance for uniformity (default: 1e-10)
 * 
 * Returns:
 *   1 if uniform, 0 if non-uniform
 */
int is_uniform_grid(double *x, int n, double *h_avg, double tolerance);

/* Get recommended method based on grid analysis
 * 
 * Parameters:
 *   analysis - Pointer to GridAnalysis structure
 * 
 * Returns:
 *   String with recommendation (static, do not free)
 */
const char* get_grid_recommendation(GridAnalysis *analysis);

/* Print grid analysis report
 * 
 * Parameters:
 *   analysis     - Pointer to GridAnalysis structure
 *   verbose      - Verbosity level (0=summary, 1=detailed, 2=full)
 *   prefix       - String prefix for each line (e.g., "# ")
 */
void print_grid_analysis(GridAnalysis *analysis, int verbose, const char *prefix);

/* Calculate effective number of points for regularization
 * 
 * For highly non-uniform grids, the effective number of points
 * may be different from the actual number.
 * 
 * Parameters:
 *   analysis - Pointer to GridAnalysis structure
 * 
 * Returns:
 *   Effective number of points (may be fractional)
 */
double effective_grid_points(GridAnalysis *analysis);

/* Free allocated memory for GridAnalysis structure
 * 
 * Parameters:
 *   analysis - Pointer to GridAnalysis structure to be freed
 */
void free_grid_analysis(GridAnalysis *analysis);

/* Utility functions */

/* Calculate optimal window size for given grid
 * 
 * Parameters:
 *   analysis    - Pointer to GridAnalysis structure
 *   min_window  - Minimum allowed window size
 *   max_window  - Maximum allowed window size
 * 
 * Returns:
 *   Recommended window size (always odd)
 */
int optimal_window_size(GridAnalysis *analysis, int min_window, int max_window);

/* Check if adaptive methods should be used
 * 
 * Parameters:
 *   analysis - Pointer to GridAnalysis structure
 * 
 * Returns:
 *   1 if adaptive methods recommended, 0 otherwise
 */
int should_use_adaptive(GridAnalysis *analysis);

#endif /* GRID_ANALYSIS_H */
