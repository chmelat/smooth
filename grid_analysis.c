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
#define THRESH_CV_UNIFORM       0.01  /* Tolerance for practical uniformity */
#define THRESH_CV_NOTICEABLE    0.2   /* CV level where non-uniformity is noted */
#define THRESH_CV_HIGH          0.5   /* CV level suggesting adaptive methods */
#define THRESH_CV_SEVERE        1.0   /* CV level indicating unreliable data */

#define UNIFORMITY_DECAY_FACTOR 2.0   /* For score calculation: exp(-CV * factor) */

#define DROPOUT_TOL        0.25  /* |r - round(r)| tolerance, in units of h_base */
#define DROPOUT_FIT_MIN    0.90  /* Min integer_fit to treat the grid as regular */
#define DROPOUT_MIN_SPACES 10    /* Below this the median is not meaningful */
#define REGIME_SHIFT_MIN   1.5   /* Half-to-half median ratio that counts as a shift */
#define JUMP_RATIO         2.0   /* Neighbouring-spacing ratio that counts as a jump */

/* qsort comparator for the median. Returns via two comparisons rather than
 * (int)(da - db), which truncates to 0 for spacings less than 1.0 apart — the
 * common case here. */
static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/* Median of an already-sorted range. */
static double median_sorted(const double *v, int len) {
    return (len % 2) ? v[len/2] : 0.5 * (v[len/2 - 1] + v[len/2]);
}

/* Helper for safe string appending */
static void append_warning(char *buffer, size_t size, const char *msg) {
    size_t len = strlen(buffer);
    if (len < size - 1) {
        snprintf(buffer + len, size - len, "%s", msg);
    }
}

GridAnalysis* analyze_grid(const double *x, int n)
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
        fprintf(stderr, "ERROR: Memory allocation failed in analyze_grid\n");
        return NULL;
    }
    
    analysis->n_points = n;

    /* Handle trivial cases */
    if (n < 2) {
        analysis->is_uniform = 1;
        analysis->reliability_warning = 0;
        analysis->uniformity_score = 1.0;
        return analysis;
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
     * Calculates: Min, Max, STD
     */
    for (i = 0; i < n-1; i++) {
        double h_curr = x[i+1] - x[i];

        /* Check for non-monotonic data */
        if (h_curr <= 0) {
            fprintf(stderr, "ERROR: Non-monotonic x data at index %d (dx=%g)\n", i, h_curr);
            free_grid_analysis(analysis);
            return NULL;
        }

        /* 1. Min/Max updates */
        if (h_curr < analysis->h_min) analysis->h_min = h_curr;
        if (h_curr > analysis->h_max) analysis->h_max = h_curr;

        /* 2. Standard Deviation accumulation */
        double dev = h_curr - analysis->h_avg;
        sum_sq_diff += dev * dev;
    }

    /* Finalize statistics
     *
     * Bessel-corrected sample STD: divide sum_sq_diff by (N-1), where N is
     * the number of samples — here N = n-1 spacings, so the denominator
     * is n-2. For n=2 (single interval) variability is zero by definition.
     */
    if (n > 2) {
        analysis->h_std = sqrt(sum_sq_diff / (n - 2));
    } else {
        analysis->h_std = 0.0;
    }
    analysis->ratio_max_min = analysis->h_max / analysis->h_min;
    analysis->cv = analysis->h_std / analysis->h_avg;
    
    /* Determine uniformity
     * Threshold 0.01 chosen to match practical interpretation: grids with CV <= 0.01
     * are "practically uniform" and behave well with all methods (including Savgol)
     */
    analysis->is_uniform = (analysis->cv <= THRESH_CV_UNIFORM) ? 1 : 0;
    
    if (analysis->is_uniform) {
        analysis->uniformity_score = 1.0;
    } else {
        analysis->uniformity_score = exp(-analysis->cv * UNIFORMITY_DECAY_FACTOR);
    }
    
    /* Dropout detection — see grid_analysis.h for what the fields mean.
     * A lost sample in fixed-rate logging leaves a gap of k*h_base for integer
     * k >= 2, a signature neither cv (overall spread) nor the cluster detector
     * (which needs a >100x jump straddling h_avg) can see.
     *
     * h_base is the MEDIAN, not the mean: the mean is contaminated by the very
     * gaps being looked for. The median holds the true period up to 45% sample
     * loss; beyond that it jumps to 2*h_base and the estimate collapses.
     *
     * Skipped silently on short grids or if the scratch allocation fails; the
     * fields stay zero from calloc. A diagnostic extra must never turn a working
     * analysis into a failure.
     *
     * ponytail: qsort is O(n log n) on top of an O(n) analysis; quickselect
     * would be O(n) if profiling ever shows this matters. */
    if (n - 1 >= DROPOUT_MIN_SPACES) {
        double *hs = (double *)malloc((size_t)(n - 1) * sizeof(double));
        if (hs != NULL) {
            for (i = 0; i < n-1; i++) hs[i] = x[i+1] - x[i];
            qsort(hs, (size_t)(n - 1), sizeof(double), cmp_double);

            int m = n - 1;
            analysis->h_base = median_sorted(hs, m);

            if (analysis->h_base > 0.0) {
                int near = 0;
                /* Scan the ORIGINAL spacings — hs is sorted in place and its
                 * order no longer matches the grid. */
                for (i = 0; i < n-1; i++) {
                    double r = (x[i+1] - x[i]) / analysis->h_base;
                    double k = floor(r + 0.5);
                    if (fabs(r - k) <= DROPOUT_TOL) {
                        near++;
                        if (k >= 2.0) {
                            int missing = (int)k - 1;
                            analysis->n_gaps++;
                            analysis->n_missing += missing;
                            if (missing > analysis->max_run)
                                analysis->max_run = missing;
                        }
                    }
                }
                analysis->integer_fit = (double)near / (double)m;

                /* Regime shift: the median spacing of the first half of the
                 * record against that of the second. This measures directly
                 * what "the sampling regime changed" means, and unlike a count
                 * of spacings sitting at the base period it does not confuse
                 * heavy scattered loss with a rate change — measured, random
                 * 40% loss puts only 55% of spacings at the base period while a
                 * 1 Hz -> 10 Hz change puts 50% there, an overlap rather than a
                 * margin. The half-to-half ratio is 1.00 against 10.00.
                 * hs is re-filled because the sort above destroyed grid order. */
                int half = m / 2;
                analysis->regime_shift = 1.0;
                if (half >= 1) {
                    for (i = 0; i < m; i++) hs[i] = x[i+1] - x[i];
                    qsort(hs, (size_t)half, sizeof(double), cmp_double);
                    qsort(hs + half, (size_t)(m - half), sizeof(double), cmp_double);
                    double m_first  = median_sorted(hs, half);
                    double m_second = median_sorted(hs + half, m - half);
                    if (m_first > 0.0 && m_second > 0.0) {
                        analysis->regime_shift = (m_first > m_second)
                                               ? m_first / m_second
                                               : m_second / m_first;
                    }
                }

                /* Largest local jump between neighbouring spacings. Purely
                 * local — no global reference — because the events this
                 * describes (rate change, interruption, burst) are local. */
                analysis->max_jump = 1.0;
                for (i = 0; i < n-2; i++) {
                    double a = x[i+1] - x[i], b = x[i+2] - x[i+1];
                    double rho = (a > b) ? a / b : b / a;
                    if (rho > JUMP_RATIO) analysis->n_jumps++;
                    if (rho > analysis->max_jump) {
                        analysis->max_jump = rho;
                        analysis->max_jump_x = x[i+1];
                    }
                }

                /* BOTH conditions are required. A shifted median alone would
                 * misclassify a smoothly graded mesh, whose halves differ by a
                 * factor of 13000 while no two neighbouring spacings differ by
                 * more than 1.1x — a gradual trend, not a regime change. */
                analysis->multi_regime =
                    (analysis->regime_shift > REGIME_SHIFT_MIN &&
                     analysis->max_jump > JUMP_RATIO) ? 1 : 0;

                /* A changed sampling regime is not lost samples: without this
                 * condition a 1 Hz -> 10 Hz record claimed 2241 missing
                 * samples and 18.2% coverage when nothing was missing at all
                 * (v5.11.51 behaviour). */
                analysis->has_dropouts =
                    (analysis->integer_fit >= DROPOUT_FIT_MIN &&
                     analysis->n_missing > 0 &&
                     !analysis->multi_regime) ? 1 : 0;
            }
            free(hs);
        }
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
    
    /* Reliability warning text is shown at every verbosity level — callers that
     * gate on reliability_warning rely on the message being visible (e.g.
     * smooth.c calls with verbose=0 inside an `if (reliability_warning)`). */
    if (analysis->reliability_warning) {
        /* warning_msg is multi-line. Prefix EVERY line: with a single prefix
         * the continuation lines landed on stdout as bare text and corrupted
         * the data table for any numeric consumer (fixed in v5.11.47).
         * Only the first line carries the "WARNING: " label. */
        const char *line = analysis->warning_msg;
        const char *label = "WARNING: ";
        while (*line) {
            const char *nl = strchr(line, '\n');
            printf("%s%s%.*s\n", prefix, label,
                   nl ? (int)(nl - line) : (int)strlen(line), line);
            label = "";
            if (nl == NULL) break;
            line = nl + 1;
        }
    }

    if (verbose >= 1) {
        /* Detailed report */
        printf("%s  Standard deviation: %.6e\n", prefix, analysis->h_std);
        /* Printed only when detected: a clean grid should not gain a line
         * saying so, and a genuinely non-uniform grid must not be given a
         * meaningless missing-sample count. */
        /* Exactly one characterization per grid, most specific first. Nothing
         * is printed for a clean uniform grid or a smoothly non-uniform one. */
        if (analysis->multi_regime) {
            printf("%s  Sampling: mixed regimes - median spacing changes %.3gx "
                   "across the record, largest jump %.3gx at x = %.6e\n",
                   prefix, analysis->regime_shift,
                   analysis->max_jump, analysis->max_jump_x);
        } else if (analysis->has_dropouts) {
            printf("%s  Missing samples: %d in %d gap(s), longest run %d "
                   "(base period %.6e, largest jump %.3gx at x = %.6e)\n",
                   prefix, analysis->n_missing, analysis->n_gaps,
                   analysis->max_run, analysis->h_base,
                   analysis->max_jump, analysis->max_jump_x);
        } else if (analysis->max_jump > JUMP_RATIO) {
            printf("%s  Largest spacing jump: %.3gx at x = %.6e\n",
                   prefix, analysis->max_jump, analysis->max_jump_x);
        }
        printf("%s  Recommendation: %s\n", prefix, get_grid_recommendation(analysis));
    }
    
}

/* Free allocated memory */
void free_grid_analysis(GridAnalysis *analysis)
{
    free(analysis);
}
