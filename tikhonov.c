/* Tikhonov regularization for data smoothing
 * Second derivative penalty via (D²)ᵀ W D² (pentadiagonal Gram matrix)
 * V5.5/2026-09-30/ GCV range 1e14*h_avg^3 (32 points); solve on y minus its LS
 *                  line (offset-proof); failed sweep candidates skipped quietly
 *                  with one note; valid flag instead of a 1e20 sentinel;
 *                  near-duplicate x warning (audit B1)
 * V5.4/2026-07-26/ Removed the L-curve method and the n>20000 branch it was
 *                  reachable from; one log-spaced GCV sweep for all n
 * V5.3/2026-06-10/ GCV: analytical trace for all n (removed n>5000 shortcut with
 *                  wrong asymptotic exponent); warn when optimal lambda is pinned
 *                  to the search-range edge
 * V5.2/2026-06-01/ L-curve: use explicit valid[] flag instead of a 0.0 sentinel
 *                  (log(data_term) can legitimately be 0.0)
 * V5.1/2026-05-30/ A1: single integral-measure discretization (removed AVERAGE/LOCAL branch + CV=0.15 switch)
 * V5.0/2026-02-07/ Corrected penalty from 1st to 2nd order: (D²)ᵀWD² pentadiagonal matrix
 * V4.7/2025-11-28/ Fixed boundary condition asymmetry in Local Spacing Method
 * V4.6/2025-11-28/ Fixed critical bugs: discretization consistency, functional computation
 * V4.5/2025-11-21/ Refactored memory management (goto pattern)
 * V4.4/2025-10-13/ Fixed discretization for non-uniform grids
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tikhonov.h"
#include "grid_analysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* LAPACK function declarations */
extern void dpbsv_(char *uplo, int *n, int *kd, int *nrhs, 
                   double *ab, int *ldab, double *b, int *ldb, int *info);

/* Build band matrix: A = I + lambda * (D2)^T W D2 (pentadiagonal Gram matrix) */
static void build_band_matrix(const double *x, int n, double lambda, double *AB, int ldab, int kd)
{
    int j;

    /* Clear matrix memory explicitly just to be safe, 
     * mostly handled by calloc in caller but useful if reused */
    memset(AB, 0, ldab * n * sizeof(double));

    /* Identity matrix (Data fidelity term) */
    for (j = 0; j < n; j++) {
        /* Diagonal element is at row index 'kd' in LAPACK band storage */
        AB[kd + j*ldab] = 1.0;
    }

    if (lambda <= 0.0 || n < 3) return;

    /* Penalty discretizes the integral lambda * integral (u'')^2 dx via the
     * Gram matrix of the second-derivative operator D2:
     *   Matrix = sum_k w_k * d_k^T d_k   (symmetric, pentadiagonal, kd=2)
     * Row k of D2 (interior point k, 1 <= k <= n-2):
     *   (D2 u)_k = (2/(h_l+h_r)) * [u_{k-1}/h_l - u_k*(1/h_l+1/h_r) + u_{k+1}/h_r]
     * Integration weight: w_k = (h_l + h_r)/2.
     * On a uniform grid (h_l = h_r = h) this reduces to the [1,-4,6,-4,1]*lambda/h^3
     * stencil. Natural BCs (D2 u = 0 at the endpoints) are implicit: D2 has no
     * rows for boundary points, so the grid spacing carries through consistently
     * for uniform and non-uniform grids alike. */
    for (int k = 1; k <= n-2; k++) {
        double h_l = x[k] - x[k-1];
        double h_r = x[k+1] - x[k];
        double h_sum = h_l + h_r;
        double w_k = h_sum / 2.0;

        double a = 2.0 / (h_sum * h_l);        /* coeff of u_{k-1} */
        double b = -2.0 / (h_l * h_r);         /* coeff of u_k */
        double c = 2.0 / (h_sum * h_r);        /* coeff of u_{k+1} */

        double lw = lambda * w_k;

        /* Accumulate d_k^T * w_k * d_k (upper triangle only) */
        /* Diagonal */
        AB[kd + (k-1)*ldab] += lw * a * a;
        AB[kd + k*ldab]     += lw * b * b;
        AB[kd + (k+1)*ldab] += lw * c * c;

        /* 1st superdiagonal */
        AB[(kd-1) + k*ldab]     += lw * a * b;
        AB[(kd-1) + (k+1)*ldab] += lw * b * c;

        /* 2nd superdiagonal */
        AB[(kd-2) + (k+1)*ldab] += lw * a * c;
    }
}

static void compute_derivatives(const double *x, const double *y_smooth, int n, double *y_deriv)
{
    if (n < 2) {
        if (n == 1) y_deriv[0] = 0.0;
        return;
    }
    
    if (n == 2) {
        double slope = (y_smooth[1] - y_smooth[0]) / (x[1] - x[0]);
        y_deriv[0] = slope;
        y_deriv[1] = slope;
        return;
    }
    
    /* First derivative by undetermined coefficients on three points: the weights
     * reproduce u, u' and u'' exactly, leaving an O(h^2)u''' error on any
     * spacing. The index-symmetric difference used before was second order only
     * for h_l == h_r; on a non-uniform grid its leading error is
     * (h_r - h_l)/2 * u''. On a uniform grid the interior weights below reduce
     * to the classical -1/(2h), 0, +1/(2h).
     * x is strictly increasing (validated in tikhonov_smooth), so every spacing
     * is > 0. n >= 3 here, so x[2] and x[n-3] are in range. */

    /* Left boundary: one-sided three-point (second order) */
    {
        double h0 = x[1] - x[0];
        double h1 = x[2] - x[1];
        y_deriv[0] = -(2.0*h0 + h1) / (h0 * (h0 + h1)) * y_smooth[0]
                   +       (h0 + h1) / (h0 * h1)       * y_smooth[1]
                   -              h0 / (h1 * (h0 + h1)) * y_smooth[2];
    }

    /* Interior: three-point, spacing-aware */
    for (int i = 1; i < n-1; i++) {
        double h_l = x[i]   - x[i-1];
        double h_r = x[i+1] - x[i];
        y_deriv[i] = -h_r / (h_l * (h_l + h_r)) * y_smooth[i-1]
                   + (h_r - h_l) / (h_l * h_r)  * y_smooth[i]
                   +  h_l / (h_r * (h_l + h_r)) * y_smooth[i+1];
    }

    /* Right boundary: one-sided three-point (second order) */
    {
        double hA = x[n-1] - x[n-2];
        double hB = x[n-2] - x[n-3];
        y_deriv[n-1] =  (2.0*hA + hB) / (hA * (hA + hB)) * y_smooth[n-1]
                     -       (hA + hB) / (hA * hB)       * y_smooth[n-2]
                     +              hA / (hB * (hA + hB)) * y_smooth[n-3];
    }
}

/* Compute data fidelity, regularization, and total functional values.
 * Uses the same integral-measure discretization as build_band_matrix. */
static void compute_functional(const double *x, const double *y, const double *y_smooth, int n, double lambda,
                              double *data_term, double *reg_term, double *total_functional)
{
    /* Data fidelity term: ||y - u||² */
    *data_term = 0.0;
    for (int i = 0; i < n; i++) {
        double residual = y[i] - y_smooth[i];
        *data_term += residual * residual;
    }
    
    /* Regularization term: λ||D²u||² */
    *reg_term = 0.0;
    
    if (lambda > 0.0 && n >= 3) {
        /* Regularization term lambda * integral (u'')^2 dx, discretized with the
         * same Gram-matrix weighting as build_band_matrix (interior points only;
         * natural BCs implicit). Uniform and non-uniform grids share one path. */
        for (int i = 1; i < n-1; i++) {
            double h_left = x[i] - x[i-1];
            double h_right = x[i+1] - x[i];
            double h_sum = h_left + h_right;

            /* Second derivative with (possibly non-uniform) spacing */
            double d2u = (2.0 / h_sum) * (
                y_smooth[i-1] / h_left -
                y_smooth[i] * (1.0/h_left + 1.0/h_right) +
                y_smooth[i+1] / h_right
            );
            /* Weight by local interval length for integration */
            *reg_term += d2u * d2u * h_sum / 2.0;
        }

        *reg_term *= lambda;
    }
    
    *total_functional = *data_term + *reg_term;
}

/* Solver behind tikhonov_smooth(). quiet = 1 suppresses the dpbsv failure
 * message: inside the GCV sweep a failure at large lambda is an expected,
 * counted outcome, not an error. */
static TikhonovResult* smooth_impl(const double *x, const double *y, int n, double lambda,
                                   const GridAnalysis *grid_info, int quiet)
{
    /* Initialize pointers to NULL for safe cleanup */
    TikhonovResult *result = NULL;
    double *AB = NULL;
    double *b = NULL;

    /* Input Validation */
    if (x == NULL || y == NULL || n < 1 || lambda < 0) {
        fprintf(stderr, "ERROR: Invalid input parameters\n");
        return NULL;
    }

    if (grid_info == NULL) {
        fprintf(stderr, "ERROR: Grid info not available\n");
        return NULL;
    }

    /* Sanity check on x monotonicity */
    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1]) {
            fprintf(stderr, "ERROR: x array must be strictly increasing\n");
            return NULL;
        }
    }

    /* --- ALLOCATION --- */

    /* 1. Result Structure */
    result = (TikhonovResult *)malloc(sizeof(TikhonovResult));
    if (!result) {
        fprintf(stderr, "ERROR: Memory allocation failed (struct)\n");
        goto error;
    }
    result->y_smooth = NULL;
    result->y_deriv = NULL;

    /* 2. Data Arrays */
    result->n = n;
    result->lambda = lambda;
    result->y_smooth = (double *)malloc(n * sizeof(double));
    result->y_deriv = (double *)malloc(n * sizeof(double));

    if (!result->y_smooth || !result->y_deriv) {
        fprintf(stderr, "ERROR: Memory allocation failed (arrays)\n");
        goto error;
    }

    /* 3. Solver Arrays */
    /* Band matrix storage: 
     * LDA = KD + 1 = 2 (1 superdiagonal + 1 diagonal) 
     * Dimensions: ldab * n 
     */
    int kd = 2;
    int ldab = kd + 1;
    
    AB = (double *)calloc(ldab * n, sizeof(double));
    b = (double *)malloc(n * sizeof(double));

    if (!AB || !b) {
        fprintf(stderr, "ERROR: Memory allocation failed (solver buffers)\n");
        goto error;
    }

    /* --- CALCULATION --- */

    /* Prepare RHS: y minus its least-squares line. Lines are in the null
     * space of D2 (the 3-point stencil is exact for quadratics on any grid),
     * so smoothing y - line and adding the line back is the same smoother.
     * It keeps dpbsv accurate: its error grows ~1.6e-15 * (lambda/h^3) * |b|,
     * and |y| includes any offset (1e5 Pa lost 3e-2 at lambda/h^3 = 1e10;
     * detrended <= 1e-7 up to 1e14). x is centred for epoch-sized values. */
    double x_mean = 0.0, y_mean = 0.0, sxx = 0.0, sxy = 0.0, slope = 0.0;
    for (int i = 0; i < n; i++) { x_mean += x[i]; y_mean += y[i]; }
    x_mean /= n;
    y_mean /= n;
    for (int i = 0; i < n; i++) {
        sxx += (x[i] - x_mean) * (x[i] - x_mean);
        sxy += (x[i] - x_mean) * (y[i] - y_mean);
    }
    if (sxx > 0.0) slope = sxy / sxx;  /* n == 1: sxx = 0, subtract the mean only */
    for (int i = 0; i < n; i++)
        b[i] = y[i] - (y_mean + slope * (x[i] - x_mean));

    /* Build System Matrix */
    build_band_matrix(x, n, lambda, AB, ldab, kd);
    
    /* Solve using LAPACK dpbsv */
    char uplo = 'U'; /* Upper triangle of symmetric band matrix */
    int nrhs = 1;
    int info;
    
    dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
    
    if (info != 0) {
        if (!quiet) {
            fprintf(stderr, "ERROR: LAPACK dpbsv failed (info=%d)\n", info);
            if (info > 0) {
                fprintf(stderr, "Leading minor of order %d not positive definite. "
                                "I + lambda*(D2)^T W D2 is SPD in exact arithmetic; this is "
                                "numerical ill-conditioning (lambda too large for this grid).\n", info);
            }
        }
        goto error;
    }
    
    /* Add the line back */
    for (int i = 0; i < n; i++)
        result->y_smooth[i] = b[i] + y_mean + slope * (x[i] - x_mean);
    
    /* Post-processing */
    compute_derivatives(x, result->y_smooth, n, result->y_deriv);
    
    compute_functional(x, y, result->y_smooth, n, lambda,
                      &result->data_term, &result->regularization_term,
                      &result->total_functional);
    
    /* Cleanup temporary buffers (Success path) */
    free(AB);
    free(b);

    return result;

    /* --- ERROR HANDLER --- */
error:
    if (AB) free(AB);
    if (b) free(b);
    free_tikhonov_result(result); /* Handles internal NULLs safely */
    return NULL;
}

TikhonovResult* tikhonov_smooth(const double *x, const double *y, int n, double lambda,
                                const GridAnalysis *grid_info)
{
    /* Near-duplicate x: at a gap g the penalty coefficients grow as 1/g and the
     * solve loses accuracy. Measured against 50-digit arithmetic (gap in units
     * of h): 1e-3 -> <= 5e-5 up to lambda/h^3 = 1e8; 1e-4 -> 3e-6 .. 5e-3;
     * 1e-5 -> 2e-4 .. 4e-2. Hence the 1e-4 threshold. ponytail: a warning, not
     * a fix -- merging such samples would change the data, so that is left to
     * the user. */
    if (grid_info != NULL && grid_info->h_min < 1e-4 * grid_info->h_avg) {
        printf("# WARNING: near-duplicate x: smallest spacing %.3e is %.1e x h_avg; the Tikhonov\n"
               "#          solve loses accuracy near it. Merge or drop near-duplicate samples.\n",
               grid_info->h_min, grid_info->h_min / grid_info->h_avg);
    }
    return smooth_impl(x, y, n, lambda, grid_info, 0);
}

/* GCV score for a single lambda. Returns 0 and sets *score on success, 1 if
 * the solve failed (ill-conditioned at this lambda), 2 if the GCV denominator
 * collapsed (tr(H) ~ n). A status flag, not a sentinel score: a finite marker
 * such as 1e20 is a real GCV value once y is large enough. */
static int gcv_score(const double *x, const double *y, int n, double lambda,
                     const GridAnalysis *grid_info, double *score)
{
    TikhonovResult *result;
    double rss = 0.0;
    double trace_H;

    result = smooth_impl(x, y, n, lambda, grid_info, 1);

    if (result == NULL) {
        fprintf(stderr, "# lambda=%9.3e: solve failed (ill-conditioned), skipped\n", lambda);
        return 1;
    }
    
    for (int i = 0; i < n; i++) {
        double residual = y[i] - result->y_smooth[i];
        rss += residual * residual;
    }
    
    double h_avg = grid_info->h_avg;

    /* Analytical trace from the uniform-grid eigenvalue model of the penalty
     * (approximate for non-uniform grids). Penalty matrix K = sum_k w_k d_k^T d_k
     * uses the integral measure, so on a uniform grid its eigenvalues are
     * h_avg * (4 sin^2(theta/2)/h^2)^2. Null space of D2 is 2-dimensional
     * (constants + linear), so trace starts at 2.0. The O(n) sum is the same
     * order as the band solve, so it is used for every n. */
    trace_H = 2.0;
    for (int k = 1; k <= n-2; k++) {
        double theta = M_PI * k / n;
        double sin_half = sin(theta / 2.0);
        double ev1 = 4.0 * sin_half * sin_half / (h_avg * h_avg);
        double eigenval = ev1 * ev1 * h_avg;
        trace_H += 1.0 / (1.0 + lambda * eigenval);
    }

    /* Standard GCV */
    double denom = 1.0 - trace_H / n;
    int status = 2;
    if (denom > 1e-8) {
        *score = (rss / n) / (denom * denom);
        status = 0;
        fprintf(stderr, "# lambda=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f (%.2f), GCV=%9.3e\n",
                lambda, result->total_functional, rss, trace_H, trace_H / n, *score);
    }

    free_tikhonov_result(result);
    return status;
}

/* Enhanced lambda selection with multiple methods */
double find_optimal_lambda_gcv(const double *x, const double *y, int n, const GridAnalysis *grid_info)
{
    double best_lambda = 0.01;
    double best_gcv = 0.0;
    if (grid_info == NULL) {
        fprintf(stderr, "ERROR: Grid info not available\n");
        return best_lambda;
    }

    /* Search range of the log-spaced GCV sweep, scaled by h^3.
     *
     * Lambda is dimensional, so a fixed range cannot fit every dataset. The
     * penalty eigenvalues used below are 16 sin^4(theta/2) / h^3, so the
     * product lambda * eigenvalue -- the only thing the smoother 1/(1+lambda*
     * eigenvalue) actually sees -- is dimensionless exactly when lambda
     * carries h^3. Scaling the bounds that way makes the search grid-scale
     * invariant instead of merely wide.
     *
     * Lambda does not depend on the y amplitude (the minimizer is linear in
     * y). What sets it is how many samples a signal feature spans: measured,
     * the optimum grows as P^4 for a period of P samples -- lambda/h^3 ~ 8 at
     * P = 20, 6e5 at P = 500, 2e10 at P = 10000 -- and GCV tracks that
     * optimum closely wherever the range lets it. The old upper bound 1e6*h^3
     * pinned every signal slower than ~500 samples per period (audit B1).
     *
     * ponytail: the upper bound 1e14*h^3 is set by dpbsv, not by the data.
     * On a uniform grid the factorization fails near 1e16*h^3 and, with the
     * detrended RHS (tikhonov_smooth), the error stays <= 3e-8 up to 1e14
     * (measured against 40-digit arithmetic). On non-uniform grids the
     * conditioning follows the smallest LOCAL spacing, so the top candidates
     * may fail there; they are skipped and counted instead of shrinking the
     * range for the whole record -- an h_min-scaled bound let one
     * near-duplicate sample cut it below the old 1e6*h^3. Signals slower
     * than ~1e5 samples per period would need another solver. The lower
     * bound keeps its margin for grids with h_min << h_avg. */
    const double h3 = grid_info->h_avg * grid_info->h_avg * grid_info->h_avg;
    const double lambda_min = 1e-8 * h3;
    const double lambda_max = 1e14 * h3;

    if (n < 3) {
        fprintf(stderr, "Warning: Too few points for GCV (n=%d)\n", n);
        return best_lambda;
    }

    /* The per-lambda search trace goes to stderr: it is progress output, not
     * something to preserve alongside the smoothed data. The chosen lambda and
     * the caveats about its reliability stay on stdout (see the diagnostic
     * output convention in the smooth-dev-tasks skill). The '#' prefix is kept
     * so a caller merging the streams with 2>&1 still sees valid comments. */
    fprintf(stderr, "# GCV optimization for n=%d points (see README: Generalized Cross Validation)\n", n);
    fprintf(stderr, "# Grid CV = %.3f (integral-measure penalty discretization)\n",
            grid_info->cv);

    /* Both caveats qualify the reliability of the lambda that is saved with the
     * data, so both are stdout — and both print once. The ratio note used to
     * live in the per-lambda GCV score, where the sweep repeated it verbatim
     * once per trial lambda (21 identical lines on a 60-point mesh). */
    if (grid_info->cv > 0.2) {
        printf("# WARNING: Highly non-uniform grid detected. Trace approximation less accurate.\n");
    }
    if (grid_info->ratio_max_min > 2.0) {
        printf("# Note: Trace(H) approximation less accurate for non-uniform grid (ratio=%.2f)\n",
               grid_info->ratio_max_min);
    }
    
    /* Log-spaced GCV search over [lambda_min, lambda_max].
     *
     * This sweep is the whole search. A sub-grid refinement pass used to run
     * on top of it, gated on n <= 5000 -- a threshold left over from the trace
     * shortcut removed in v5.11.39 (TK1), which made identical data land on a
     * lambda differing by ~28% either side of n = 5000. Measured before
     * removing it: refinement cost 16% (n=8k) to 36% (n=100k) of total runtime
     * and moved the smoothed output by at most 0.21-0.33% of the data range
     * (RMS 0.03-0.05%), because the GCV curve is flat near its minimum -- the
     * objective itself improved by under 0.1%. That is well below the noise
     * being smoothed, so the sweep grid is resolution enough and every n now
     * follows the same path. Audit TK8. */
    /* ~0.7 decades per step over the 22 decades. Measured on 8 datasets with
     * known ground truth: finer than 0.7 gained 1-2% RMSE for 22% more
     * runtime, coarser (1.1) lost 5-11%. */
    int n_points = 32;
    int have_best = 0, n_failed = 0;
    double top_ok = 0.0;  /* largest candidate the solver handled */

    for (int i = 0; i < n_points; i++) {
        double log_lambda = log10(lambda_min) + (log10(lambda_max) - log10(lambda_min)) * i / (n_points - 1);
        double lambda_test = pow(10.0, log_lambda);
        double gcv;

        int status = gcv_score(x, y, n, lambda_test, grid_info, &gcv);
        if (status == 1) n_failed++;
        if (status != 1) top_ok = lambda_test;
        if (status != 0) continue;

        if (!have_best || gcv < best_gcv) {
            best_gcv = gcv;
            best_lambda = lambda_test;
            have_best = 1;
        }
    }

    if (n_failed > 0) {
        printf("# Note: %d of %d GCV candidates skipped: the solve is ill-conditioned at large lambda on this grid.\n",
               n_failed, n_points);
    }

    /* Every candidate failed (dpbsv or a collapsed GCV denominator): there is
     * no GCV answer. Return the least-smoothing lambda, which cannot distort
     * the data, and say so. */
    if (!have_best) {
        printf("# WARNING: GCV failed for every candidate lambda; returning lambda = %.3e (no smoothing).\n",
               lambda_min);
        printf("#          Set lambda manually (-l <value>).\n");
        return lambda_min;
    }

    /* Upper edge: the optimum may lie beyond it, so say so. Lower edge: GCV
     * has a flat limit as lambda -> 0 (RSS and (1 - tr/n)^2 both go as
     * lambda^2), so on a near-uniform grid a minimum there means "do not
     * smooth" -- noise-free data, output already equal to the input -- and is
     * a note. On a non-uniform grid (the trace caveats above were printed)
     * the h_avg trace model may be what drove GCV there, so it stays a
     * warning. The upper edge is the largest candidate the solver handled. */
    if (best_lambda >= top_ok * 0.99) {
        printf("# WARNING: optimal lambda = %.3e lies at the upper edge of the %s [%.0e, %.0e].\n",
               best_lambda, n_failed ? "solvable search range" : "search range", lambda_min, top_ok);
        printf("#          The true optimum may lie beyond it; consider setting lambda manually (-l <value>).\n");
    } else if (best_lambda <= lambda_min * 1.01) {
        if (grid_info->cv > 0.2 || grid_info->ratio_max_min > 2.0)
            printf("# WARNING: optimal lambda = %.3e lies at the lower edge of the search range;\n"
                   "#          on this grid that may reflect the trace approximation (see above). Consider -l <value>.\n",
                   best_lambda);
        else
            printf("# Note: GCV prefers no smoothing (lambda at the lower search edge); output ~ input.\n");
    }

    fprintf(stderr, "# Optimal lambda: %.6e (GCV=%.3e)\n", best_lambda, best_gcv);
    return best_lambda;
}

void free_tikhonov_result(TikhonovResult *result)
{
    if (result != NULL) {
        if (result->y_smooth) free(result->y_smooth);
        if (result->y_deriv) free(result->y_deriv);
        free(result);
    }
}
