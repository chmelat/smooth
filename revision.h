/*
 * revision.h - define the version number and date
 *
 * Version History
 * ---------------
 * v5.11.43 (current): Code audit v5.11.41 fix A2 — parse_timestamp accepted
 *           calendar-impossible dates. The field validation checked numeric
 *           ranges only (day 1-31, month 1-12, ...), so dates like 2025-02-31 or
 *           2025-04-31 passed and timegm() silently normalized them forward
 *           (Feb 31 -> Mar 3), accepting a timestamp with a shifted date. After
 *           timegm() (which normalizes the struct tm in place) the code now
 *           verifies the normalized tm_year/tm_mon/tm_mday still match the
 *           requested date and returns -1 otherwise — no custom leap-year/month
 *           -length logic needed. Valid leap days (2024-02-29) still parse. New
 *           test parse_timestamp_nonexistent_date. 115 tests pass, zero valgrind
 *           leaks. See doc/code-audit-v5.11.41.md. Closes A2.
 * v5.11.42: Code audit v5.11.41 fix A1 — timestamp-mode x/y
 *           desynchronization. In -T mode parse_input fills timestamp_strings[]
 *           and y[] in parallel by row, then convert_timestamps_to_relative
 *           compacted only x[] (dropping invalid timestamps) while y[] kept its
 *           original row order and was merely truncated. Any invalid timestamp
 *           on a row with an otherwise-valid y therefore shifted y against x and
 *           silently dropped trailing valid points (e.g. timestamp of row 3
 *           paired with y of row 2). convert_timestamps_to_relative now takes an
 *           optional double *y_inout and compacts it in lockstep with x (safe
 *           in place: valid_count <= i). New test convert_compacts_parallel_y;
 *           existing callers pass NULL. 114 tests pass, zero valgrind leaks. See
 *           doc/code-audit-v5.11.41.md. Closes A1.
 * v5.11.41: Butterworth audit #15 — fc-adaptive filtfilt padding.
 *           calculate_pad_length previously returned a fixed 14 samples
 *           regardless of fc, but the filter transient decays as r^k with r the
 *           slowest pole radius, so its length ~1/(1-r) grows without bound as
 *           fc shrinks (~35 samples at fc=0.02, ~6900 at fc=1e-4). For fc<~0.05
 *           the transient leaked past the padding as an undocumented edge
 *           artifact. Padding is now sized to PAD_DECAY_FACTOR/(1-r_max)
 *           (5 time constants, ~0.7% residual), floored at the old 14 and
 *           capped at n-1; when the cap engages a "# WARNING ... too short to
 *           absorb the filter transient" line is printed. The pole-radius math
 *           is factored into biquad_pole_radius()/max_pole_radius() (shared with
 *           check_pole_stability); design + pole check now run before
 *           allocation. run_filtfilt_trial (auto-cutoff trials) gets the same
 *           fc-aware padding. New tests: small_cutoff_preserves_edges (ramp,
 *           fc=0.02, n=600 — boundaries within 0.05; old padding left ~0.6) and
 *           small_cutoff_short_data_clamps (clamp path stays finite). 113 tests
 *           pass, zero valgrind leaks. Closes #15 and #12.
 * v5.11.40: Butterworth audit round 2026-06-10 simple fixes
 *           (#8, #16, #17, #18, #19, #3-comment). (#8) compute_biquad_ic no
 *           longer silently zeroes the IC on a singular system — it warns to
 *           stderr so an edge transient is not emitted unannounced. (#16)
 *           estimate_cutoff_frequency now falls back to the largest (weakest)
 *           fc candidate when no candidate satisfies the discrepancy principle,
 *           instead of the more aggressive AUTO_CUTOFF_FALLBACK=0.2, and warns.
 *           (#17) NOISE_MAD_NORMALIZATION corrected 1.6528553 -> 1.6521808 to
 *           actually equal sqrt(6)*0.6745. (#18) removed the dead public type
 *           ButterworthCoeffs. (#19) corrected calculate_pad_length and
 *           residual_std comments to match the code. (#3) inlined the buf[0]
 *           DC-gain invariant proof into apply_cascade. No behaviour change on
 *           valid input except the #16 fallback path; 111 tests pass, zero
 *           valgrind leaks.
 * v5.11.39: Tikhonov GCV fixes. (1) compute_gcv_score_robust used a
 *           "fast" trace approximation n/(1+sqrt(lambda/h^3)) for n>5000, but
 *           its large-lambda asymptotic exponent is wrong (~s^-1/2 instead of
 *           the correct ~s^-1/4 from the eigenvalue integral), underestimating
 *           tr(H) by a growing factor and biasing GCV toward undersmoothing on
 *           large datasets. The O(n) analytical eigenvalue sum — same order as
 *           the band solve itself — is now used for every n. (2) Lambda is
 *           dimensional (scales with h^3 and y amplitude), so the fixed search
 *           range [1e-8, 1] cannot fit every data scale; find_optimal_lambda_gcv
 *           now warns when the selected lambda is pinned to the range edge and
 *           suggests setting -l manually. README size-tier table updated.
 * v5.11.38: Deep-audit fixes (S3 + T1). (S3) smooth.c guarded the
 *           Tikhonov "Data/Total ratio ... Regularization/Total ratio" line
 *           against a degenerate functional: with lambda=0 the fit is exact, so
 *           data_term=reg_term=total_functional=0 and the ratios printed as
 *           0/0 = nan. The line is now skipped when total_functional==0.
 *           (T1) tikhonov.c find_lambda_lcurve used rss_vals[i]==0.0 as an
 *           "invalid point" sentinel, but rss_vals[i]=log(data_term) is
 *           legitimately 0.0 when data_term==1.0 — a valid L-curve point could
 *           be discarded. Replaced with an explicit valid[] flag (freed in
 *           lcurve_cleanup). L-curve is only used for n>20000. No behaviour
 *           change on typical input; 111 tests pass, zero valgrind leaks.
 * v5.11.37: Polyfit cleanup. (1) The dgelss workspace-query call now
 *           checks its info return and bails to cleanup_error instead of
 *           trusting a possibly-garbage work_query (the per-window solves
 *           already checked info; the query did not). (2) The "Worst Vandermonde
 *           condition number" note is relabelled "Worst effective condition
 *           number (after rcond truncation)" — it reports s_max/s_min over the
 *           RETAINED singular values, bounded by 1/rcond, not the raw Vandermonde
 *           conditioning. The dead s_min<=0 fallback (1e16) is dropped since the
 *           smallest retained s.v. is always positive. (3) polyfit.h poly_degree
 *           range reworded ("0 to DPMAX, which is 12"). (4) Stale test comment
 *           in test_polyfit.c corrected (the NULL there is from window_size>n,
 *           not an over-maximum degree). No behaviour change; 111 tests pass.
 * v5.11.36: Savgol cleanup. (1) Boundary-point handling now treats a
 *           coefficient or allocation failure as a hard error — it frees and
 *           returns NULL instead of silently substituting the raw input y[i]
 *           with a zero derivative. This matches the central-point path and the
 *           "reject loudly" rule (the silent path could only trigger on OOM or
 *           an ill-conditioned boundary solve, but the asymmetry was wrong).
 *           (2) savgol.h doc fixes: example usage now passes grid_info (the old
 *           snippets dropped the argument and would not compile); poly_degree
 *           range corrected to 0..12 with the >6 instability warning noted; the
 *           tikhonov_smooth / find_optimal_lambda_gcv signatures in the
 *           non-uniform-grid recommendations updated to include grid_info. No
 *           behaviour change on valid uniform-grid input; 111 tests pass.
 * v5.11.35: Butterworth cleanup. (1) Reordered the grid-uniformity
 *           CV check ahead of auto-cutoff estimation so a non-uniform grid is
 *           rejected immediately instead of after running several trial
 *           filtfilts and emitting stdout diagnostics. (2) Corrected the
 *           derivative comments: the 5-point stencils are O(h^4) only on a
 *           uniform grid and degrade toward first order (O(CV)) as CV
 *           approaches the 0.15 cap — the old "adds only O(CV*h^2)" claim was
 *           optimistic. No behaviour change on uniform grids; 111 tests pass.
 * v5.11.34: Audit V5.0 A1 — unify Tikhonov discretization.
 *           The penalty now uses ONE integral-measure scheme for all grids:
 *           the Gram matrix (D2)^T W D2 with weights w_k=(h_l+h_r)/2. The
 *           AVERAGE branch (lambda/h^4 unweighted sum), the CV=0.15 method
 *           switch, the DiscretizationMethod enum and select_discretization_
 *           method() are removed; build_band_matrix and compute_functional
 *           collapse to a single path (grid_info param dropped from both). On
 *           a uniform grid the matrix reduces to [1,-4,6,-4,1]*lambda/h^3, so
 *           lambda now scales as ~lambda/h^3 (units Length^3) and reg_term
 *           carries the integration measure; GCV trace eigenvalues rescaled by
 *           h_avg. NOT backward compatible: a given numeric -l value smooths
 *           differently than in <=v5.11.33. Tests: removed two branch-specific
 *           tests (uniform_grid_average_method, discretization_threshold_cv015),
 *           added average_branch_integral_measure. README/CLAUDE.md updated.
 *           111 tests pass, zero leaks.
 * v5.11.33: Audit V5.0 A4+A5+A6 — Tikhonov GCV cleanup.
 *           (A4) Removed the ad-hoc over-fitting penalty
 *           `GCV *= exp(10*(tr(H)/n - 0.7))` from `compute_gcv_score_robust`;
 *           `-l auto` now minimizes textbook GCV. Empirically the penalty
 *           changed the selected lambda by at most one grid step (the GCV
 *           denominator (1 - tr(H)/n)^2 already penalizes tr(H)->n), so the
 *           practical effect is negligible while the code and the diagnostic
 *           label ("pGCV" -> "GCV") become simpler. README "Enhanced GCV"
 *           subsection dropped. (A5) `dpbsv` info>0 message no longer
 *           suggests "Try larger lambda" (I + lambda*(D2)^T D2 is SPD for
 *           lambda>=0, so info>0 means ill-conditioning, not small lambda).
 *           (A6) `compute_gcv_score_robust` reuses `grid_info->ratio_max_min`
 *           instead of recomputing h_min/h_max, and the h4/scale dead
 *           computation moved into the large-n branch where it is used.
 *           No test changes; 112 tests pass.
 * v5.11.32: Audit V5.0 A2 — Tikhonov L-curve corner now uses
 *           the regularization seminorm ||D^2 u||^2 instead of the
 *           lambda-scaled penalty. `find_lambda_lcurve` divides
 *           `regularization_term` (which already includes lambda) back
 *           out before log-transforming the L-curve axes; leaving lambda
 *           in added a monotone log(lambda) ramp that shifted the detected
 *           corner. Affects auto-lambda only for n>20000 where the L-curve
 *           cross-check runs. 112 tests pass.
 * v5.11.31: Audit v5.11.22 B6 — `savgol_coefficients` made
 *           `static` and removed from `savgol.h`. Six internal callsites
 *           inside `savgol.c`, no tests, no external consumers; same
 *           pattern as B5. Substantive docstring (parameter semantics,
 *           special-case behaviour for `deriv_order > poly_degree`)
 *           preserved by moving it above the definition in `savgol.c`.
 *           No behavioural change; B6 closed.
 * v5.11.30: Audit v5.11.22 B5 — `estimate_cutoff_frequency`
 *           made `static` and removed from `butterworth.h`. Function had
 *           a single internal callsite (`butterworth_filtfilt`) and no
 *           tests/external consumers; declaration in the header was dead
 *           API surface. No behavioural change; B5 closed.
 * v5.11.29: Butterworth `-f` default flipped to `auto`.
 *           Without `-f`, the program now invokes
 *           `estimate_cutoff_frequency` (Morozov's discrepancy principle)
 *           instead of using the fixed fc=0.2. The selected fc continues
 *           to appear in the output header (already printed via
 *           `result->cutoff_freq`). Numeric `-f <value>` still overrides
 *           and now explicitly sets `auto_cutoff = 0`. Help text and
 *           README updated. No API change; tests unaffected (they call
 *           `butterworth_filtfilt` directly with explicit args).
 * v5.11.28: Audit v5.11.22 B4 (complete) — parser extracted.
 *           The input parser (~290 lines: line-overflow detection,
 *           normal-mode tokenizer with placeholder semantics,
 *           timestamp-mode tokenizer with logical-column model,
 *           timestamp-to-relative conversion) moves out of `main()` into
 *           a new `parser.c`/`parser.h` module exposing `parse_input()` /
 *           `free_parse_result()`. `main()` shrinks from 599 to 315
 *           lines; `smooth.c` from 754 to 466. The MAX_LINE / MAX_COLS
 *           / BUF constants and the `math.h` / `errno.h` includes move
 *           with the parser. No behavioural change to the binary; B4
 *           closed.
 * v5.11.27: Audit v5.11.22 B4 (partial) — main() goto-cleanup.
 *           The 11 duplicated parser-error cleanup blocks in `smooth.c`
 *           collapse into a single `cleanup:` label at the end of `main()`.
 *           `main()` shrinks from 642 to 599 lines; `smooth.c` from 797 to
 *           754. Two pre-existing exit-time leaks fixed as a side effect:
 *           the `n < sp` early-out and all four "<method> failed!" branches
 *           used to call `exit()` without freeing `x`/`y`/`grid_info`
 *           (previously masked by OS reclaim). Parser extraction into a
 *           dedicated `parser.c` module (the second half of B4) is deferred.
 * v5.11.26: Audit v5.11.22 B3 — Tikhonov h_avg de-duplication.
 *           `select_discretization_method` no longer carries a NULL-grid_info
 *           fallback (the production callsite always passes the result of
 *           `analyze_grid`); the helper now only reads `grid_info->cv` and
 *           drops its `x`, `n` parameters. The three open-coded
 *           `(x[n-1]-x[0])/(n-1)` recomputations in `build_band_matrix`,
 *           `compute_functional`, and `compute_gcv_score_robust` are replaced
 *           with `grid_info->h_avg`. `tikhonov_smooth` and
 *           `find_optimal_lambda_gcv` now reject NULL `grid_info` with
 *           `ERROR: Grid info not available` (matches savgol_smooth contract).
 *           Test 11 rewritten: now asserts NULL grid_info returns NULL.
 * v5.11.25: Audit v5.11.22 B2 — const-correctness on public API.
 *           `tikhonov_smooth`, `find_optimal_lambda_gcv`, `savgol_smooth`,
 *           and `polyfit_smooth` now declare `x`, `y`, and (where applicable)
 *           `grid_info` as `const`, matching butterworth_filtfilt. const
 *           propagated through static helpers (build_band_matrix,
 *           compute_functional, compute_derivatives, compute_gcv_score_robust,
 *           find_lambda_lcurve, select_discretization_method, polyfit
 *           boundary-fallback helpers). Source pointer of the `dcopy_` extern
 *           is now `const double *`. No behavioural change.
 * v5.11.24: Audit v5.11.22 B1 — uniform `ERROR:` prefix for hard
 *           failures (function returns NULL or process exit) across all
 *           non-butterworth modules (decomment, grid_analysis, polyfit,
 *           savgol, tikhonov, smooth). Title-Case `Error:` retired; the
 *           convention documented in CLAUDE.md and butterworth.c is now
 *           applied project-wide. `Warning:` (continuing) untouched.
 * v5.11.23: Audit v5.11.22 A1 + C10. (A1) print_grid_analysis
 *           no longer hides reliability_warning behind verbose>=1; the
 *           warning text is now emitted at every verbosity level so that
 *           callers gating on `if (reliability_warning)` actually see the
 *           message (previously smooth.c printed "# Grid analysis warnings:"
 *           with no warning body for non-uniform grids). (C10) README/code
 *           drift around `-f auto`: README updated to document Morozov's
 *           discrepancy principle (implemented since v5.11.3) instead of
 *           the obsolete "currently returns 0.2" placeholder text.
 * v5.11.22: Audit B3 — document Tikhonov size-dependent algorithm
 *           tiers (n<=5000, 5000<n<=20000, n>20000) in README. Covers trace
 *           estimator switch (eigenvalue sum vs n/(1+sqrt(scale)) approximation)
 *           and refinement step. No algorithm change.
 * v5.11.21: Strict whitespace tokenization in normal-mode parser.
 *           Each whitespace-separated token is one logical column (so an
 *           ISO 8601 timestamp `2026-04-29T11:40:00` counts as one column,
 *           not three). Tokens that strtod cannot fully consume, or that
 *           parse to NaN/Inf, are treated as placeholders; rows where the
 *           selected x or y column lands on a placeholder are skipped, with
 *           a `# Skipped N row(s) ...` summary on stdout. Same NaN/Inf
 *           rejection applied to the y-token in timestamp mode.
 * v5.11.20: Audit B9 — fixed parser limits (MAX_LINE=4096,
 *           MAX_COLS=100) now enforced explicitly. Lines exceeding the line
 *           buffer or rows with more than MAX_COLS columns/tokens cause a
 *           hard error with a clear message instead of silently splitting
 *           lines (corrupting line numbering) or truncating columns
 *           (yielding misleading "only N columns" errors).
 * v5.11.19: Audit C2 — polyfit condition/rank diagnostics now
 *           aggregate across all interior windows instead of reporting only
 *           the first window. Tracks the worst observed condition number and
 *           the count of rank-deficient windows; emits one Note per category
 *           at the end of smoothing. No algorithm change.
 * v5.11.18: Audit B7 — polyfit no longer silently substitutes raw
 *           y[i] (and dy=0) when dgelss returns info != 0 for a window.
 *           A counter tracks per-window SVD failures; if any occurred,
 *           a stderr Warning summarises count, total, and percentage at
 *           the end of smoothing. No algorithm change.
 * v5.11.17: Audit C6 — `make test-valgrind` now fails (exit 1) on
 *           definite/indirect leaks and memory errors. Fixed 4 pre-existing
 *           leaks in butterworth tests (missing `free(grid)` in
 *           higher_cutoff_less_smoothing, invalid_cutoff_frequency,
 *           too_few_points, null_inputs). Baseline now zero leaks.
 * v5.11.16: Audit C4 — bump y/y' output precision from 6 to 8
 *           significant figures (`%12.8lG`), matching the existing x format.
 *           Diagnostic header values (lambda, J, fc, ...) keep `%.6lG`.
 * v5.11.15: Audit B2 — clarify that Tikhonov `-l auto` minimizes
 *           penalized GCV (standard GCV plus exp penalty when tr(H)/n > 0.7),
 *           not textbook GCV. Diagnostic output now labels the score `pGCV`
 *           and points to README "Enhanced GCV". No algorithm change.
 * v5.11.14: Audit B11 — Butterworth diagnostic output convention
 *           documented in module header; pole-stability and fc-near-Nyquist
 *           warnings reclassified from stderr to stdout '#' (preserved with
 *           saved data); large-dataset RAM warning kept on stderr (operational,
 *           not data).
 * v5.11.13: Audit B15 — `-k N:M` now works in `-T` mode: N selects timestamp
 *           logical column (default 1), M selects y column (default 2).
 *           Tokenizer-based parser replaces sscanf with split-on-dot workaround.
 * v5.11.12: Fix audit B10 — uniform `decomment` for stdin and files via new
 *           decomment_stream(FILE*, name); strips '#' comments (incl. inline)
 *           consistently from both inputs.
 * v5.11.7:  Butterworth derivative support (`-d` flag) via 5-point O(h^4)
 *           stencils on filtered output, 106 tests.
 * v5.11.6:  Butterworth cosmetic cleanups — drop unused `x` param from
 *           estimate_cutoff_frequency, clarify "Effective sample rate" label,
 *           adaptive MB/GB memory format.
 * v5.11.5:  Rename CUTOFF_FREQ_STABILITY_WARN to _INEFFECTIVE_WARN, clarify
 *           warning text and document practical fc range.
 * v5.11.4:  Butterworth explicit minimum fc (FC_MIN_PRACTICAL = 1e-4) to reject
 *           numerically ill-conditioned inputs.
 * v5.11.3:  Butterworth auto cutoff via Morozov's discrepancy principle
 *           (noise-aware fc selection).
 * v5.11.2:  Butterworth pole-stability check (warns when poles approach unit
 *           circle).
 * v5.11.1:  Fix DST corruption in timestamps (timegm() instead of mktime()),
 *           103 tests.
 * v5.11.0:  True 2nd-order Tikhonov penalty (D^2)^T W D^2, pentadiagonal matrix.
 * v5.10.1:  Butterworth biquad cascade rewrite, analytical IC.
 * v5.7.1:   Added polyfit unit tests, small bug fixes.
 * v5.6:     First unity tests added.
 * v5.5:     Butterworth filter added, Unix filter support, centralized grid
 *           analysis.
 * v5.4:     Tikhonov hybrid discretization (auto-switch at CV=0.15), GCV
 *           improvements.
 * v5.3:     Savitzky-Golay grid uniformity enforcement.
 * v5.2:     Grid analysis module with `-g` flag.
 * v5.1:     Optional derivative output with `-d` flag.
 * v5.0:     Complete modularization.
 */
#define VERSION "5.11.43"
#define REVDATE "2026-06-16"
