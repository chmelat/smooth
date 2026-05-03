/*
 * revision.h - define the version number and date
 *
 * Version History
 * ---------------
 * v5.11.29 (current): Butterworth `-f` default flipped to `auto`.
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
#define VERSION "5.11.29"
#define REVDATE "2026-05-03"
