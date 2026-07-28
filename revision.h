/*
 * revision.h - define the version number and date
 *
 * Version History
 * ---------------
 * v5.11.53 (current): Keep the GCV sweep trace out of the data stream, and
 *           make all program output ASCII (audit findings B8+B9; plan 007).
 *           Three separable defects sharing one set of lines.
 *           (1) BEHAVIOUR CHANGE: the per-lambda GCV search trace moved from
 *           stdout to stderr — the sweep rows, its header, the grid-CV line,
 *           the refinement header and the "Optimal lambda" conclusion. They
 *           keep their '#' prefix, so merging with 2>&1 still yields valid
 *           comments. `smooth -m 2 -l auto data > out.dat` no longer captures
 *           the table; `2>> out.dat` recovers the old behaviour. On stdout
 *           stay only the lines that qualify the lambda saved with the data:
 *           the edge-of-range and highly-non-uniform warnings. Measured on
 *           test_data.dat: 28 stdout comment lines -> 7. This also removes one
 *           of the three places the chosen lambda was printed to stdout
 *           (tikhonov.c's "Optimal lambda"; smooth.c's two lines are
 *           unchanged). No CLI flag was added — compute_gcv_score_robust()
 *           had a `verbose` parameter and lost it in v5.11.43 because both
 *           call sites passed 1; stream selection gets the same result with
 *           no new surface.
 *           (2) The "Trace(H) approximation less accurate" note lived inside
 *           compute_gcv_score_robust(), so the sweep repeated it verbatim once
 *           per trial lambda — 21 identical lines on a 60-point alternating
 *           mesh, where the program emitted 62 diagnostic lines for 60 data
 *           rows. Hoisted into find_optimal_lambda_gcv() next to the cv>0.2
 *           warning: both non-uniformity caveats now print once, on stdout.
 *           compute_gcv_score_robust() is now free of output entirely, the
 *           state v5.11.43 was reaching for. Its `ratio` local went with it.
 *           (3) Seven output strings carried non-ASCII glyphs: five `λ` in
 *           tikhonov.c, `×` in smooth.c's Butterworth cutoff line, and an em
 *           dash in butterworth.c's auto-cutoff warning. All ASCII now
 *           (`lambda`, `*`, `;`); source comments are untouched. Program
 *           output contains zero non-ASCII bytes on either stream.
 *           tikhonov.h:77 claimed the detailed information goes to stdout;
 *           corrected. The convention block in butterworth.c and the
 *           smooth-dev-tasks skill gain the stderr "# ..." category.
 *           No numerics change: data rows, the selected lambda, the sweep
 *           range and the refinement rule are all identical. 138 tests pass.
 * v5.11.52: Replace the cluster detector with sampling-regime
 *           detection, and correct v5.11.51. Two problems with one root
 *           cause — a local phenomenon judged against a global statistic.
 *           First, v5.11.51's dropout gate could not tell lost samples from
 *           a changed sampling regime. On a record where the logger switched
 *           from 1 Hz to 10 Hz it reported "2241 missing samples ...
 *           coverage 18.2%" when nothing was missing at all: both halves
 *           have spacings that are integer multiples of the median, so
 *           integer_fit reached 1.00 and the gate passed. That case is
 *           common — it happens whenever a logger is reconfigured or two
 *           runs are concatenated. `coverage` compounded it by framing a
 *           deliberate 1 h pause as 88% data loss; it is gone, replaced by
 *           the location of the largest spacing jump, which the detector
 *           never reported and which is what a user actually needs.
 *           Second, the cluster detector never worked. It required
 *           h_prev < 0.1*h_avg && h_curr > 10*h_avg, i.e. neighbouring
 *           spacings differing by more than 100x while straddling the
 *           global mean, so an isolated 200x gap scored 0, dense blocks
 *           needed 10+ points each (five blocks of five scored 0 even with
 *           a 1000x jump), and only trailing edges counted, so k blocks
 *           yielded k-1. Across the whole scenario set used to plan this it
 *           fired once. It had zero test coverage. Removed outright rather
 *           than repaired: max_jump answers the same question better, with
 *           no global reference and with a location attached. New fields:
 *           max_jump, max_jump_x, n_jumps, regime_shift, multi_regime.
 *           A record is "mixed regimes" only when the median spacing of its
 *           first half differs from that of its second AND a real local jump
 *           exists. Both halves matter: a smoothly graded mesh shifts by a
 *           factor of 13000 between halves while no two neighbouring
 *           spacings differ by more than 1.1x, and calling that a regime
 *           change would be wrong. The first design of this used the
 *           fraction of spacings sitting at the base period instead, and was
 *           discarded on measurement — it cannot separate heavy scattered
 *           loss from a rate change (random 40% loss puts 55% of spacings at
 *           the base period, a rate change 50%: an overlap, not a margin).
 *           The half-to-half median ratio is 1.00 for dropouts of any
 *           severity against 10.00 for a rate change. Caught by the
 *           regression test from 005, which is why that test existed.
 *           Output is now exactly one characterization per grid, printed
 *           only when there is something to say: clean and smoothly
 *           non-uniform grids print nothing. Verified across ten scenarios
 *           plus both real data files; all.dat still reports 58 missing in
 *           48 gaps. Five tests (grid_analysis 12 -> 17, suite 133 -> 138).
 *           No new allocation — the existing scratch array is reused for the
 *           half medians. Zero leaks.
 * v5.11.51: Detect missing samples in a nominally regular grid.
 *           Fixed-rate logging drops samples, and when it does the gap left
 *           behind is an integer multiple of the base period: one lost
 *           sample gives 2*h_base, three consecutive give 4*h_base. Neither
 *           existing signal can see that. CV measures overall spread — the
 *           repo's own all.dat sits at 0.078, below every warning threshold
 *           in the file — and the cluster detector needs neighbouring
 *           spacings differing by more than 100x while straddling h_avg
 *           (an AND of h_prev < 0.1*h_avg and h_curr > 10*h_avg), whereas a
 *           dropout gap is 2x. Measured before the change, `-g` on all.dat
 *           reported "Grid is nearly uniform - standard methods work well,
 *           Detected clusters: 0" for a record missing 58 samples across 48
 *           gaps; pt.dat is missing 1106 samples, 20.5% of the record, and
 *           was equally silent. Not a correctness bug — the smoothing was
 *           fine — but a user had no way to learn a fifth of the data was
 *           gone. New: h_base (MEDIAN spacing), integer_fit, n_gaps,
 *           n_missing, max_run, has_dropouts. Median rather than mean
 *           because the mean is contaminated by the very gaps being looked
 *           for (pt.dat: mean 1.258 s against a true period of 1.0 s);
 *           measured, the median holds up to 45% sample loss and collapses
 *           at 50% when it jumps to 2*h_base. integer_fit gates the report
 *           and separates a regular grid with holes from genuinely
 *           non-uniform data with a wide margin: 100% for regular-with-
 *           dropouts against 55% (Poisson), 49% (geometric grading) and 0%
 *           (alternating mesh), so the 0.90 threshold is not finely tuned.
 *           Two further measured limits, documented rather than papered
 *           over: clock jitter above +-15% under-reports, and an unbalanced
 *           alternating mesh whose spacings are exact multiples of each
 *           other is a false positive by construction (30x0.1 with 10x0.4
 *           reports 150 missing when nothing is). Deliberately additive:
 *           every existing field keeps its value, the only output change is
 *           one new `verbose >= 1` line printed solely when dropouts are
 *           found, and reliability_warning is NOT set — smooth.c:248 gates
 *           a whole output block on it, so setting it would make every
 *           normal run on data with a single dropped sample print a warning
 *           it does not print today. Verified: `-g` on all.dat differs from
 *           v5.11.50 by exactly one added line; a graded mesh prints none.
 *           Five tests (grid_analysis 7 -> 12, suite 128 -> 133); replacing
 *           the median with the mean fails two of them. First malloc inside
 *           analyze_grid beyond the struct itself — it sits after the main
 *           loop, so the non-monotonic early return cannot leak it (checked
 *           separately under valgrind). Zero leaks, 1054 allocs/1054 frees.
 * v5.11.50: Parse CRLF input files. The `-T` tokenizer at
 *           parser.c:105-115 treated ' ', '\t' and '\n' as separators but
 *           not '\r', so on a file with CRLF line endings the carriage
 *           return stayed attached to the last whitespace token. strtod
 *           then failed the `endptr != y_tok_end` check at :163 and every
 *           row was counted into skipped_nonnumeric, ending the run with
 *           "ERROR: No valid data points found" plus a message blaming the
 *           data ("non-numeric or NaN/Inf value in column 2") rather than
 *           the line endings. Total failure, not degradation: a 4286-row
 *           file produced zero points. Any file that had passed through a
 *           Windows editor or a Windows logger was rejected wholesale, and
 *           the diagnostic pointed the user at the wrong thing. The plain
 *           numeric branch never had the bug — it lists '\r' explicitly at
 *           :211-232 — which is exactly why this stayed hidden: the two
 *           tokenizers in the same function disagreed about one character.
 *           Fixed at the root rather than per-branch: the line terminator,
 *           CR included, is now stripped once after fgets (parser.c:70-81),
 *           before the comment strip and after the truncation check, which
 *           still needs the raw '\n'. One guard covers both tokenizers and
 *           any future one. Verified: output on LF input is byte-identical
 *           to v5.11.49 across all four methods plus -g and -T (14
 *           input/flag combinations, zero differing bytes); CRLF edge cases
 *           all parse — missing final newline, blank lines, full-line and
 *           inline comments, doubled CR, trailing whitespace. Three
 *           regression tests added (parser 14 -> 17, suite 125 -> 128); the
 *           two -T ones fail against the old parser with 0 rows instead of
 *           5, while the numeric one passes there by design, guarding that
 *           the shared strip did not disturb the branch that was already
 *           correct. Zero valgrind leaks.
 * v5.11.49: Make the Tikhonov `-d` output second order on
 *           non-uniform grids. compute_derivatives() used a difference that
 *           is symmetric in index but not in coordinate,
 *           (u_{i+1}-u_{i-1})/(x_{i+1}-x_{i-1}), whose Taylor expansion leaves
 *           a leading error (h_r - h_l)/2 * u'' that vanishes only for
 *           h_l == h_r. On an alternating 0.1/0.4 mesh with u = sin x and
 *           lambda = 0 (so y_smooth == y exactly and the differentiation error
 *           is isolated) theory and measurement agree to three figures:
 *           predicted 1.499e-01, measured 1.494e-01. This was the worst place
 *           for a uniform-grid assumption to hide: Savgol rejects above
 *           CV 0.05 and Butterworth above CV 0.15, and savgol.c's rejection
 *           message points the user at "-m 2 -l auto (Works correctly with
 *           non-uniform grids)" — so non-uniform data lands precisely here.
 *           The v5.11 line made the *penalty* grid-correct (weighted Gram
 *           matrix (D2)^T W D2, integral measure, v5.11.34); the derivative it
 *           emitted never got that care, leaving `-m 0` (polyfit, which
 *           differentiates its local fit analytically) more accurate than the
 *           method the tool recommends. Fix: all three formulas replaced by
 *           three-point weights from undetermined coefficients — chosen to
 *           reproduce u, u' and u'' exactly, so the error is O(h^2)u''' on any
 *           spacing — one-sided at each boundary, spacing-aware in the
 *           interior, same O(n) cost. On a uniform grid the interior weights
 *           degenerate to the classical -1/(2h), 0, +1/(2h). Measured:
 *           interior max error 1.494e-01 -> 6.642e-03 (22x); uniform-mesh
 *           boundary on sin(x+1) 1.10e-01 -> 7.76e-03 (14x), i.e. the boundary
 *           gain is unconditional — the old two-point one-sided difference was
 *           first order even on uniform grids. (Where u'' happens to be 0 at
 *           an endpoint the old formula's leading term vanishes by coincidence
 *           and the new one is not better there; that is the only case it does
 *           not win.) Printed output: across 2970 uniform-grid interior points
 *           spanning spacings 0.01-3.7 and amplitudes 1e-3-1e3, zero lines
 *           differ under %12.8lG — the new interior form is not bitwise
 *           identical (~70 ULP) but is identical to 8 significant digits, so
 *           the byte-identity release check still holds for interior points.
 *           The two boundary points do change, by design, on every input.
 *           No h_l == h_r fast path: floating-point accumulation means a
 *           nominally uniform mesh does not have exactly equal spacings, so
 *           such a branch would fire unpredictably. build_band_matrix() and
 *           compute_functional() are untouched — this changes only the
 *           reported derivative, not what is minimized. Two regression tests
 *           (tikhonov 25 -> 27, suite 123 -> 125) assert exactness for a
 *           quadratic on alternating and uniform meshes over all points
 *           including both boundaries; both FAIL against the old
 *           implementation (-1.7 and -1.25 vs expected -2). Zero valgrind
 *           leaks.
 * v5.11.48: Report malformed timestamp rows in `-T` mode instead
 *           of dropping them silently. The `else` branch at parser.c:141
 *           (now :142) fired whenever the timestamp column was present but
 *           the second whitespace token was missing — a bare `2026-01-01`,
 *           a stray label, or a truncated final line — and simply
 *           `continue`d without incrementing any counter or printing
 *           anything. It was the only unaccounted row-drop in the file; its
 *           sibling non-numeric-y drop (parser.c:161-164) already counted
 *           into skipped_nonnumeric. Invisible and misleading: with no
 *           message at all, the resulting gap in the x grid was reported
 *           back to the user as a spurious "# WARNING: Significant spacing
 *           variation detected: CV = 0.38", blaming the user's sampling for
 *           a hole the parser itself created. Fix: a new
 *           skipped_malformed_ts counter, incremented at the drop site and
 *           reported on stdout as "# Skipped N data row(s) with malformed
 *           timestamp in column M", mirroring skipped_nonnumeric's existing
 *           convention and output stream. Which rows get dropped is
 *           unchanged — a truncated final line is common in real data and
 *           must stay tolerated; only whether the user is told changed.
 *           Verified: three program outputs (-T, Tikhonov, polyfit -d) on
 *           clean input byte-identical to v5.11.47; 116 tests pass; zero
 *           valgrind leaks. The end-to-end regression guard for this path
 *           (a malformed-timestamp fixture exercising the parser) lands
 *           with the timestamp parser tests in a follow-up plan.
 * v5.11.47: Fix multi-line grid warnings corrupting the output
 *           table. print_grid_analysis() prefixed only the first line of
 *           warning_msg with the caller's "# ", so every continuation line
 *           ("Adaptive methods may improve results.", "Consider resampling
 *           data to a more uniform grid.", ...) was written to stdout as
 *           bare text in the middle of the data table, breaking any numeric
 *           consumer. All three CV bands and the cluster warning are
 *           multi-line, so this hit every non-uniform input above CV 0.2.
 *           The fix walks warning_msg by line and emits the prefix on each;
 *           only line 1 carries the "WARNING: " label, so line-1 output is
 *           byte-identical to v5.11.46. Verified over CV 0.00-4.63 plus the
 *           clustered case: zero unprefixed lines. Regression guard lives in
 *           the run-smooth driver (stdout table integrity check).
 * v5.11.46: Ponytail audit v5.11.43 — closes the report. Two
 *           findings applied, one rejected on measurement.
 *           (19) tikhonov: the two dcopy_ calls were plain array copies ->
 *           memcpy; drops the extern BLAS declaration and the `inc` variable.
 *           Bit-identical by construction (no arithmetic). Note this removes
 *           the only pure-BLAS call in the tree but not the -lblas dependency:
 *           LAPACK pulls BLAS in regardless.
 *           (16) Makefile: removed the MEMCHECK variable (a commented-out
 *           electric-fence hook that always expanded to nothing) and the
 *           `memcheck` target, which built debug and then only echoed the
 *           valgrind command line for the user to retype. `make test-valgrind`
 *           already runs valgrind for real with --error-exitcode=1.
 *           (15) REJECTED: replacing savgol's power() with pow(). The audit
 *           assumed pow() is exact for these small integer exponents; measured
 *           over the range savgol computes (j^i for j in the window, i up to
 *           2*poly_degree), 906 of 16875 values differ in the last ULP. Those
 *           moments build the normal-equations matrix of a Vandermonde system,
 *           badly conditioned by construction, so the noise would propagate
 *           into the coefficients. Only degree 11+ in a wide window is
 *           affected, but that is the worst-conditioned case. The exact
 *           integer loop stays.
 *           All 21 audit findings are now resolved: 18 applied across
 *           v5.11.44-46, 3 rejected (7 documented -g output, 15 above, 21
 *           replacing Unity with assert.h).
 *           Six program outputs re-verified byte-for-byte against v5.11.45.
 *           116 tests pass, zero valgrind leaks.
 * v5.11.45: Ponytail audit v5.11.43 — the `shrink:` and `yagni:`
 *           findings (4, 5, 6, 8, 10, 11, 13, 17, 18). Pure refactoring: every
 *           one of the nine is behaviour-preserving, verified by diffing nine
 *           program outputs (all four methods, -g, -T, -h) byte-for-byte
 *           against v5.11.44 — all identical. Findings 7 (`-g` output fields)
 *           and 21 (Unity) were deliberately left in place.
 *           Butterworth: check_pole_stability() now calls max_pole_radius()
 *           instead of repeating its scan (the biquad index is dropped from the
 *           message); butterworth_filtfilt() calls run_filtfilt_trial() instead
 *           of repeating the IC/pad/cascade/reverse sequence, which also
 *           retires the y_work buffer; estimate_memory_usage() inlined and the
 *           no-op NULL guards in free_butterworth_result() dropped.
 *           Savgol: the two 39-line boundary blocks became one loop. An
 *           asymmetric window still spans window_size points and the central
 *           loop is finished by then, so c_func/c_deriv are reused as scratch —
 *           two allocations per boundary point and a duplicated cleanup path
 *           are gone.
 *           Polyfit: fill_boundary_fallback_left/right deleted. They were
 *           reachable only when dgelss fails at exactly the first or last
 *           interior window; that case now uses raw y with zero derivative,
 *           the same substitution the inline fallback already makes for every
 *           other window (the boundary points cannot be left at their calloc'd
 *           zeros, so the replacement is two short loops, not nothing).
 *           Timestamp: convert_timestamps_to_relative() hands its working x
 *           array straight to the caller instead of allocating a second one and
 *           memcpy-ing into it (over-allocated by the dropped rows, harmless
 *           since only ctx->n entries are read), and its three copy-pasted
 *           cleanup blocks collapse into one `goto fail`. original_timestamps
 *           switched to calloc so that path is safe from every failure point.
 *           Grid analysis: analyze_grid() loses the store_spacings parameter
 *           (never once passed 1), the struct loses `spacings` and
 *           `n_intervals`, and print_grid_analysis() loses its verbose>=2
 *           branch (never once reached). 73 call sites updated, nearly all in
 *           tests. `-g` output is unchanged.
 *           smooth.c: usage() folded into its only caller, help().
 *           tikhonov.c: compute_gcv_score_robust() loses its `verbose`
 *           parameter (both call sites passed 1).
 *           Net -181 lines; heap allocations across the test suite drop from
 *           1148 to 911. 116 tests pass, zero valgrind leaks.
 * v5.11.44: Ponytail audit v5.11.43 — all seven `delete:` findings
 *           (dead code and duplicated layers; see doc/ponytail-audit-v5.11.43.md).
 *           (1) Tikhonov: removed find_lambda_lcurve() and the n>20000 branch it
 *           was reachable from. That branch duplicated the GCV sweep with a
 *           hardcoded 12-lambda list, ran the L-curve on the same list, then
 *           arbitrated between two answers. One log-spaced 13-point GCV sweep
 *           over [1e-8, 1e0] now serves every n. BEHAVIOUR CHANGE: for
 *           n > 20000 the selected lambda may differ by up to one grid step;
 *           the search range and the n<=5000 refinement rule are unchanged.
 *           (2) Removed the decomment module (decomment.c/h, 149 lines) and its
 *           tmpfile() copy of the whole input. parse_input() now strips '#'
 *           comments itself (strchr, one place, both modes) and its existing
 *           blank-line checks drop comment-only lines, so they still never count
 *           as skipped data rows. smooth.c opens the file (or uses stdin)
 *           directly. Two improvements fall out: error messages now report the
 *           real input line number instead of the post-decomment temp-file line,
 *           and a data line whose truncation falls inside a trailing comment is
 *           parsed instead of rejected. New test parser_comments_are_stripped.
 *           (3) Deleted untracked scratch programs debug_test.c,
 *           test_actual_results.c, test_debug_edge.c and its compiled binary.
 *           (9) Deleted unused test helper create_grid_with_cv().
 *           (12) Deleted unused free_parse_result().
 *           (14) Deleted savgol factorial(): savgol_coefficients is static and
 *           only ever called with deriv_order 0 or 1, both giving 1.0. The
 *           parameter validation now rejects deriv_order > 1 so the hardcoded
 *           RHS cannot silently go wrong later.
 *           (20) Deleted unused macro BUTTERWORTH_NUM_COEFFS.
 *           Net -476 lines, one module fewer. 116 tests pass, zero valgrind leaks.
 * v5.11.43: Code audit v5.11.41 fix A2 — parse_timestamp accepted
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
#define VERSION "5.11.53"
#define REVDATE "2026-07-28"
