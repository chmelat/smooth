/*
 * revision.h - define the version number and date
 *
 * Version History
 * ---------------
 * v5.11.17 (current): Audit C6 — `make test-valgrind` now fails (exit 1) on
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
#define VERSION "5.11.17"
#define REVDATE "2026-04-28"
