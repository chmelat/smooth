# Smooth v5.11.0 -- Critical Code Analysis Report

**Date:** 2026-02-15
**Scope:** Complete source code review of all modules
**Files analyzed:** `smooth.c`, `tikhonov.c`, `butterworth.c`, `savgol.c`, `polyfit.c`, `grid_analysis.c`, `timestamp.c`, `decomment.c` and all headers

---

## Executive Summary

The codebase is well-structured with clean modular architecture, consistent error handling (goto cleanup pattern), and careful memory management. The analysis identified 15 issues: 2 high-severity bugs in the timestamp module (both now fixed), 3 medium-severity issues affecting numerical accuracy and input handling, and 10 low-severity issues related to API consistency and code hygiene.

---

## Issues Found

### HIGH severity

#### H1. Memory leak in timestamp cleanup -- FIXED

**File:** `smooth.c:407-412`
**Fixed in:** commit `24311e8` (saved original count before overwriting `n`)

After converting timestamps to relative time, the original `n` was overwritten before freeing the original string array. The cleanup loop used the new `n` instead of the original count, leaving unfreed strings.

---

#### H2. DST transitions corrupt relative timestamps -- FIXED

**File:** `timestamp.c:74-77`
**Fixed in:** v5.11.0+ (replaced `mktime()` with `timegm()`)

`mktime()` interpreted timestamps as **local time**. If data spanned a DST transition (e.g., spring forward), the relative time array contained a spurious $\pm 3600$ s jump. Fix: `timegm()` interprets `struct tm` as UTC -- no DST transitions possible. Added `#define _DEFAULT_SOURCE` for portability and a DST-invariance unit test (`test_parse_timestamp_dst_invariant`, 103 tests total).

---

### MEDIUM severity

#### M1. L-curve invalid point detection uses ambiguous sentinel

**File:** `tikhonov.c:482-484, 494`

When `tikhonov_smooth()` returns NULL, the L-curve code marks the point with `rss_vals[i] = 0.0`. Later, valid points are skipped if their value equals `0.0`:

```c
if (rss_vals[i-1] == 0.0 || rss_vals[i] == 0.0 || ...) {
    continue;  // skip
}
```

Since `rss_vals` stores $\log(\text{data\_term})$, a data\_term of exactly $1.0$ yields $\log(1.0) = 0.0$ -- a perfectly valid value that would be falsely skipped.

**Fix:** Use a separate boolean validity array, or use `NAN` as the sentinel (and check with `isnan()`).

---

#### M2. Stdin input bypasses comment filtering

**File:** `smooth.c:181-184`

When reading from stdin, `decomment()` is not called. Lines containing `#` characters are passed directly to the numeric parser. While the parser will likely fail to parse these lines and skip them silently, this is inconsistent with file input behavior and can cause confusing results.

**Fix:** Either apply `decomment`-style filtering inline during parsing, or document that stdin input must not contain comments.

---

#### M3. GCV eigenvalue approximation inaccurate for near-threshold grids

**File:** `tikhonov.c:410-416`

The GCV trace computation uses analytical eigenvalues derived for a **uniform** grid:

$$\text{ev}_k = \left(\frac{4 \sin^2(\theta_k/2)}{h_{\text{avg}}^2}\right)^2$$

This formula is applied for all grids with $\text{CV} < 0.15$ (the AVERAGE discretization branch). For grids with $\text{CV}$ in the range $0.10$--$0.14$, the eigenvalue approximation can be significantly off, leading to suboptimal $\lambda$ selection.

**Impact:** Automatic lambda selection may over-smooth or under-smooth data with moderately non-uniform grids.

---

### LOW severity

#### L1. Unnecessary validation of window parameters for Tikhonov/Butterworth

**File:** `smooth.c:166-169`

The window size validation (`sp` must be odd and $\geq 3$) runs unconditionally, even when the selected method (Tikhonov, Butterworth) does not use window parameters. With default values this is harmless, but if a user specifies `-m 2 -n 4`, they get a confusing error about window size.

**Fix:** Move validation inside the relevant method branches.

---

#### L2. `estimate_cutoff_frequency()` is a stub

**File:** `butterworth.c:206-210`

```c
double estimate_cutoff_frequency(const double *x, const double *y, int n)
{
    (void)x; (void)y; (void)n;
    return 0.2;
}
```

When the user specifies `-f auto`, they receive a hardcoded value of $0.2$ without any indication that auto-selection is not implemented. The output says "Auto-selected cutoff frequency: fc = 0.2000" which is misleading.

---

#### L3. Inconsistent CV thresholds across modules

Threshold constants are defined locally in each module without a shared source of truth:

| Module | Constant | Value | Purpose |
|--------|----------|-------|---------|
| `tikhonov.c` | `CV_THRESHOLD` | 0.15 | Discretization switch |
| `savgol.c` | `UNIFORMITY_CV_THRESHOLD` | 0.05 | Input rejection |
| `butterworth.c` | `UNIFORMITY_CV_THRESHOLD` | 0.15 | Input rejection |
| `butterworth.c` | `UNIFORMITY_CV_WARNING` | 0.05 | Warning |
| `grid_analysis.c` | `THRESH_CV_UNIFORM` | 0.01 | `is_uniform` flag |

If the uniformity semantics change, multiple files must be updated independently.

---

#### L4. Missing `const` qualifiers on input arrays

Functions `tikhonov_smooth()`, `polyfit_smooth()`, and `savgol_smooth()` accept `double *x` and `double *y` parameters that are never modified. These should be `const double *` for type safety and to enable compiler optimizations.

---

#### L5. Savitzky-Golay boundary derivatives use global `h_avg`

**File:** `savgol.c:370, 429`

Derivative coefficients at boundary points are divided by `h_avg` (the global average spacing), not the local spacing. For grids with $0.01 < \text{CV} < 0.05$ (accepted by Savitzky-Golay), this introduces a small but unnecessary error at the edges.

---

#### L6. Memory leaks on method failure in `smooth.c`

**File:** `smooth.c:462-465, 517-519, 571-572, 612-613`

When a smoothing method returns NULL, the program calls `exit(EXIT_FAILURE)` without freeing `grid_info`, `x`, `y`, or `ts_ctx`. The OS reclaims memory on exit, so this is not a runtime leak, but it clutters valgrind output and is inconsistent with the careful cleanup elsewhere.

---

#### L7. Timestamp parser is fragile for edge cases

**File:** `smooth.c:217-252`

The `sscanf`-based timestamp parser uses a fragile multi-branch strategy to handle space vs. `T` separators. The logic for `parsed == 3` with `time_str[0] == '.'` (line 242) attempts to handle cases where `sscanf` splits subseconds, but this can misparse unusual inputs. A dedicated parsing function (already in `timestamp.c`) would be more robust here.

---

#### L8. `polyfit.c` rebuilds SVD at every point

**File:** `polyfit.c:257-327`

For each central point, a full SVD decomposition (`dgelss_`) is computed. Adjacent windows share $n-1$ data points but the Vandermonde matrix is rebuilt and decomposed from scratch. This is $O(n \cdot p^3)$ where caching could reduce it to $O(n \cdot p^2)$ for uniform grids. Current approach is correct but suboptimal for large datasets.

---

#### L9. Header example code is outdated

**File:** `tikhonov.h:101-104`

The usage example in the header shows function signatures without `grid_info` parameter:

```c
TikhonovResult *result = tikhonov_smooth(x, y, n, 0.1);       // wrong
double optimal_lambda = find_optimal_lambda_gcv(x, y, n);      // wrong
```

These calls are missing the `GridAnalysis*` parameter added in later versions.

---

#### L10. No protection against `lambda == 0` with auto-lambda

**File:** `tikhonov.c:533`

`find_optimal_lambda_gcv()` initializes `best_lambda = 0.01`, but the search range is $[10^{-8}, 10^0]$. If all GCV scores are $10^{20}$ (all calls fail), the function returns `0.01` without indication of failure. A more robust approach would return a status code or log a warning.

---

## Summary Table

| ID | Severity | Module | Description |
|----|----------|--------|-------------|
| H1 | ~~HIGH~~ FIXED | `smooth.c` | Memory leak in timestamp string cleanup |
| H2 | ~~HIGH~~ FIXED | `timestamp.c` | DST transitions corrupt relative time |
| M1 | MEDIUM | `tikhonov.c` | L-curve sentinel value collides with valid data |
| M2 | MEDIUM | `smooth.c` | Stdin bypasses comment filtering |
| M3 | MEDIUM | `tikhonov.c` | GCV eigenvalues inaccurate for $0.10 \leq \text{CV} < 0.15$ |
| L1 | LOW | `smooth.c` | Irrelevant window validation for Tikhonov/Butterworth |
| L2 | LOW | `butterworth.c` | Auto-cutoff returns hardcoded value |
| L3 | LOW | multiple | CV thresholds not centralized |
| L4 | LOW | multiple | Missing `const` on input array parameters |
| L5 | LOW | `savgol.c` | Boundary derivatives use `h_avg` instead of local step |
| L6 | LOW | `smooth.c` | Missing cleanup before `exit()` on method failure |
| L7 | LOW | `smooth.c` | Fragile timestamp parser logic |
| L8 | LOW | `polyfit.c` | Redundant SVD decomposition per window |
| L9 | LOW | `tikhonov.h` | Outdated function signatures in example |
| L10 | LOW | `tikhonov.c` | No failure indication from GCV search |

---

## Overall Assessment

The codebase demonstrates solid engineering practices: consistent modular design, proper use of LAPACK for numerical stability, and defensive memory management. Both high-severity issues (H1, H2) have been fixed. The most actionable remaining issue is M1 (simple sentinel change). The remaining issues are either cosmetic or affect edge cases with limited practical impact.
