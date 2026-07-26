# Ponytail audit v5.11.43 — over-engineering scan

**Date:** 2026-07-26
**Scope:** whole tree (14561 lines: 9 source modules + tests)
**Focus:** over-engineering and complexity only. Correctness bugs, security and
performance are explicitly out of scope — route those to a normal review pass.

This report **lists** findings. **Status:** the seven `delete:` findings (1, 2, 3, 9, 12, 14, 20) were applied
in **v5.11.44**; the `shrink:` and `yagni:` findings (4, 5, 6, 8, 10, 11, 13,
17, 18) in **v5.11.45**. Findings 7 and 21 were reviewed and deliberately
rejected. Still open: `stdlib:` 15 and 19, `native:` 16.

Tags: `delete:` dead code / speculative feature. `stdlib:` hand-rolled thing the
standard library ships. `native:` code doing what the platform already does.
`yagni:` abstraction with one implementation, config nobody sets, layer with one
caller. `shrink:` same logic, fewer lines.

---

## Ranked findings

| # | Tag | Location | Cut | Stav |
|---|-----|----------|-----|------|
| 1 | delete | tikhonov.c:333-421, 450-470 | -105 | DONE v5.11.44 |
| 2 | delete | decomment.c, decomment.h, smooth.c:196-206 | -145 | DONE v5.11.44 |
| 3 | delete | untracked scratch files in repo root | -178 | DONE v5.11.44 |
| 4 | shrink | polyfit.c:101-143, 353-362 | -43 | DONE v5.11.45 |
| 5 | shrink | savgol.c:341-420 | -35 | DONE v5.11.45 |
| 6 | yagni | grid_analysis.c:61-68, 94-96, 245-253 | -25 | DONE v5.11.45 |
| 7 | yagni | grid_analysis.c:139-142, 190-213 | -25 | zamítnuto (výstup -g zůstává) |
| 8 | shrink | timestamp.c:126-217 | -25 | DONE v5.11.45 |
| 9 | delete | tests/grid_helpers.c:97-124 | -24 | DONE v5.11.44 |
| 10 | shrink | butterworth.c:182-210 | -15 | DONE v5.11.45 |
| 11 | shrink | butterworth.c:676-698 | -15 | DONE v5.11.45 |
| 12 | delete | parser.c:312-322, parser.h:46 | -13 | DONE v5.11.44 |
| 13 | shrink | butterworth.c:97-101, 696-698, 720-731 | -12 | DONE v5.11.45 |
| 14 | delete | savgol.c:46-54, 143 | -10 | DONE v5.11.44 |
| 15 | stdlib | savgol.c:35-43 | -9 | open |
| 16 | native | Makefile:31-33, 48, 67-69 | -8 | open |
| 17 | yagni | smooth.c:415-421 | -7 | DONE v5.11.45 |
| 18 | yagni | tikhonov.c:271-330 | -4 | DONE v5.11.45 |
| 19 | stdlib | tikhonov.c:32, 224, 247 | -3 | open |
| 20 | delete | butterworth.h:16 | -1 | DONE v5.11.44 |
| 21 | yagni | tests/unity.c, unity.h, unity_internals.h | -4595 | zamítnuto (Unity zůstává) |

**Pozn.:** řádky a čísla řádků v sekcích níže odpovídají stavu v5.11.43. Po
provedení `delete:` nálezů se čísla řádků v `tikhonov.c`, `savgol.c`,
`parser.c` a `smooth.c` posunula — u zbývajících nálezů je nutné je dohledat
znovu.

---

### 1. `delete:` Tikhonov L-curve method and the `n > 20000` branch

`find_lambda_lcurve()` is reachable only from the `n > 20000` branch of
`find_optimal_lambda_gcv()`. That branch duplicates the GCV sweep with a
hardcoded 12-lambda list, then runs the L-curve on the same list and arbitrates
between two answers ("GCV and L-curve disagree - using more conservative
choice"). Two selection criteria, one of which is unreachable for the data sizes
this tool is normally used on.

**Replacement:** the standard log-spaced GCV sweep in the `else` branch, for all
`n`. The `n <= 5000` refinement guard already scales the work.

**Location:** `tikhonov.c:333-421` (function), `tikhonov.c:450-470` (branch).

---

### 2. `delete:` The whole `decomment` module

A separate module plus a `tmpfile()` round-trip of the entire input, purely to
strip `#` comments and blank lines, before the parser reads it back line by
line. The parser already reads with `fgets` and already skips blank lines
(`parser.c:89`, `parser.c:187`).

**Replacement:** two lines at the top of the `parse_input` loop —

```c
char *hash = strchr(line, '#');
if (hash) *hash = '\0';
```

— and `smooth.c` opens the file (or uses `stdin`) directly. Removes a module, a
temporary file, and a full copy of the input.

**Location:** `decomment.c` (116), `decomment.h` (33), `smooth.c:196-206`.
Also drops `decomment.c`/`decomment.h` from `SRC`/`HEAD` in the Makefile.

---

### 3. `delete:` Untracked scratch programs in the repo root

`debug_test.c`, `test_actual_results.c` and `test_debug_edge.c` are three
variants of the same throwaway main() with three copies of a local
`create_uniform_grid()`. `test_debug_edge` is a 72 kB compiled binary sitting in
the tree. `savgol.c~`, `smooth.c~`, `makefile~`, `revision.h~` are stale editor
backups (gitignored, but present on disk).

**Replacement:** nothing — `tests/` covers these cases.

---

### 4. `shrink:` Polyfit boundary fallbacks

`fill_boundary_fallback_left()` and `fill_boundary_fallback_right()` are two
near-identical 21-line functions, reachable only when `dgelss` fails at exactly
the first or the last interior window (any other failure is already handled
inline by the `fallback_count` path at `polyfit.c:281-287`).

**Replacement:** one function parameterised by direction, or drop straight to
`y[k]` with `dy = 0` as the inline path already does.

**Location:** `polyfit.c:101-143`, call sites `polyfit.c:353-362`.

---

### 5. `shrink:` Savgol left/right boundary loops

`savgol.c:341-379` and `savgol.c:382-420` are the same 39 lines twice —
allocate, compute coefficients, convolve, free, with an identical five-line
error-cleanup block. The only difference is how `left_pts` is derived.

**Replacement:** one loop over the `2*offset` boundary indices with

```c
int left_pts = (i < offset) ? i : window_size - 1 - (n - 1 - i);
```

---

### 6. `yagni:` `store_spacings` / `spacings` / `n_intervals` / `verbose >= 2`

`analyze_grid()` is called with `store_spacings = 0` at every call site in the
program and in all tests. The only consumer of `spacings` is the `verbose >= 2`
branch of `print_grid_analysis()`, and no caller ever passes verbose above 1
(`smooth.c` passes 0 and 1). `n_intervals` is used only inside that dead branch.

**Replacement:** drop the parameter, the struct field, the allocation, and the
branch.

**Location:** `grid_analysis.c:61-68, 94-96, 245-253`; `grid_analysis.h:25, 34,
39`. Touches `tests/test_grid_analysis.c:276` (the `n_intervals` assertion).

---

### 7. `yagni:` `uniformity_score` and `get_grid_recommendation()`

`uniformity_score` is `exp(-2 * CV)` — a monotone restatement of CV that carries
no new information. `get_grid_recommendation()` then re-derives CV bands a third
time from that score (thresholds 0.8 / 0.5 / 0.2), alongside the
`THRESH_CV_*` ladder already in the same file and the per-method thresholds in
savgol/butterworth/tikhonov.

**Replacement:** branch on `cv` directly at the single print site.

**Location:** `grid_analysis.c:139-142` (score), `grid_analysis.c:190-213`
(recommendation), `grid_analysis.h:27, 45`.

**Note:** this touches a documented policy point — CLAUDE.md's "Grid uniformity
philosophy" table. Update it in the same change.

---

### 8. `shrink:` Timestamp array plumbing

`convert_timestamps_to_relative()` fills `x_temp` (size `n`), then mallocs a
second array of size `valid_count`, memcpys into it, and frees `x_temp` — plus
an extra failure branch for the second allocation. Three copy-pasted cleanup
blocks free the same four things.

**Replacement:** hand `x_temp` straight to `*x_out` (over-allocated by the
number of dropped rows, which is harmless — the caller only reads `ctx->n`
elements), and route all failures through one `goto fail`.

**Location:** `timestamp.c:126-217`.

---

### 9. `delete:` `create_grid_with_cv()`

Zero callers. Implements a trial-and-error CV fitter with an "empirical
relationship" fudge factor (`target_cv * 2.5`) and returns the CV it happened to
land on.

**Location:** `tests/grid_helpers.c:97-124`, `tests/grid_helpers.h`.

---

### 10. `shrink:` `check_pole_stability()` re-implements `max_pole_radius()`

Both scan `NUM_BIQUADS` sections for the largest `biquad_pole_radius()`.
`max_pole_radius()` is defined ten lines above.

**Replacement:** call `max_pole_radius()`; drop `worst_section` from the
diagnostic (with 2 sections it adds nothing actionable).

**Location:** `butterworth.c:182-210`.

---

### 11. `shrink:` Main filtfilt path duplicates `run_filtfilt_trial()`

`butterworth_filtfilt()` repeats the exact sequence already in
`run_filtfilt_trial()`: compute ICs, pad, cascade, reverse, cascade, reverse,
extract.

**Replacement:** call `run_filtfilt_trial()` after validation. The `sections`
and `pad_len` recomputation inside it is O(1) compared to the filter passes.

**Location:** `butterworth.c:676-698` vs `butterworth.c:452-475`.

---

### 12. `delete:` `free_parse_result()`

Zero callers — `smooth.c` frees `pr.x`, `pr.y` and `ts_ctx` by hand in its
`cleanup:` block (`smooth.c:363-368`), and `parse_input()` already frees its own
partial allocations on failure.

**Location:** `parser.c:312-322`, `parser.h:46`.

---

### 13. `shrink:` Butterworth small change

- `estimate_memory_usage()` is a `static inline` function with one caller
  (`butterworth.c:634`) — inline the expression.
- `free_butterworth_result()` guards each `free()` with a NULL check;
  `free(NULL)` is a no-op.
- The extraction loop at `butterworth.c:696-698` is a `memcpy` — the same file
  already uses `memcpy` for exactly this at line 472.

---

### 14. `delete:` Savgol `factorial()`

Called once, as `factorial(deriv_order)`, and `savgol_coefficients()` is only
ever called with `deriv_order` 0 or 1 (`savgol.c:298, 306, 354-355, 395-396`).
Both return 1.0.

**Replacement:** `B[j] = 1.0;` at `savgol.c:143`.

---

### 15. `stdlib:` Savgol `power()`

Hand-rolled integer power loop. `pow()` from `math.h` is already available and
`-lm` is already linked. Exponents here are small non-negative integers, where
`pow()` is exact on any IEEE-754 implementation.

**Location:** `savgol.c:35-43`, call sites `savgol.c:133, 161`.

---

### 16. `native:` Makefile memcheck scaffolding

`MEMCHECK` is a commented-out electric-fence hook threaded into `LIB`, and the
`memcheck` target builds debug and then only echoes the command the user should
type. `test-valgrind` already runs valgrind for real, with
`--error-exitcode=1`.

**Location:** `Makefile:31-33, 48, 67-69` (and the `memcheck` lines in `help`).

---

### 17. `yagni:` `usage()` in smooth.c

Single caller: `help()`, defined immediately below it (`smooth.c:469`).

---

### 18. `yagni:` `compute_gcv_score_robust(..., verbose)`

All four call sites pass `1`.

**Location:** `tikhonov.c:271-330`, call sites at 456, 479, 495.

---

### 19. `stdlib:` BLAS `dcopy_` for two plain array copies

`dcopy_(&n, y, &inc, b, &inc)` and the mirror call are `memcpy`. This is the
only pure-BLAS call in the tree (LAPACK's `dpbsv`/`dposv`/`dgelss` still need
`-llapack -lblas`, so no dependency is removed — just the `extern` declaration
and the `inc` variable).

**Location:** `tikhonov.c:32, 224, 247`.

---

### 20. `delete:` `BUTTERWORTH_NUM_COEFFS`

Defined, never referenced.

**Location:** `butterworth.h:16`.

---

### 21. `yagni:` Vendored Unity (4595 lines)

`unity.c` + `unity.h` + `unity_internals.h` are 4595 lines supporting what the
suite actually uses: `TEST_ASSERT_DOUBLE_WITHIN`, a handful of int/pointer
assertions, `RUN_TEST`, and a pass/fail tally.

**Replacement:** `assert.h` plus a ~20-line runner.

**Ranked last deliberately.** This is the one finding where the framework may be
earning its keep: per-test isolation via `setjmp`, "continue after failure"
reporting, and a stable output format. If the suite's value is in *seeing which
of 113 tests failed in one run*, keep Unity. If it is in *failing the build*,
`assert.h` is the lazier tool.

---

## Net

```
net: -702 lines, -0 deps possible
     (-5297 lines and -1 vendored dep if Unity goes too)
```

No third-party runtime dependency can be dropped: LAPACK is load-bearing in
three of the four methods, and BLAS comes along with it regardless of finding 19.

## Notes on applying

- Findings 1, 2, 6, 7 change observable behaviour or documented policy — they
  need a `revision.h` bump and a CLAUDE.md update (grid uniformity table for 6/7,
  module list for 2).
- Findings 6 and 9 touch `tests/`; re-run `make test` and `make test-valgrind`.
- Everything else is behaviour-preserving.
