# CLAUDE.md

Design and maintenance notes for working on `smooth`. User-facing documentation
(installation, CLI options, method tutorials, examples) lives in `README.md` —
do not duplicate it here.

**Current version:** see `revision.h`. The full version history is the comment
block at the top of `revision.h`.

**Audit reports, analysis writeups, comparison studies:** `doc/`. Reference
these from commit messages rather than duplicating their content here.

## Documentation Guidelines

### README.md equation format

The README uses **GitHub LaTeX math syntax** for all mathematical content:

- Inline: `$...$` (e.g. `$\lambda$`, `$\|y - u\|^2$`)
- Display: `$$...$$` for standalone equations
- Matrices: `\begin{pmatrix}...\end{pmatrix}` inside `$$...$$`

Do not use plain-text approximations of math in code blocks. Code blocks are
reserved for code, CLI examples, and program output.

### README.md character restrictions

PDF generation uses DejaVu fonts which lack box-drawing and decorative glyphs.
Restrictions:

- **Box-drawing:** use `|-`, `-`, `|`, `+-`, `=` instead of `├ ─ │ └ ━`
- **Checkmarks/warnings:** use `[OK]`, `[X]`, `[WARNING]` instead of `✓ ✗ ⚠️`
- **Arrows:** `→` is allowed in plain text; inside math use `\to` / `\rightarrow`

## Architecture

### Modular structure

```
smooth.c              # Main program, CLI parsing, I/O, output formatting
├─ polyfit.c/h        # Polynomial fitting (local least squares)
├─ savgol.c/h         # Savitzky-Golay filter (pre-computed convolution)
├─ tikhonov.c/h       # Tikhonov regularization (global variational)
├─ butterworth.c/h    # Butterworth filter (frequency-domain)
├─ grid_analysis.c/h  # Grid uniformity analysis (shared utility)
├─ timestamp.c/h      # Timestamp parsing for `-T` mode (UTC via timegm())
└─ decomment.c/h      # Input comment removal (full-line and inline `#`)
```

### Design principles

1. **One grid analysis, shared results.** `analyze_grid()` runs once at startup
   in `smooth.c`; the resulting `GridAnalysis*` is passed to every method.
   Methods do not re-analyze the grid.
2. **Result structures.** Each method returns a `*Result` struct (`PolyfitResult`,
   `SavgolResult`, `TikhonovResult`, `ButterworthResult`) with `y_smooth`,
   `y_deriv`, `n`, plus method-specific diagnostics.
3. **Caller owns memory.** Methods allocate; the caller frees via the matching
   `free_*_result()`. The same pattern applies to `GridAnalysis`.
4. **Grid-aware methods.** Methods inspect `GridAnalysis*` and either adapt
   (Tikhonov) or reject (Savgol) based on uniformity — the policy is in the
   method, not in `smooth.c`.
5. **No cross-method dependencies.** Method modules never include each other;
   shared logic belongs in `grid_analysis.c` or `smooth.c`.

### LAPACK choices

| Method      | Routine  | Why |
|-------------|----------|-----|
| polyfit     | `dgelss` | SVD; tolerates rank-deficient Vandermonde near boundaries |
| savgol      | `dposv`  | Coefficient system is symmetric positive definite |
| tikhonov    | `dpbsv`  | $(D^2)^T W D^2 + I$ is pentadiagonal SPD (kd=2); $O(n)$ solve |
| butterworth | (none)   | Biquad cascade with analytical IC via Cramer's rule |

### Grid uniformity philosophy

`grid_analysis.c` computes the coefficient of variation $CV = \sigma(h)/h_{avg}$.
Each method uses CV to make a policy decision:

| CV          | Behaviour |
|-------------|-----------|
| $\le 0.01$  | `is_uniform = 1` |
| $> 0.05$    | Savgol **rejects** with detailed error (uniformity is a mathematical requirement, not implementation choice) |
| $< 0.15$    | Tikhonov uses average-coefficient discretization (uniform stencil $[1,-4,6,-4,1]/h^4$) |
| $\ge 0.15$  | Tikhonov switches to local-spacing discretization (weighted Gram matrix $\sum w_k \mathbf{d}_k^T \mathbf{d}_k$) |
| $> 0.15$    | Butterworth **rejects** (frequency analysis assumes uniform sampling) |

Polyfit tolerates any grid (local fit per window). When changing CV thresholds
or adding a new method, update **all** policy points consistently.

### Per-method design notes

These are the load-bearing design choices, not user-facing math (see README for
that):

- **Polyfit:** SVD per window with `rcond = 1e-10` to truncate small singular
  values. Asymmetric windows + polynomial extrapolation at boundaries.
  $O(n \cdot p^3)$.
- **Savgol:** Universal convolution coefficients pre-computed once via moment
  conditions. Translation invariance is the whole point — same coefficients
  applied at every interior point. Uniform-grid requirement enforced by CV
  check, not silently degraded.
- **Tikhonov:** True 2nd-order penalty $(D^2)^T W D^2$ (pentadiagonal Gram
  matrix), corrected in v5.11. Hybrid discretization (auto-switch at CV=0.15).
  GCV trace uses 2D null space (constants and linear functions are unpenalized).
- **Butterworth:** 4th-order low-pass split into a biquad cascade for numerical
  stability. Filtfilt (forward-backward) gives zero phase. Per-biquad analytical
  IC via Cramer's rule avoids LAPACK and `complex.h`. Auto-cutoff via Morozov's
  discrepancy principle (v5.11.3).

## Testing

Uses the **Unity** framework (vendored in `tests/`).

- 106 tests total: grid_analysis (7), polyfit (21), savgol (16), tikhonov (26),
  butterworth (20), timestamp (16). Source of truth is `tests/test_main.c`.
- Zero leaks. `make test-valgrind` exits 1 on any definite/indirect leak or
  memory error — keep it that way.

### Adding a test

1. Add the test function in the relevant `tests/test_<module>.c`. Use
   `TEST_ASSERT_DOUBLE_WITHIN()` for floating-point comparisons.
2. Declare the function and `RUN_TEST(...)` it in `tests/test_main.c`.
3. If the test exercises a new module, append the source to `TEST_SRC` and the
   object to `TEST_MODULES` in the `Makefile`.
4. Run `make test` and `make test-valgrind` before committing.

Test best practices:

- One logical assertion per test.
- Cover edge cases: `n=1`, `n=2`, NULL inputs, empty arrays.
- Always free in the test body — the leak budget is fixed (see above).

## Common development tasks

### Adding a new smoothing method

1. Create `newmethod.c` / `newmethod.h` following the existing pattern.
2. Define `NewmethodResult` with at minimum `y_smooth`, `y_deriv`, `n`.
3. Implement `newmethod_smooth(...)` accepting `const GridAnalysis*` —
   policy decisions about grid uniformity belong here.
4. Implement `free_newmethod_result()`.
5. In `smooth.c`: add `METHOD_NEWMETHOD` constant, dispatch case, output path.
6. Update `Makefile` (`SRC`, `HEAD`).
7. Bump `revision.h` (version + date + history line).
8. Add `tests/test_newmethod.c` and wire it into `tests/test_main.c` and the
   `Makefile`.

### Modifying grid analysis

`GridAnalysis` is shared state — every change ripples to all methods.

1. Update the struct in `grid_analysis.h`.
2. Populate the new field in `analyze_grid()`.
3. Update every method that should react to it.
4. Update `tests/test_grid_analysis.c`.

### Memory management

```c
NewmethodResult* newmethod_smooth(...) {
    NewmethodResult *r = malloc(sizeof(*r));
    r->y_smooth = malloc(n * sizeof(double));
    r->y_deriv  = malloc(n * sizeof(double));
    return r;
}

void free_newmethod_result(NewmethodResult *r) {
    if (r) {
        free(r->y_smooth);
        free(r->y_deriv);
        free(r);
    }
}
```

Verify with `make test-valgrind` before committing.

### Diagnostic output convention (butterworth, applies to new modules)

- `stdout` as `# ...` — info that should be preserved with the saved data
  (selected parameters, grid CV, numerical-quality warnings).
- `stderr` as `Warning: ...` — runtime/operational concerns (memory usage)
  that do not belong in the data file.
- `stderr` as `ERROR: ...` — hard failures; function returns NULL/non-zero.

## Build notes for development

- Default compiler: `clang` (override with `CC=gcc`).
- LAPACK/BLAS required: `-llapack -lblas`.
- Default library path: `~/lib` (override with `LIBDIR=/path/to/libs`).
- Production: `-O2`. Debug: `make debug` → `-g -O0`.
- Standard: C99 or later.

## Hard rules (do not break)

- Do not introduce cross-dependencies between method modules.
- Do not change grid analysis behaviour without updating every method that
  consumes the affected field.
- Do not commit without `make test` passing and no new valgrind leaks.
- Do not silently degrade a method when its mathematical preconditions are
  violated — reject loudly (Savgol, Butterworth on non-uniform grids).
