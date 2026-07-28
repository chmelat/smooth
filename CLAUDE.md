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
└─ parser.c/h         # Input table parsing, `#` comment stripping, column selection
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
| any         | Tikhonov uses one integral-measure scheme: weighted Gram matrix $\sum w_k \mathbf{d}_k^T \mathbf{d}_k$ (uniform stencil $[1,-4,6,-4,1]/h^3$ on uniform grids; no CV switch) |
| $> 0.15$    | Butterworth **rejects** (frequency analysis assumes uniform sampling) |

Polyfit tolerates any grid (local fit per window). When changing CV thresholds
or adding a new method, update **all** policy points consistently.

**Dropout detection is independent of CV** (v5.11.51). A regular grid with
missing samples leaves gaps of $k \cdot h_{base}$ for integer $k \ge 2$ — a
signature neither CV nor the cluster detector can see. It is keyed on the
**median** spacing, not `h_avg`, because the mean is contaminated by the gaps
being looked for. Purely advisory: no method reads `has_dropouts`, `n_missing`,
`h_base`, `n_gaps`, `max_run` or `integer_fit`, and `reliability_warning` is
deliberately not set from it. Keep it that way unless you intend to change what
every normal run prints.

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
  matrix), corrected in v5.11. Single integral-measure discretization (no CV
  switch; unified in v5.11.34). GCV trace uses 2D null space (constants and linear functions are unpenalized).
- **Butterworth:** 4th-order low-pass split into a biquad cascade for numerical
  stability. Filtfilt (forward-backward) gives zero phase. Per-biquad analytical
  IC via Cramer's rule avoids LAPACK and `complex.h`. Auto-cutoff via Morozov's
  discrepancy principle (v5.11.3). Padding length is fc-adaptive
  (`PAD_DECAY_FACTOR/(1-r_max)`, floored at 14, capped at n-1) so the transient
  is absorbed even for small fc (v5.11.41); a `# WARNING` is emitted when the cap
  engages.

## Testing

Uses the **Unity** framework (vendored in `tests/`).

- 138 tests total: grid_analysis (17), polyfit (21), savgol (16), tikhonov (27),
  butterworth (22), timestamp (18), parser (17). Source of truth is `tests/test_main.c`.
- Zero leaks. `make test-valgrind` exits 1 on any definite/indirect leak or
  memory error — keep it that way.

Step-by-step recipes — adding a test, adding a new smoothing method, modifying
grid analysis, the result-struct memory-management pattern, and the diagnostic
output convention — live in the `smooth-dev-tasks` skill
(`.claude/skills/smooth-dev-tasks/SKILL.md`).

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
