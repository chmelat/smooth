# CLAUDE.md

How the `smooth` project is run: principles, rules, and where things live.
Implementation details (algorithms, solver choices, thresholds, complexity) are
in `README.md` (Appendix A/B) and must not be duplicated here.

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
├─ timestamp.c/h      # Timestamp parsing for `-T` mode
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

### Grid uniformity philosophy

Each method decides for itself what a non-uniform grid means for it: adapt,
tolerate, or reject. Where uniformity is a mathematical requirement of the
method, it is enforced, not an implementation choice to relax. When changing a
threshold or adding a method, update **all** policy points consistently
(and the README grid tables).

Grid diagnostics that no method consumes (dropout / sampling-regime detection)
are **advisory only**: they never change a method's behaviour or set
`reliability_warning`. Keep it that way unless you intend to change what every
normal run prints.

## Testing

Uses the **Unity** framework (vendored in `tests/`). The source of truth for
what runs is `tests/test_main.c`.

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
