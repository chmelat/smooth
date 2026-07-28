# Plan 004: Make the Tikhonov output derivative second-order on non-uniform grids

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md` — unless a reviewer dispatched you and told you they
> maintain the index.
>
> **Drift check (run first)**: `git diff --stat c8d3384..HEAD -- tikhonov.c tests/test_tikhonov.c`
> If either file changed since this plan was written, compare the "Current state"
> excerpts against the live code before proceeding; on a mismatch, treat it as a
> STOP condition.

## Status

- **Priority**: P2
- **Effort**: S
- **Risk**: LOW
- **Depends on**: none
- **Category**: bug
- **Planned at**: commit `c8d3384`, 2026-07-27

## Why this matters

`compute_derivatives()` in `tikhonov.c` produces the `-d` output using a
difference that is symmetric in *index* but not in *coordinate*:

```c
y_deriv[i] = (y_smooth[i+1] - y_smooth[i-1]) / (x[i+1] - x[i-1]);
```

With $h_l = x_i - x_{i-1}$ and $h_r = x_{i+1} - x_i$, a Taylor expansion gives

$$\frac{u_{i+1}-u_{i-1}}{x_{i+1}-x_{i-1}} = u' + \frac{h_r-h_l}{2}u'' + O(h^2)u'''$$

The error term $\frac{h_r-h_l}{2}u''$ vanishes only when $h_l = h_r$. Measured on
an alternating 0.1/0.4 mesh with $u=\sin x$ and $\lambda=0$ (so `y_smooth == y`
exactly and the derivative error is isolated): predicted max error 1.499e-01,
measured 1.494e-01 — the theory matches to three significant figures.

This matters because Tikhonov is the method the program **recommends** when the
grid is not uniform. `savgol.c:233-234` tells the user:

```
  1. Use Tikhonov method: -m 2 -l auto
     (Works correctly with non-uniform grids)
```

Savgol rejects above CV 0.05 and Butterworth above CV 0.15, so Tikhonov is where
non-uniform data ends up. The whole v5.11 line of work made its *penalty* correct
for such grids — the weighted Gram matrix $(D^2)^TWD^2$ with an integral measure.
The derivative it emits never got that care. Today a user who needs a derivative
on an abruptly non-uniform grid gets a better answer from `-m 0` (polyfit, which
differentiates its local fit analytically) than from the method the tool points
them to.

After this plan the interior and both boundaries use second-order formulas derived
by undetermined coefficients, at the same cost and the same $O(n)$.

## Current state

### The function to change — `tikhonov.c:87-111`, complete

```c
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
    
    /* Forward difference at start */
    y_deriv[0] = (y_smooth[1] - y_smooth[0]) / (x[1] - x[0]);
    
    /* Central differences in interior */
    for (int i = 1; i < n-1; i++) {
        y_deriv[i] = (y_smooth[i+1] - y_smooth[i-1]) / (x[i+1] - x[i-1]);
    }
    
    /* Backward difference at end */
    y_deriv[n-1] = (y_smooth[n-1] - y_smooth[n-2]) / (x[n-1] - x[n-2]);
}
```

It is called once, from `tikhonov_smooth()` at `tikhonov.c:248`. Nothing else in
the tree calls it (`static`). `x` is guaranteed strictly increasing by the
validation at `tikhonov.c:174-179`, so every $h$ below is strictly positive and no
division-by-zero guard is needed.

The `n < 2` and `n == 2` early returns stay exactly as they are — the new formulas
need three points, and those two branches are precisely the cases with fewer.

### The replacement formulas

All three come from undetermined coefficients on three points: choose weights that
reproduce $u$, $u'$ and $u''$ exactly, so the leading error is $O(h^2)u'''$. For
three points the closed form is short enough that no solver is needed at runtime.

**Interior**, $1 \le i \le n-2$, with $h_l = x_i - x_{i-1}$, $h_r = x_{i+1} - x_i$:

$$u'_i \approx -\frac{h_r}{h_l(h_l+h_r)}u_{i-1} + \frac{h_r-h_l}{h_l h_r}u_i + \frac{h_l}{h_r(h_l+h_r)}u_{i+1}$$

**Left boundary**, with $h_0 = x_1 - x_0$, $h_1 = x_2 - x_1$:

$$u'_0 \approx -\frac{2h_0+h_1}{h_0(h_0+h_1)}u_0 + \frac{h_0+h_1}{h_0h_1}u_1 - \frac{h_0}{h_1(h_0+h_1)}u_2$$

**Right boundary**, with $h_A = x_{n-1}-x_{n-2}$, $h_B = x_{n-2}-x_{n-3}$:

$$u'_{n-1} \approx \frac{2h_A+h_B}{h_A(h_A+h_B)}u_{n-1} - \frac{h_A+h_B}{h_Ah_B}u_{n-2} + \frac{h_A}{h_B(h_A+h_B)}u_{n-3}$$

Sanity check for the interior formula on a uniform grid ($h_l=h_r=h$): the weights
become $-\frac{1}{2h}$, $0$, $\frac{1}{2h}$ — the classical central difference.

### Measured effect

Alternating 0.1/0.4 mesh, $u = \sin x$, $\lambda = 0$:

| | current | new | factor |
|---|---|---|---|
| interior, max abs error | 1.49e-01 | 6.64e-03 | 22x |
| left boundary (on `exp`) | 5.17e-02 | 9.72e-03 | 5x |
| left boundary (on `x^3`) | 3.10e-01 | 5.00e-02 | 6x |
| left boundary, **uniform** mesh | 1.10e-01 | 7.76e-03 | 14x |

Note the last row: the current boundary formula is a two-point one-sided
difference, first-order **even on a uniform grid**. The boundary improvement is
therefore unconditional; it does not depend on non-uniformity at all.

### The floating-point consequence you must know about

The new interior formula is mathematically identical to the old one on a uniform
grid but **not bitwise identical** — it is three multiply-adds instead of one
subtract and one divide. Measured on a uniform mesh: 17 of 18 interior points
differ, by up to 1.6e-14 relative (~70 ULP). Two further facts:

- Even a "uniform" grid built by accumulating `x += h` does not have
  $h_l = h_r$ exactly: measured 15 of 18 interior points equal, 3 not. You cannot
  branch on `h_l == h_r` to preserve bit-identity, and you should not try.
- **The printed output is unaffected.** `print_result()` in `smooth.c` uses
  `%12.8lG` — 8 significant digits. Across 5600 uniform-grid interior points
  spanning spacings 0.01–3.7 and amplitudes 1e-3–1e3, the number of points whose
  *printed* value differs is **zero**.

So the repo's byte-identity release check still holds for interior points on
uniform grids. The two boundary points will change in the printed output, because
their formula genuinely improves. Expect that; it is the point of the change.

### Conventions

- The module's comment style documents *why*, with the discretization spelled out
  — see the block at `tikhonov.c:49-59` explaining the penalty stencil. Match that
  density: state that the formulas come from undetermined coefficients and that
  they degenerate to the classical stencils on a uniform grid.
- Tests use Unity with the ARRANGE / ACT / ASSERT / CLEANUP comment structure —
  see `tests/test_tikhonov.c:96-125` for a representative example.
- `make test` must stay green and `make test-valgrind` must stay at zero leaks.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build | `make` | exit 0, no compiler warnings |
| Tests | `make test` | `125 Tests 0 Failures 0 Ignored` after this plan |
| Leak check | `make test-valgrind` | exit 0 |
| Clean rebuild | `make clean && make` | exit 0, no warnings |

## Scope

**In scope** (the only files you may modify):

- `tikhonov.c` — `compute_derivatives()` only
- `tests/test_tikhonov.c` — two new tests
- `tests/test_main.c` — two declarations, two `RUN_TEST()` calls
- `CLAUDE.md` — test counts
- `revision.h` — version bump to 5.11.49 plus a history entry

**Out of scope** (do NOT touch):

- `build_band_matrix()` and `compute_functional()` — the penalty discretization is
  correct and was settled in v5.11.34. This plan changes only the reported
  derivative, not what is minimized.
- `tikhonov.h` — the public contract does not change. (It was just corrected in
  plan 003; do not re-edit it.)
- The `n < 2` and `n == 2` branches of `compute_derivatives()`.
- Every other method's derivative. `polyfit.c` differentiates its local fit
  analytically and is already grid-agnostic; `savgol.c` and `butterworth.c` assume
  uniform spacing **by design** and enforce it with a CV reject. Do not "fix" them.
- `README.md`.

## Git workflow

- Branch: `advisor/004-tikhonov-nonuniform-derivative`
- One commit.
- Message style follows `git log`: `fix:` prefix with the version in parens, e.g.
  `fix: prefix every line of multi-line grid warnings (v5.11.47)`. Use:
  `fix: second-order Tikhonov derivative on non-uniform grids (v5.11.49)`
- Do NOT push or open a PR.

## Steps

### Step 1: Capture the baseline

```bash
mkdir -p /tmp/plan004
python3 - <<'EOF' > /tmp/plan004/alt.dat
import math
x = 0.0
print("%.10f %.12f" % (x, math.sin(x)))
for k in range(40):
    x += 0.1 if k % 2 == 0 else 0.4
    print("%.10f %.12f" % (x, math.sin(x)))
EOF
make >/dev/null
./smooth -m2 -l 0 -d /tmp/plan004/alt.dat 2>/dev/null | grep -v '^#' > /tmp/plan004/before.txt
wc -l /tmp/plan004/before.txt
make test 2>&1 | tail -3
```

**Verify**: `before.txt` has 41 lines; `make test` reports
`123 Tests 0 Failures 0 Ignored`.

`-l 0` makes the system matrix the identity, so `y_smooth == y` exactly and the
`-d` column is purely the differentiation error. Keep using it for measurement.

### Step 2: Replace the three formulas in `compute_derivatives()`

Edit only the body after the `n == 2` early return. Keep the early returns.

Target shape — compute the spacings explicitly and apply the weights:

```c
    /* First derivative by undetermined coefficients on three points: the
     * weights are chosen to reproduce u, u' and u'' exactly, leaving an
     * O(h^2)u''' error on any spacing. The index-symmetric difference used
     * before was second-order only for h_l == h_r; on a non-uniform grid its
     * leading error is (h_r - h_l)/2 * u''. On a uniform grid the interior
     * weights below reduce to the classical -1/(2h), 0, +1/(2h).
     * x is strictly increasing (validated in tikhonov_smooth), so every
     * spacing is > 0. */

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
```

The boundary blocks index `x[2]` and `x[n-3]`, which requires `n >= 3`. That holds
here: `n < 2` and `n == 2` returned earlier.

**Verify**:

```bash
make 2>&1 | grep -iE "warning|error"; echo "(nothing above = clean)"
```

### Step 3: Confirm the accuracy actually improved

```bash
./smooth -m2 -l 0 -d /tmp/plan004/alt.dat 2>/dev/null | grep -v '^#' > /tmp/plan004/after.txt
python3 - <<'EOF'
import math
def err(path):
    m = 0.0
    for line in open(path):
        p = line.split()
        if len(p) < 3: continue
        m = max(m, abs(float(p[2]) - math.cos(float(p[0]))))
    return m
b, a = err('/tmp/plan004/before.txt'), err('/tmp/plan004/after.txt')
print("before max |error| = %.3e" % b)
print("after  max |error| = %.3e" % a)
print("improvement        = %.0fx" % (b/a))
EOF
```

**Verify**: `before` is ~1.5e-01, `after` is ~7e-03, improvement roughly 20x. If
`after` is not at least 5x better than `before`, a weight is wrong — STOP.

### Step 4: Add the discriminating unit test

Append to `tests/test_tikhonov.c`, following the ARRANGE / ACT / ASSERT / CLEANUP
structure of the existing tests (model on `tests/test_tikhonov.c:96-125`).

**Test 1 — `test_tikhonov_derivative_exact_for_quadratic_nonuniform`**

This is the test with teeth. The three-point weights reproduce $u''$ exactly, so
for a quadratic ($u'''=0$) the new formulas are exact to machine precision, while
the old ones are off by $\frac{h_r-h_l}{2}\cdot 2c$.

- ARRANGE: `x` alternating spacing 0.1 / 0.4, `n = 21` starting at 0.0.
  `y[i] = a + b*x[i] + c*x[i]*x[i]` with `a=1.0, b=-2.0, c=3.0`.
  `lambda = 0.0` so `y_smooth == y` exactly and only the derivative is under test.
- ACT: `analyze_grid(x, n)` then `tikhonov_smooth(x, y, n, 0.0, grid)`.
- ASSERT: for **every** `i` including `0` and `n-1`,
  `TEST_ASSERT_DOUBLE_WITHIN(1e-9, b + 2.0*c*x[i], result->y_deriv[i])`.
- CLEANUP: `free_tikhonov_result`, `free_grid_analysis`.

Measured margins, so you know the tolerance is honest rather than lucky:

| | old formula error | new formula error |
|---|---|---|
| interior, max | 9.00e-01 | 6.75e-14 |
| left boundary | 3.00e-01 | 1.78e-15 |
| right boundary | 1.20e+00 | 1.71e-13 |

A tolerance of `1e-9` passes the new code with four orders of margin and fails the
old code by nine. Assert over all points — dropping the boundaries would lose half
the value, and the existing linear test at `tests/test_tikhonov.c:116` already
skips them.

**Test 2 — `test_tikhonov_derivative_uniform_grid_still_exact`**

Guards against a regression in the common case while the weights are being
rearranged.

- ARRANGE: uniform `x[i] = i * 0.25`, `n = 21`, same quadratic coefficients,
  `lambda = 0.0`.
- ACT: same.
- ASSERT: same expression over all `i`, tolerance `1e-9`.
- CLEANUP: same.

Do **not** write a test that asserts the uniform-grid derivative is bit-identical
to the old values. It is not, by ~70 ULP, and such a test would be wrong.

### Step 5: Register the tests

In `tests/test_main.c`: add the two declarations to the tikhonov declaration block
and the two `RUN_TEST()` calls to the tikhonov runner block (the calls near
`tests/test_main.c:345-354`).

**Verify**:

```bash
grep -c "RUN_TEST(test_tikhonov" tests/test_main.c
make test 2>&1 | tail -3
```

→ the count rises by 2; the suite reports `125 Tests 0 Failures 0 Ignored`.

### Step 6: Prove the new tests would have caught the old bug

Temporarily restore the old body of `compute_derivatives()`, rebuild, and confirm
the new tests fail. Then restore your fix. This is the check that proves the tests
have teeth rather than merely passing.

```bash
cp tikhonov.c /tmp/plan004/tikhonov.c.fixed
git stash push -- tikhonov.c        # back to the old implementation
make test 2>&1 | grep -E "quadratic_nonuniform|Tests .* Failures"
git stash pop                        # restore the fix
make test 2>&1 | tail -3
```

**Verify**: with the old code, `test_tikhonov_derivative_exact_for_quadratic_nonuniform`
**FAILS**; after `git stash pop`, all 125 pass again. If the test passes against the
old implementation, it is not testing what it claims — STOP and report.

### Step 7: Leak check and full sweep

```bash
make clean && make
make test
make test-valgrind
```

**Verify**: zero compiler warnings, `125 Tests 0 Failures 0 Ignored`,
valgrind exits 0 with no definite/indirect leaks.

### Step 8: Version bump and history entry

In `revision.h`: `VERSION` to `5.11.49`, `REVDATE` to today, and a new entry at the
top of the history block (demote the current `(current)` marker). The entry should
carry: the wrong leading error term $\frac{h_r-h_l}{2}u''$ and that theory matched
measurement to three figures; that Tikhonov is the method the tool recommends for
non-uniform grids, so this was the worst place for the assumption to hide; the
three replacement formulas and where they come from; the measured 22x interior and
14x uniform-grid boundary improvements; that printed output is unchanged for
interior points on uniform grids (5600 points checked) while the two boundary
points do change; and that the penalty discretization was not touched.

**Verify**: `./smooth -h 2>&1 | grep -i version` reports `5.11.49`.

### Step 9: Update the test counts in `CLAUDE.md`

`123 tests total` → `125 tests total`, `tikhonov (25)` → `tikhonov (27)`.

**Verify**: `grep -c "125 tests total" CLAUDE.md` → `1`.

## Test plan

Two new tests, both in `tests/test_tikhonov.c`, described in step 4.

The design rests on one property worth stating plainly: **the three-point weights
are exact for quadratics.** That converts a fuzzy "is it more accurate?" question
into an exact assertion with a nine-order-of-magnitude margin between pass and
fail. It is why the test uses a quadratic rather than `sin`, and why `lambda = 0`
— any smoothing would mix penalty bias into the residual and force a loose
tolerance.

Deliberately not tested: smoothly graded meshes. Measured during planning, a
geometric grading of ratio 1.02–1.30 shows **no** improvement (factor ~1x) despite
CV up to 0.88, because there $h_r - h_l$ is comparable to the new formula's own
truncation error. A test on a graded mesh would have no teeth. The predictor of
damage is local asymmetry $|h_r-h_l|/h$, not CV — which is also why an alternating
mesh (CV 0.608, 22x) and a graded mesh (CV 0.622, 1x) behave completely
differently at nearly the same CV.

## Done criteria

ALL must hold:

- [ ] `make clean && make` exits 0 with zero compiler warnings
- [ ] Step 3 reports at least a 5x reduction in max derivative error
- [ ] `make test` → `125 Tests 0 Failures 0 Ignored`
- [ ] Step 6 shows the new quadratic test FAILS against the old implementation and passes after restoring the fix
- [ ] `make test-valgrind` exits 0
- [ ] `grep -c 'define VERSION "5.11.49"' revision.h` → `1`, with a matching history entry
- [ ] `grep -c "125 tests total" CLAUDE.md` → `1` and `grep -c "tikhonov (27)" CLAUDE.md` → `1`
- [ ] `git diff -- build_band_matrix compute_functional` shows no change to those functions (inspect `git diff tikhonov.c` by eye: only `compute_derivatives` differs)
- [ ] `git status --porcelain` lists only `tikhonov.c`, `tests/test_tikhonov.c`, `tests/test_main.c`, `CLAUDE.md`, `revision.h`

## STOP conditions

Stop and report back (do not improvise) if:

- `tikhonov.c:87-111` does not match the "Current state" excerpt.
- Step 3 shows less than a 5x improvement, or the error gets worse — a weight is
  transcribed wrong.
- Step 6's quadratic test passes against the old implementation.
- Any existing test among the 123 starts failing. In particular
  `test_tikhonov_linear_exact_null_space` (`tests/test_tikhonov.c:96-125`) must
  still pass: a linear function is exactly differentiated by both the old and new
  formulas, so a failure there means a transcription error, not a tolerance issue.
- You find yourself editing `build_band_matrix()` or `compute_functional()`.
- `make test-valgrind` reports any definite or indirect leak.

## Maintenance notes

- **What a reviewer should scrutinise**: step 6. Without it, two tests that pass
  prove nothing. Also check the assertions cover `i == 0` and `i == n-1` — the
  boundaries are where the largest improvement is, and the pre-existing linear
  test deliberately skips them.
- **Boundary improvement is unconditional.** The old two-point one-sided
  difference was first-order even on uniform grids, so `-d` output changes at the
  first and last point for *every* input, not just non-uniform ones. Anyone
  diffing release output should expect exactly two changed lines per run.
- **Do not add a `h_l == h_r` fast path** to preserve bit-identity on uniform
  grids. Floating-point accumulation means the spacings are not exactly equal even
  on a nominally uniform mesh (measured: 15 of 18), so such a branch would fire
  unpredictably and make the output depend on how `x` was constructed.
- **Still first-order, deliberately**: nothing here changes `compute_functional()`,
  which evaluates $\int (u'')^2$ with its own stencil. That is the penalty's
  discretization and is consistent with `build_band_matrix()`; the two must stay
  in step with each other, not with this derivative.
- **Rejected alternative**: a spline-consistent derivative,
  $u'_i = \frac{u_{i+1}-u_i}{h_r} - \frac{h_r(2M_i+M_{i+1})}{6}$, using the $M$
  that the penalty already computes and the natural BC $M_0 = 0$. Conceptually
  attractive — the Tikhonov solution is the discrete analogue of a cubic smoothing
  spline — but measured 2–4x *worse* at the boundary than the plain three-point
  formula (3.12e-02 vs 9.72e-03 on `exp`), at the cost of coupling the derivative
  to the penalty discretization. Not worth it.
- **Rejected alternative**: Fornberg's algorithm for general stencil weights. For
  three points a closed form exists; a generic solver would be pure overhead.
