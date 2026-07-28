# Plan 005: Detect missing samples in a nominally regular grid

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md` — unless a reviewer dispatched you and told you they
> maintain the index.
>
> **Drift check (run first)**: `git diff --stat d53266f..HEAD -- grid_analysis.c grid_analysis.h tests/test_grid_analysis.c smooth.c`
> If any of those changed since this plan was written, compare the "Current state"
> excerpts against the live code before proceeding; on a mismatch, treat it as a
> STOP condition.

## Status

- **Priority**: P2
- **Effort**: M
- **Risk**: LOW (purely additive fields; no existing field or output line changes)
- **Depends on**: none
- **Category**: feature (diagnostic)
- **Planned at**: commit `d53266f`, 2026-07-28

## Why this matters

Data logged at a fixed rate loses samples. When it does, the gap left behind is
an integer multiple of the base period: one lost sample gives $2h_0$, three
consecutive give $4h_0$. This is a *different* signature from general grid
non-uniformity, and nothing in the program reports it today.

Measured on the repo's own real data (`all.dat`, 14506 rows of 1 Hz logging):

```
$ ./smooth -T -g all.dat
#   CV = 0.078
#   Grid type: NON-UNIFORM
#   Detected clusters: 0
#   Recommendation: Grid is nearly uniform - standard methods work well
```

That data is missing **58 samples across 48 separate gaps** (max run 4, coverage
99.6%). The tool says "standard methods work well" and reports zero clusters.
`pt.dat` is worse: **1106 missing samples, 20.5% of the record gone**, again
invisible.

Neither existing signal can catch this, and neither should be stretched to:

- **CV** measures overall spread. `all.dat` sits at 0.078 — below every warning
  threshold in the file.
- **The cluster detector** (`grid_analysis.c:90-97`) requires
  `h_prev < 0.1*h_avg && h_curr > 10*h_avg`, i.e. neighbouring spacings differing
  by more than 100x while straddling `h_avg`. A dropout gap is 2x, occasionally
  5x. It will never fire. (See the TK5 re-analysis in `plans/README.md` for the
  full characterization of that detector's limits.)

This is a diagnostic gap, not a correctness bug: the smoothing results are fine.
But a user smoothing `pt.dat` has a fifth of the record missing and no way to
learn that from the tool.

## Current state

### `GridAnalysis` — `grid_analysis.h:11-26`, complete

```c
typedef struct {
    double h_min;           /* Minimum spacing */
    double h_max;           /* Maximum spacing */
    double h_avg;           /* Average spacing */
    double h_std;           /* Standard deviation of spacing */
    double ratio_max_min;   /* Ratio h_max/h_min */
    double cv;              /* Coefficient of variation (std/avg) */
    int is_uniform;         /* 1 if CV <= 0.01 (practical uniformity), 0 otherwise */
    int reliability_warning; /* 1 if reliability concerns exist */
    char warning_msg[512];  /* Warning message buffer */

    /* Additional statistics */
    int n_points;           /* Number of points */
    int n_clusters;         /* Number of detected clusters */
    double uniformity_score; /* Score 0-1, where 1 is perfectly uniform */
} GridAnalysis;
```

The struct is `calloc`'d at `grid_analysis.c:43` (so new fields default to zero)
and `free_grid_analysis()` at `:245-248` is a plain `free(analysis)` — there are
no nested allocations to track. **Keep it that way**: the temp array this plan
introduces must be freed inside `analyze_grid()`, never stored in the struct.

### `analyze_grid()` main loop — `grid_analysis.c:69-100`

The loop is single-pass and deliberately stores nothing:

```c
    for (i = 0; i < n-1; i++) {
        double h_curr = x[i+1] - x[i];
        if (h_curr <= 0) { /* non-monotonic -> free + return NULL */ }
        if (h_curr < analysis->h_min) analysis->h_min = h_curr;
        if (h_curr > analysis->h_max) analysis->h_max = h_curr;
        double dev = h_curr - analysis->h_avg;
        sum_sq_diff += dev * dev;
        /* cluster detection using h_prev ... */
        h_prev = h_curr;
    }
```

A median needs the spacings materialized, so this plan adds an `O(n)` temp array.
That is a deliberate departure from the "single-pass, no allocation" property the
comment at `:59-61` advertises — see "Design decisions" below.

### Existing thresholds — `grid_analysis.c:14-21`

```c
#define THRESH_CV_UNIFORM       0.01
#define THRESH_CV_NOTICEABLE    0.2
#define THRESH_CV_HIGH          0.5
#define THRESH_CV_SEVERE        1.0
#define CLUSTER_RATIO_SMALL     0.1
#define CLUSTER_RATIO_LARGE     10.0
#define UNIFORMITY_DECAY_FACTOR 2.0
```

### Output path — `grid_analysis.c:235-240`

```c
    if (verbose >= 1) {
        printf("%s  Standard deviation: %.6e\n", prefix, analysis->h_std);
        printf("%s  Detected clusters: %d\n", prefix, analysis->n_clusters);
        printf("%s  Recommendation: %s\n", prefix, get_grid_recommendation(analysis));
    }
```

`smooth.c:242` calls `print_grid_analysis(grid_info, 1, "# ")` for `-g`, so a
`verbose >= 1` line does appear in `-g` output. `smooth.c:248-252` calls it with
`verbose = 0` but only inside `if (grid_info->reliability_warning)`.

## The algorithm

```
h[i] = x[i+1] - x[i]                          (i = 0 .. n-2)
h_base = median(h)                            <- MEDIAN, not mean
for each h[i]:
    r = h[i] / h_base
    k = round(r)
    if |r - k| <= DROPOUT_TOL:
        near++                                 (counts toward integer_fit)
        if k >= 2:
            n_gaps++
            n_missing += k - 1
            max_run = max(max_run, k - 1)
integer_fit  = near / (n-1)
has_dropouts = (integer_fit >= DROPOUT_FIT_MIN) && (n_missing > 0)
```

Two choices carry the whole design.

**Median, not mean.** The mean is contaminated by exactly the gaps being looked
for. On `pt.dat` the mean spacing is 1.258 s while the true base period is 1.0 s;
the median returns 1.0 s. Measured, the median holds the correct `h_base` up to
45% sample loss. (The existing cluster detector's use of `h_avg` is the same
mistake and is why it needs blocks of 10+ points to fire at all.)

**`integer_fit` as the gate.** It separates "regular grid with dropouts" from
"genuinely non-uniform data", where counting missing samples is meaningless:

| data | `integer_fit` |
|---|---|
| regular with dropouts (up to 40% loss) | 100% |
| random points (Poisson) | 55% |
| geometric grading r=1.02 | 49% |
| alternating 0.1/0.4 | 0% |

The margin between 100% and 55% is wide, so `DROPOUT_FIT_MIN = 0.90` is not a
finely tuned number.

### Measured limits — document these, do not try to fix them

- **Loss >= 50% breaks it.** The median jumps to $2h_0$. Measured: 45% loss gives
  an exact count, 50% loss reports 111 missing instead of 1035. Sharp edge, not a
  gradual degradation.
- **Clock jitter above ~±15%.** Exact through ±15%; at ±20% it under-reports
  (47 of 59). Irrelevant for millisecond timestamps, relevant for hand-logged data.
- **False positive on an unbalanced alternating mesh.** A mesh that is mostly
  $h$ with an occasional exact multiple is indistinguishable from dropouts by
  construction. Measured: 30×0.1 alternating with 10×0.4 yields `h_base = 0.1`,
  `integer_fit = 100%`, and a report of 150 missing samples when nothing is
  missing. A balanced 20:20 mesh is correctly rejected (median lands between the
  two values, `integer_fit = 0%`). `integer_fit` cannot catch this case — it is a
  genuine limit of the method. Say so in the header comment; do not add a
  heuristic to paper over it.

## Design decisions

1. **Additive only.** No existing field changes meaning, no existing output line
   changes text or position. `all.dat`'s current `-g` output must keep every line
   it has today, with the new line appended in the `verbose >= 1` block.
2. **Do NOT set `reliability_warning`.** `smooth.c:248` gates a whole extra output
   block on it. Setting it here would make every normal (non-`-g`) run on data
   with a single dropped sample print a warning block it does not print today —
   a behaviour change well beyond "add a diagnostic". Report at `verbose >= 1`
   only. If a warning policy is wanted later, that is a separate decision with its
   own plan.
3. **Do NOT touch `n_clusters`.** Different question, different semantics. Adding
   dropouts to that counter would corrupt a field that already has a documented
   (if strict) meaning.
4. **Graceful degradation on allocation failure.** If the temp array cannot be
   allocated, leave the dropout fields zeroed and return the analysis anyway.
   `analyze_grid()` returning NULL is reserved for invalid input; a diagnostic
   extra must never turn a working run into a failure.
5. **`qsort` from the standard library**, on a copy of the spacings. `O(n log n)`
   on top of the existing `O(n)`. For the 10000-point test this is microseconds;
   do not hand-roll quickselect for an asymptotics win nobody will measure.
   Mark it: `/* ponytail: qsort is O(n log n) on top of an O(n) analysis;
   quickselect would be O(n) if profiling ever shows this matters. */`

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build | `make` | exit 0, no compiler warnings |
| Tests | `make test` | `130 Tests 0 Failures 0 Ignored` after this plan |
| Leak check | `make test-valgrind` | exit 0 |
| Real-data check | `./smooth -T -g all.dat` | reports 58 missing samples |

## Scope

**In scope** (the only files you may modify):

- `grid_analysis.h` — six new struct fields + doc comment
- `grid_analysis.c` — thresholds, median helper, detection, one report line
- `tests/test_grid_analysis.c` — five new tests
- `tests/test_main.c` — five declarations, five `RUN_TEST()` calls
- `CLAUDE.md` — test counts, and a row in the grid-philosophy section
- `revision.h` — version bump to 5.11.50 plus a history entry

**Out of scope** (do NOT touch):

- `smooth.c` — no CLI flag, no new call site. The `-g` path already passes
  `verbose = 1` and picks the new line up for free.
- `reliability_warning`, `warning_msg`, `n_clusters`, `cv`, `is_uniform`,
  `uniformity_score` — every existing field keeps its current value on every
  input.
- The cluster detector at `grid_analysis.c:90-97`. It is separately wrong-ish and
  separately documented; fixing it is not this plan.
- Any smoothing method. None of them read the new fields, so none need changes.
- `README.md`.

## Git workflow

- Branch: `advisor/005-dropout-detector`
- One commit.
- Message style follows `git log`: `feat:` prefix with the version in parens:
  `feat: detect missing samples in nominally regular grids (v5.11.50)`
- Do NOT push or open a PR.

## Steps

### Step 1: Capture the baseline

```bash
make >/dev/null
./smooth -T -g all.dat > /tmp/plan005-before.txt 2>&1
cat /tmp/plan005-before.txt
make test 2>&1 | tail -3
```

**Verify**: output contains `CV = 0.078` and `Detected clusters: 0`; `make test`
reports `125 Tests 0 Failures 0 Ignored`. Keep this file — step 7 diffs against it.

### Step 2: Add the struct fields

In `grid_analysis.h`, append to `GridAnalysis` after `uniformity_score`:

```c
    /* Dropout detection: a nominally regular grid with missing samples leaves
     * gaps that are integer multiples of the base period. Distinct from
     * n_clusters (abrupt spacing change) and from cv (overall spread).
     * All zero when has_dropouts == 0. */
    double h_base;       /* Median spacing — robust base period estimate */
    double integer_fit;  /* Fraction of spacings within DROPOUT_TOL of k*h_base */
    int    n_gaps;       /* Gaps that are an integer multiple k >= 2 of h_base */
    int    n_missing;    /* Estimated missing samples, sum of (k-1) */
    int    max_run;      /* Longest run of consecutive missing samples */
    int    has_dropouts; /* 1 if integer_fit >= DROPOUT_FIT_MIN and n_missing > 0 */
```

**Verify**: `make 2>&1 | grep -iE "warning|error"` prints nothing.

### Step 3: Add thresholds and the median helper

In `grid_analysis.c`, after the existing `#define` block:

```c
#define DROPOUT_TOL        0.25  /* |r - round(r)| tolerance, in units of h_base */
#define DROPOUT_FIT_MIN    0.90  /* Min integer_fit to treat the grid as regular */
#define DROPOUT_MIN_SPACES 10    /* Below this the median is not meaningful */
```

and a comparison function for `qsort` next to `append_warning`:

```c
static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}
```

Note the comparator returns via two comparisons rather than `(int)(da - db)`,
which truncates to 0 for spacings closer than 1.0 apart — the common case here.

### Step 4: Compute the detection

In `analyze_grid()`, after the statistics are finalized (after the `cv` assignment
at `:114`) and **before** the warning block at `:128`:

```c
    /* Dropout detection — see grid_analysis.h. Skipped silently on short grids
     * or if the scratch allocation fails; the fields stay zero from calloc.
     * ponytail: qsort is O(n log n) on top of an O(n) analysis; quickselect
     * would be O(n) if profiling ever shows this matters. */
    if (n - 1 >= DROPOUT_MIN_SPACES) {
        double *hs = (double *)malloc((size_t)(n - 1) * sizeof(double));
        if (hs != NULL) {
            for (i = 0; i < n-1; i++) hs[i] = x[i+1] - x[i];
            qsort(hs, (size_t)(n - 1), sizeof(double), cmp_double);

            int m = n - 1;
            analysis->h_base = (m % 2) ? hs[m/2] : 0.5 * (hs[m/2 - 1] + hs[m/2]);

            if (analysis->h_base > 0.0) {
                int near = 0;
                for (i = 0; i < n-1; i++) {
                    double r = (x[i+1] - x[i]) / analysis->h_base;
                    double k = floor(r + 0.5);
                    if (fabs(r - k) <= DROPOUT_TOL) {
                        near++;
                        if (k >= 2.0) {
                            int missing = (int)k - 1;
                            analysis->n_gaps++;
                            analysis->n_missing += missing;
                            if (missing > analysis->max_run)
                                analysis->max_run = missing;
                        }
                    }
                }
                analysis->integer_fit = (double)near / (double)m;
                analysis->has_dropouts =
                    (analysis->integer_fit >= DROPOUT_FIT_MIN &&
                     analysis->n_missing > 0) ? 1 : 0;
            }
            free(hs);
        }
    }
```

Scan the *original* spacings for the ratio test, not the sorted copy — `hs` is
sorted in place and its order no longer matches the grid.

**Verify**:

```bash
make 2>&1 | grep -iE "warning|error"; echo "(nothing above = clean)"
```

### Step 5: Report it

In `print_grid_analysis()`, inside the existing `verbose >= 1` block, after the
`Detected clusters` line:

```c
        if (analysis->has_dropouts) {
            printf("%s  Missing samples: %d in %d gap(s), longest run %d "
                   "(base period %.6e, coverage %.1f%%)\n",
                   prefix, analysis->n_missing, analysis->n_gaps,
                   analysis->max_run, analysis->h_base,
                   100.0 * analysis->n_points /
                       (analysis->n_points + analysis->n_missing));
        }
```

Print nothing when `has_dropouts == 0` — a clean grid should not gain a line
saying so, and a genuinely non-uniform grid must not be given a meaningless
missing-sample count.

### Step 6: Confirm it works on the real data

```bash
./smooth -T -g all.dat | grep -i "missing"
```

**Verify**: reports `58` missing samples in `48` gaps, longest run `4`,
coverage `99.6%`. If the numbers differ, the algorithm is transcribed wrong —
STOP. (These were measured independently from the timestamps during planning.)

Also confirm the negative case:

```bash
python3 -c "
h=0.1; x=0.0
for k in range(200):
    print('%.10f %.6f'%(x, x*x)); x+=h; h*=1.02" > /tmp/plan005-graded.dat
./smooth -g /tmp/plan005-graded.dat | grep -ci missing
```

**Verify**: `0` — a geometrically graded mesh has `integer_fit ~= 0.49`, below the
gate, and must print no missing-sample line.

### Step 7: Confirm nothing else changed

```bash
./smooth -T -g all.dat > /tmp/plan005-after.txt 2>&1
diff /tmp/plan005-before.txt /tmp/plan005-after.txt
```

**Verify**: the *only* difference is the one added `Missing samples:` line. Any
change to `CV`, `Grid type`, `Uniformity score`, `Detected clusters` or
`Recommendation` means an existing field was disturbed — STOP.

### Step 8: Add the tests

Append to `tests/test_grid_analysis.c`, following the ARRANGE / ACT / ASSERT /
CLEANUP structure of the existing tests (model on `tests/test_grid_analysis.c:45-75`).

1. **`test_grid_dropouts_detected`** — uniform base 1.0, 40 points, drop indices
   10, 20, 21, 22 before building `x`. ASSERT `has_dropouts == 1`,
   `n_missing == 4`, `n_gaps == 2`, `max_run == 3`,
   `TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, h_base)`.
2. **`test_grid_dropouts_none_when_complete`** — perfectly uniform 40 points.
   ASSERT `has_dropouts == 0`, `n_missing == 0`, and `integer_fit == 1.0` (every
   spacing is exactly `1*h_base`).
3. **`test_grid_dropouts_rejects_graded_mesh`** — geometric grading r=1.05, 60
   points. ASSERT `has_dropouts == 0`. This is the discriminator test; without it
   the gate is untested.
4. **`test_grid_dropouts_median_survives_heavy_loss`** — base 1.0, 200 nominal
   points, deterministically keep every point except those where `i % 5 < 2`
   (40% loss, no `rand()`). ASSERT `h_base` within `1e-12` of 1.0 and
   `has_dropouts == 1`. Guards the median-vs-mean choice: the mean spacing here
   is ~1.67, so a mean-based implementation fails this test.
5. **`test_grid_dropouts_short_grid_guard`** — 6 points, uniform. ASSERT
   `has_dropouts == 0` and `h_base == 0.0` (below `DROPOUT_MIN_SPACES`, detection
   skipped), and that `analyze_grid` still returns non-NULL with `cv` computed.

All five: CLEANUP with `free_grid_analysis(result)`.

### Step 9: Register the tests

In `tests/test_main.c`: five declarations in the grid_analysis declaration block,
five `RUN_TEST()` calls in the grid_analysis runner block.

**Verify**:

```bash
grep -c "RUN_TEST(test_grid" tests/test_main.c
make test 2>&1 | tail -3
```

→ count rises by 5; suite reports `130 Tests 0 Failures 0 Ignored`.

### Step 10: Prove the tests have teeth

Temporarily change `analysis->h_base` to use the mean instead of the median
(`analysis->h_base = analysis->h_avg;`), rebuild, and confirm
`test_grid_dropouts_median_survives_heavy_loss` **FAILS**. Then restore.

```bash
make test 2>&1 | grep -E "median_survives|Tests .* Failures"
```

**Verify**: that test fails with the mean and passes with the median. If it passes
with the mean, the test does not test what it claims — STOP and report.

### Step 11: Leak check

```bash
make clean && make
make test
make test-valgrind
```

**Verify**: zero warnings, `130 Tests 0 Failures 0 Ignored`, valgrind exits 0.
Pay attention here — this plan adds the first `malloc`/`free` pair inside
`analyze_grid()` beyond the struct itself, and the early-`return NULL` path for
non-monotonic data at `:76-80` must not leak `hs`. (It cannot, if `hs` is
allocated after the main loop as specified. If you moved the allocation earlier,
you own that leak.)

### Step 12: Version bump

`revision.h`: `VERSION` to `5.11.50`, `REVDATE` to today, new entry at the top of
the history block (demote the current `(current)` marker). Carry: that dropouts
leave integer-multiple gaps, which neither CV nor the cluster detector can see;
the measured `all.dat` (58 missing across 48 gaps, previously reported as "nearly
uniform, clusters: 0") and `pt.dat` (1106 missing, 20.5%) figures; the
median-not-mean choice and why (`pt.dat` mean 1.258 vs true 1.0); `integer_fit` as
the gate with its measured margin (100% vs 55/49/0%); the three measured limits
(>=50% loss, jitter >±15%, unbalanced alternating mesh); and that no existing
field or output line changed and `reliability_warning` was deliberately not set.

**Verify**: `./smooth -h 2>&1 | grep -i version` reports `5.11.50`.

### Step 13: Update `CLAUDE.md`

- `125 tests total` → `130 tests total`, `grid_analysis (7)` → `grid_analysis (12)`.
- In the "Grid uniformity philosophy" section, add a line noting that dropout
  detection is independent of CV and advisory only (no method reads it).

**Verify**: `grep -c "130 tests total" CLAUDE.md` → `1`.

## Test plan

Five tests, described in step 8. The design rests on two properties that convert
vague questions into exact assertions:

- **Dropout counts are exact integers.** Building `x` by deleting known indices
  from a uniform grid means `n_missing`, `n_gaps` and `max_run` have single
  correct values, not tolerances.
- **The median/mean distinction is falsifiable.** At 40% loss the mean spacing is
  1.67x the true period, so test 4 separates the two implementations by a wide
  margin. That is what step 10 verifies.

Deliberately not tested: the false-positive case (unbalanced alternating mesh).
It is a documented limit of the method, and a test asserting the wrong-but-expected
answer would ossify it. Note it in the header comment instead.

## Done criteria

ALL must hold:

- [ ] `make clean && make` exits 0 with zero compiler warnings
- [ ] `./smooth -T -g all.dat` reports 58 missing in 48 gaps, longest run 4
- [ ] A geometrically graded mesh reports **no** missing-sample line
- [ ] Step 7 diff shows exactly one added line and no other change
- [ ] `make test` → `130 Tests 0 Failures 0 Ignored`
- [ ] Step 10 shows the median test FAILS with a mean-based `h_base`
- [ ] `make test-valgrind` exits 0, no definite or indirect leaks
- [ ] `grep -c 'define VERSION "5.11.50"' revision.h` → `1`, with a history entry
- [ ] `grep -c "130 tests total" CLAUDE.md` → `1`
- [ ] `git diff smooth.c` is empty
- [ ] `git status --porcelain` lists only `grid_analysis.c`, `grid_analysis.h`,
      `tests/test_grid_analysis.c`, `tests/test_main.c`, `CLAUDE.md`, `revision.h`

## STOP conditions

Stop and report back (do not improvise) if:

- The "Current state" excerpts do not match the live code.
- Step 6 gives a count other than 58 missing / 48 gaps / run 4 on `all.dat`.
- Step 7's diff shows any change beyond the single added line.
- The graded mesh in step 6 reports missing samples — the gate is not working.
- Step 10's median test passes against a mean-based implementation.
- Any of the existing 125 tests starts failing.
- `make test-valgrind` reports any leak — the new `malloc` is the first thing to
  look at.
- You find yourself editing `smooth.c`, `reliability_warning`, or `n_clusters`.

## Maintenance notes

- **What a reviewer should scrutinise**: step 7 (nothing else changed) and step 10
  (the median test has teeth). Everything else is routine.
- **Why `reliability_warning` is untouched**: `smooth.c:248` gates a whole output
  block on it, so setting it would change what every normal run prints for data
  with a single dropped sample. That is a policy decision deserving its own plan,
  not a side effect of adding a diagnostic.
- **`pt.dat` has CRLF line endings** and currently fails to parse at all
  (`ERROR: No valid data points found` — the `\r` makes every y value
  non-numeric). Unrelated to this plan and not fixed here, but it is why the
  verification steps use `all.dat`. Worth its own finding.
- **Do not merge this into the cluster detector.** They answer different
  questions: "is there an abrupt spacing change" vs "is this a regular grid with
  holes". A grid can have either, both, or neither.
- **Rejected alternative**: gating on `ratio_max_min` instead of a median-based
  integer test. It cannot distinguish one 5x gap in an otherwise perfect grid
  (a dropout) from a smoothly graded mesh with the same extremes (not a dropout),
  and gives no count of what is missing.
- **Rejected alternative**: autocorrelation or FFT to find the base period.
  Correct, and far more machinery than a median needs to be for a grid whose
  spacings are already integer multiples by construction.
