# Plan 006: Replace the cluster detector with a sampling-regime detector

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md` — unless a reviewer dispatched you and told you they
> maintain the index.
>
> **Drift check (run first)**: `git diff --stat e314dfb..HEAD -- grid_analysis.c grid_analysis.h tests/test_grid_analysis.c README.md`
> If any of those changed since this plan was written, compare the "Current state"
> excerpts against the live code before proceeding; on a mismatch, treat it as a
> STOP condition.

## Status

- **Priority**: P1
- **Effort**: M
- **Risk**: LOW (diagnostic output only; no method reads any field involved)
- **Depends on**: 005 (this corrects it)
- **Category**: bug + feature
- **Planned at**: commit `e314dfb`, 2026-07-28

## Why this matters

Two problems, one fix.

### 1. The dropout detector shipped in v5.11.51 misreports common real data

Measured against the binary at `e314dfb`:

| Input | Reported today | Reality |
|---|---|---|
| Logger reconfigured 1 Hz to 10 Hz mid-record | `Missing samples: 2241 in 249 gap(s) ... coverage 18.2%` | **Nothing is missing.** The sampling rate changed. |
| Instrument restarted after a 1 h pause | `... coverage 12.2%` | Record is complete on both sides of a deliberate gap |
| Event-triggered burst sampling | `Missing samples: 9500 ... coverage 3.1%` | That is the intended sampling design |

The first is flatly wrong and is the most likely of the three in practice — it
happens whenever a logger is reconfigured or two runs are concatenated. The
detector cannot tell "samples lost from one regular record" from "the sampling
regime changed", because in both cases the spacings are integer multiples of the
median and `integer_fit` reaches 1.00.

`coverage` is the other half of the damage: it presents a deliberate pause as
88% data loss.

### 2. The cluster detector has never worked

`grid_analysis.c:102-109` requires `h_prev < 0.1*h_avg && h_curr > 10*h_avg`,
i.e. neighbouring spacings differing by more than **100x** while straddling the
global mean. Measured consequences: an isolated 200x gap in an otherwise uniform
grid reports `0`; dense blocks need 10+ points each to register (five blocks of
five score `0` even with a 1000x jump); only trailing edges count, so $k$ blocks
yield $k-1$. Across the entire scenario set used to plan this, it fired exactly
once (burst sampling, 19).

It has **zero test coverage** — `grep cluster tests/*.c` returns nothing.

Both problems have the same root: **a local phenomenon judged against a global
statistic.** The fix is the same lesson 005 already applied when it chose the
median over the mean.

## Current state

### The cluster detector — `grid_analysis.c:102-109`, complete

```c
        /* 3. Cluster Detection */
        /* Logic: A cluster boundary is a small gap followed by a large gap */
        if (i > 0) {
            if (h_prev < CLUSTER_RATIO_SMALL * analysis->h_avg && 
                h_curr > CLUSTER_RATIO_LARGE * analysis->h_avg) {
                analysis->n_clusters++;
            }
        }
```

`h_prev` exists in `analyze_grid()` solely to feed this; once the detector goes,
so does that variable and the `CLUSTER_RATIO_*` defines.

### Everything that touches `n_clusters`

```
grid_analysis.h:24      field declaration
grid_analysis.h:29      mentioned in the dropout-detection comment
grid_analysis.c:107     the increment above
grid_analysis.c:225,230 warning text "N abrupt spacing changes detected"
grid_analysis.c:301     report line "Detected clusters: %d"
README.md:1403          struct documentation in Appendix B
```

No method reads it. No test asserts on it.

### The dropout gate — `grid_analysis.c` around `:140` and `:306`

The block introduced by 005 computes `h_base` (median), `integer_fit`, `n_gaps`,
`n_missing`, `max_run` and sets `has_dropouts` when
`integer_fit >= DROPOUT_FIT_MIN && n_missing > 0`. The report at `:306` prints
`Missing samples: ... (base period ..., coverage ...%)`.

That gate is what this plan tightens; the median and the integer-multiple test
are correct and stay.

## The design

Three quantities, all computed without any global reference.

**Local jump ratio** — for each interior pair of spacings,

$$\rho_i = \frac{\max(h_i, h_{i+1})}{\min(h_i, h_{i+1})}$$

Record the maximum and the $x$ at which it occurs. Measured: a clean uniform grid
and a slowly drifting clock both give exactly 1.0; a graded mesh at ratio 1.10
gives 1.1; a rate change gives 10; a restart gives 3600.

**Mode fraction** — the share of spacings sitting at the base period itself,
$|h_i/h_{base} - 1| \le 0.25$. This separates one regime from several:

| Data | mode fraction |
|---|---|
| clean uniform | 100% |
| `all.dat` (real, 58 dropouts) | 99.7% |
| restart after 1 h pause | 99.8% |
| burst sampling | 93.6% |
| `pt.dat` (real, 20.5% loss) | 94.7% |
| **rate change 1 Hz to 10 Hz** | **50.1%** |

The bimodal case is isolated with a wide margin. `MODE_MIN = 0.85`.

**The classification.** Exactly one characterization per grid, in this order:

```
multi_regime = (mode_fraction < MODE_MIN) && (max_jump > JUMP_RATIO)
has_dropouts = integer_fit >= FIT_MIN && n_missing > 0 && !multi_regime

if      multi_regime  -> "Sampling: mixed regimes ..."
else if has_dropouts  -> "Missing samples: ..."
else if max_jump > JUMP_RATIO -> "Largest spacing jump ..."
else                  -> print nothing
```

**Both conditions are required for `multi_regime`.** A low mode fraction alone is
not enough: a smoothly drifting clock scores 34% and a graded mesh 3%, yet
neither has a jump (1.0 and 1.1), and calling them "mixed regimes" would be
wrong. This was caught by prototyping — the first version of this design got it
wrong. `JUMP_RATIO = 2.0`.

### Validation of the full decision logic

Every case below was run through a prototype of the exact logic above, and every
number the verification steps expect was taken from that prototype rather than
predicted:

| Input | Verdict |
|---|---|
| `all.dat` (real) | DROPOUTS: 58 in 48 gaps, longest 4, max jump 5.0x |
| `pt.dat` (real) | DROPOUTS: 1106 in 227 gaps, longest 13, max jump 14.0x |
| clean uniform | (nothing printed) |
| clock drift 0.3% | smoothly non-uniform, no jump |
| geometric grading r=1.10 | smoothly non-uniform, no jump |
| **rate change 1 Hz to 10 Hz** | **MIXED REGIMES, mode 50%, jump 10.0x at x = 249** |
| restart, 1 h pause | DROPOUTS: 3599 in 1 gap, max jump 3600x at x = 249 |
| burst sampling | DROPOUTS: 9500 in 19 gaps, max jump 501x at x = 36.05 |
| Poisson points | MIXED REGIMES, mode 18% |
| alternating 0.1/0.4 | MIXED REGIMES, mode 0%, jump 4.0x |

Restart and burst stay classified as dropouts. That is deliberate: the counts are
factually right (those samples were not taken), and with the jump location now
reported the user can see it is one interruption rather than scattered loss.
Only `coverage` was misleading, and it goes.

## Design decisions

1. **Delete `n_clusters` outright**, along with `CLUSTER_RATIO_SMALL`,
   `CLUSTER_RATIO_LARGE`, the `h_prev` variable, the warning at `:225-230` and
   the report line at `:301`. It has never worked, nothing reads it, and the user
   has explicitly waived backward compatibility for this diagnostic. Do not keep
   it "just in case".
2. **Drop `coverage` from the output.** It framed deliberate gaps as data loss.
   Replace its slot with the location of the largest jump, which is the thing a
   user actually needs — *where* in the record something happened. The current
   detector reports no location at all.
3. **Keep the median and the integer-multiple test from 005.** They are correct.
   Only the gate changes.
4. **Still advisory.** No method may read the new fields, and
   `reliability_warning` stays untouched — `smooth.c:248` gates a whole output
   block on it.
5. **No new allocation.** The jump scan and the mode count both fit in the loop
   that 005 already added over the spacings array.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build | `make` | exit 0, no compiler warnings |
| Tests | `make test` | `138 Tests 0 Failures 0 Ignored` after this plan |
| Leak check | `make test-valgrind` | exit 0 |
| Real-data check | `./smooth -T -g all.dat` | still reports 58 missing in 48 gaps |

## Scope

**In scope** (the only files you may modify):

- `grid_analysis.h` — remove `n_clusters`, add five fields
- `grid_analysis.c` — remove the cluster detector, add the jump/mode scan and the
  three-way report
- `tests/test_grid_analysis.c` — five new tests
- `tests/test_main.c` — five declarations, five `RUN_TEST()` calls
- `README.md` — struct in Appendix B, the `-g` example, the detection subsection
- `CLAUDE.md` — test counts, grid-philosophy note
- `revision.h` — version bump to 5.11.52 plus a history entry

**Out of scope** (do NOT touch):

- `smooth.c` — no CLI flag, no new call site.
- `reliability_warning`, `warning_msg`, `cv`, `is_uniform`, `uniformity_score`,
  `h_min`, `h_max`, `h_avg`, `h_std`, `ratio_max_min` — unchanged on every input.
- `h_base`, `integer_fit`, `n_gaps`, `n_missing`, `max_run` — the 005
  computation stays as it is; only `has_dropouts` gains a condition.
- Any smoothing method.

## Git workflow

- Branch: `advisor/006-sampling-regime-detector`
- One commit.
- Message: `feat: replace cluster detector with sampling-regime detection (v5.11.52)`
- Do NOT push or open a PR.

## Steps

### Step 1: Capture the baseline

```bash
make >/dev/null
./smooth -T -g all.dat > /tmp/plan006-before.txt 2>&1
python3 -c "
xs=[i*1.0 for i in range(250)]
t=xs[-1]
xs+=[t+(i+1)*0.1 for i in range(250)]
for v in xs: print('%.6f %.3f'%(v, v*0.1))" > /tmp/plan006-ratechange.dat
./smooth -g /tmp/plan006-ratechange.dat | grep -E "Missing|clusters"
make test 2>&1 | tail -3
```

**Verify**: the rate-change file reports `Missing samples: 2241 in 249 gap(s)`
and `Detected clusters: 0` — the bug this plan fixes. `make test` reports
`133 Tests 0 Failures 0 Ignored`.

### Step 2: Remove the cluster detector

Delete from `grid_analysis.c`: the `CLUSTER_RATIO_SMALL` and
`CLUSTER_RATIO_LARGE` defines, the detection block at `:102-109`, the `h_prev`
declaration and its `h_prev = h_curr;` assignment, the warning block at
`:225-230`, and the report line at `:301`. Delete the `n_clusters` field from
`grid_analysis.h` and fix the comment at `:29` that references it.

**Verify**:

```bash
grep -rn "n_clusters\|CLUSTER_RATIO\|h_prev" grid_analysis.c grid_analysis.h; echo "(nothing above = clean)"
make 2>&1 | grep -iE "warning|error"; echo "(nothing above = clean build)"
```

An unused-variable warning here means `h_prev` was left behind.

### Step 3: Add the new fields

In `grid_analysis.h`, replacing the deleted `n_clusters`:

```c
    /* Sampling regime. All local: no global mean is involved, because the
     * events these describe are local. max_jump is the largest ratio between
     * neighbouring spacings, mode_fraction the share of spacings sitting at the
     * base period itself. Advisory only; no method reads them. */
    double max_jump;      /* max(h_i,h_i+1)/min(h_i,h_i+1) over the grid; 1.0 if uniform */
    double max_jump_x;    /* x coordinate where max_jump occurs */
    int    n_jumps;       /* Neighbouring pairs whose ratio exceeds JUMP_RATIO */
    double mode_fraction; /* Share of spacings within DROPOUT_TOL of 1*h_base */
    int    multi_regime;  /* 1 if the record mixes two or more sampling regimes */
```

### Step 4: Compute them

Inside the block 005 added (the one guarded by `n - 1 >= DROPOUT_MIN_SPACES`,
after `h_base` is known), extend the existing scan and add the jump scan:

```c
                /* Mode fraction: share of spacings at the base period itself.
                 * One regime with holes keeps this near 1; two regimes split it. */
                int at_mode = 0;
                for (i = 0; i < n-1; i++) {
                    double r = (x[i+1] - x[i]) / analysis->h_base;
                    if (fabs(r - 1.0) <= DROPOUT_TOL) at_mode++;
                }
                analysis->mode_fraction = (double)at_mode / (double)m;

                /* Largest local jump between neighbouring spacings. */
                analysis->max_jump = 1.0;
                for (i = 0; i < n-2; i++) {
                    double a = x[i+1] - x[i], b = x[i+2] - x[i+1];
                    double rho = (a > b) ? a / b : b / a;
                    if (rho > JUMP_RATIO) analysis->n_jumps++;
                    if (rho > analysis->max_jump) {
                        analysis->max_jump = rho;
                        analysis->max_jump_x = x[i+1];
                    }
                }

                analysis->multi_regime =
                    (analysis->mode_fraction < MODE_MIN &&
                     analysis->max_jump > JUMP_RATIO) ? 1 : 0;
```

with the new thresholds next to the existing `DROPOUT_*` ones:

```c
#define MODE_MIN    0.85  /* Below this AND with a jump: several sampling regimes */
#define JUMP_RATIO  2.0   /* Neighbouring-spacing ratio that counts as a jump */
```

and `has_dropouts` gains one condition:

```c
                analysis->has_dropouts =
                    (analysis->integer_fit >= DROPOUT_FIT_MIN &&
                     analysis->n_missing > 0 &&
                     !analysis->multi_regime) ? 1 : 0;
```

**Both** parts of `multi_regime` are required. A smoothly drifting clock scores
mode 34% with max_jump 1.0; without the jump condition it would be misreported
as mixed regimes.

### Step 5: The three-way report

Replace the deleted `Detected clusters` line and rewrite the missing-sample line,
inside the existing `verbose >= 1` block:

```c
        if (analysis->multi_regime) {
            printf("%s  Sampling: mixed regimes - only %.0f%% of spacings at the "
                   "base period, largest jump %.3gx at x = %.6e\n",
                   prefix, 100.0 * analysis->mode_fraction,
                   analysis->max_jump, analysis->max_jump_x);
        } else if (analysis->has_dropouts) {
            printf("%s  Missing samples: %d in %d gap(s), longest run %d "
                   "(base period %.6e, largest jump %.3gx at x = %.6e)\n",
                   prefix, analysis->n_missing, analysis->n_gaps,
                   analysis->max_run, analysis->h_base,
                   analysis->max_jump, analysis->max_jump_x);
        } else if (analysis->max_jump > JUMP_RATIO) {
            printf("%s  Largest spacing jump: %.3gx at x = %.6e\n",
                   prefix, analysis->max_jump, analysis->max_jump_x);
        }
```

Note `coverage` is gone from the second line, replaced by the jump location.

**Verify**:

```bash
make 2>&1 | grep -iE "warning|error"; echo "(clean)"
```

### Step 6: Confirm the fix on the case that motivated it

```bash
./smooth -g /tmp/plan006-ratechange.dat | grep -E "Sampling|Missing"
```

**Verify**: reports `Sampling: mixed regimes - only 50% of spacings at the base
period, largest jump 10x at x = 2.490000e+02` and **no** `Missing samples` line.
If it still claims missing samples, the gate is not wired — STOP.

Then confirm the real data is unaffected:

```bash
./smooth -T -g all.dat | grep -i "missing"
```

**Verify**: still `58 in 48 gap(s), longest run 4`, now ending with
`largest jump 5x at x = ...` instead of `coverage 99.6%`.

### Step 7: Confirm the smooth cases stay silent

```bash
python3 -c "
h=0.1; x=0.0
for k in range(200): print('%.10f %.6f'%(x, x*x)); x+=h; h*=1.10" > /tmp/plan006-graded.dat
./smooth -g /tmp/plan006-graded.dat | grep -cE "Sampling|Missing|jump"
python3 -c "
for i in range(500): print('%.1f %.3f'%(i*1.0, i*0.01))" > /tmp/plan006-uniform.dat
./smooth -g /tmp/plan006-uniform.dat | grep -cE "Sampling|Missing|jump"
```

**Verify**: both print `0`. A graded mesh has a low mode fraction but no jump, and
a uniform grid has neither — neither may produce a line. If the graded mesh
reports mixed regimes, the `&&` in `multi_regime` was written as `||` — STOP.

### Step 8: Add the tests

Append to `tests/test_grid_analysis.c`, ARRANGE / ACT / ASSERT / CLEANUP as the
existing tests do.

1. **`test_grid_regime_rate_change_detected`** — 250 points at spacing 1.0 then
   250 at 0.1. ASSERT `multi_regime == 1`, `has_dropouts == 0`,
   `mode_fraction` within 0.05 of 0.50, `max_jump` within 0.1 of 10.0. This is
   the test with teeth: it fails against v5.11.51, which reports 2241 missing.
2. **`test_grid_regime_graded_mesh_is_not_mixed`** — geometric grading r=1.10,
   60 points. ASSERT `multi_regime == 0` and `max_jump < 2.0`, with
   `mode_fraction < 0.85`. Pins that a low mode fraction alone does not trigger
   the classification — the exact mistake the prototype made.
3. **`test_grid_regime_uniform_has_no_jump`** — 40 uniform points. ASSERT
   `max_jump` within 1e-12 of 1.0, `n_jumps == 0`, `multi_regime == 0`,
   `mode_fraction` within 1e-12 of 1.0.
4. **`test_grid_regime_isolated_gap_located`** — 20 points at spacing 1.0 from
   0.0, then a 200.0 gap, then 20 more at spacing 1.0. ASSERT `max_jump` within
   1.0 of 200.0 and `max_jump_x` within `1e-12` of **19.0**. This is the case the
   old cluster detector scored `0` on; assert the location, since reporting
   *where* is the point. Note this fixture also sets `has_dropouts == 1` (a 200x
   gap is an integer multiple, so 199 "missing" samples) — that is the restart
   case in miniature and is correct; do not assert against it.

   **Location semantics**: `max_jump_x` is `x[i+1]`, the point *between* the two
   spacings being compared — i.e. the last sample before a gap. With the fixture
   above the compared spacings are `h[18] = 1.0` and `h[19] = 200.0`, so the
   reported location is `x[19] = 19.0`. Measured, not assumed.
5. **`test_grid_regime_dropouts_still_reported`** — the fixture from
   `test_grid_dropouts_detected` (drop 10, 20, 21, 22 from 40 uniform points).
   ASSERT `has_dropouts == 1`, `multi_regime == 0`, `n_missing == 4` unchanged.
   Guards that tightening the gate did not break the 005 behaviour.

### Step 9: Register the tests

Five declarations and five `RUN_TEST()` calls in the grid_analysis blocks of
`tests/test_main.c`.

**Verify**: `make test` reports `138 Tests 0 Failures 0 Ignored`.

### Step 10: Prove the rate-change test has teeth

```bash
git stash push -- grid_analysis.c grid_analysis.h
make test 2>&1 | grep -E "rate_change|Tests .* Failures"
git stash pop
make test 2>&1 | tail -3
```

**Verify**: against v5.11.51 the suite fails to build or
`test_grid_regime_rate_change_detected` fails. (A build failure is acceptable
here — the new fields do not exist in the stashed header. If so, note it and
verify the behaviour manually instead: `./smooth -g /tmp/plan006-ratechange.dat`
on the stashed build must print the wrong `Missing samples: 2241` line.)

### Step 11: Leak check

```bash
make clean && make
make test
make test-valgrind
```

**Verify**: zero warnings, `138 Tests 0 Failures 0 Ignored`, valgrind exits 0.
No new allocation is introduced by this plan, so the alloc/free count should
match v5.11.51 exactly.

### Step 12: Version bump

`revision.h`: `VERSION` to `5.11.52`, `REVDATE` to today, new top entry. Carry:
that v5.11.51's dropout gate could not distinguish lost samples from a changed
sampling regime and reported 2241 phantom missing samples with 18.2% coverage on
a rate change; that `coverage` framed deliberate gaps as data loss and is gone;
that the cluster detector was removed rather than repaired, with the measured
reasons (needs >100x neighbouring ratio straddling `h_avg`, blocks of 10+ points,
counts trailing edges only, fired once across the whole scenario set, zero test
coverage); the new local quantities and the two-part `multi_regime` condition,
including that a low mode fraction alone would misclassify a drifting clock; and
the measured mode-fraction margin (50.1% against 93.6% and above).

### Step 13: Update `README.md` and `CLAUDE.md`

- `README.md` Appendix B struct: drop `n_clusters`, add the five new fields.
- `README.md` `-g` example: the `Detected clusters: 0` line goes; the
  `Missing samples` line ends with the jump location instead of coverage. Re-run
  the command and paste real output — do not hand-edit.
- `README.md` "Missing Sample Detection": document the regime gate and that
  `coverage` was replaced. Add the rate-change case to the limits.
- `CLAUDE.md`: `133 tests total` to `138`, `grid_analysis (12)` to `(17)`; update
  the grid-philosophy note to mention regime detection alongside dropouts.

**Verify**: `grep -c "138 tests total" CLAUDE.md` returns `1`;
`grep -c "n_clusters" README.md` returns `0`.

## Test plan

Five tests, step 8. Two carry the design:

- **Test 1** is the corrective. It fails against v5.11.51 by construction.
- **Test 2** pins the two-part `multi_regime` condition. Without it the `&&`
  could silently become `||` and every graded mesh would be misreported — the
  exact error the prototype made before measurement caught it.

Deliberately not tested: burst sampling and instrument restart. Both stay
classified as dropouts, which is a judgement call rather than a fact, and a test
would ossify it. They are documented in the history entry instead.

## Done criteria

ALL must hold:

- [ ] `make clean && make` exits 0 with zero compiler warnings
- [ ] Rate-change file reports mixed regimes and **no** missing-sample line
- [ ] `all.dat` still reports 58 missing in 48 gaps, longest run 4
- [ ] Graded mesh and uniform grid produce no regime/missing/jump line at all
- [ ] `make test` -> `138 Tests 0 Failures 0 Ignored`
- [ ] Step 10 shows the rate-change test cannot pass against v5.11.51
- [ ] `make test-valgrind` exits 0, alloc/free count unchanged from v5.11.51
- [ ] `grep -rn "n_clusters\|CLUSTER_RATIO" grid_analysis.c grid_analysis.h README.md` is empty
- [ ] `git diff smooth.c` is empty
- [ ] `grep -c 'define VERSION "5.11.52"' revision.h` -> `1`

## STOP conditions

Stop and report back (do not improvise) if:

- The "Current state" excerpts do not match the live code.
- The rate-change file still reports missing samples after step 5.
- A graded mesh or a uniform grid produces any line in step 7.
- `all.dat`'s missing count changes from 58 / 48 / 4.
- Any of the existing 133 tests starts failing — in particular the five
  `test_grid_dropouts_*` tests from 005 must all still pass.
- You find yourself editing `smooth.c`, `reliability_warning`, or any method.
- Removing `h_prev` produces an unused-variable warning you cannot resolve
  without touching the statistics loop.

## Maintenance notes

- **What a reviewer should scrutinise**: step 7. The two-part `multi_regime`
  condition is the one piece of this design that a reasonable person would get
  wrong, and a graded mesh reporting "mixed regimes" is the symptom.
- **Why the cluster detector was deleted, not fixed.** It answered a question
  ("is there an abrupt spacing change") that `max_jump` now answers better, with
  no threshold tied to a global mean and with a location attached. Keeping both
  would mean two fields for one question.
- **Restart and burst remain "dropouts" on purpose.** The counts are factually
  right and the jump location now shows they are interruptions rather than
  scattered loss. If a future report wants to separate them, the discriminator is
  `n_gaps` small with `max_run` large — do not reach for a new threshold on
  `coverage`, which is exactly the metric this plan removes.
- **Rejected alternative**: keeping `coverage` but suppressing it when
  `multi_regime`. It would still read as 12.2% data loss on a deliberate pause,
  which is the actual complaint. The location of the largest jump is more
  informative in every case measured.
- **Rejected alternative**: clustering the spacings (k-means, kernel density) to
  count regimes properly. Correct, and far more machinery than a mode fraction
  needs to be for a diagnostic that only has to separate "one regime" from "more
  than one".
