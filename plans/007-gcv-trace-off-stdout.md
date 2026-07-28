# Plan 007: Keep the GCV sweep trace out of the data stream, and make program output ASCII

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md` — unless a reviewer dispatched you and told you they
> maintain the index.
>
> **Drift check (run first)**: `git diff --stat fbf4130..HEAD -- tikhonov.c tikhonov.h smooth.c butterworth.c tests/test_main.c Makefile`
> If any of those changed since this plan was written, compare the "Current state"
> excerpts against the live code before proceeding; on a mismatch, treat it as a
> STOP condition.

## Status

- **Priority**: P2
- **Effort**: S–M
- **Risk**: LOW-MEDIUM (output format changes; no numerics change)
- **Depends on**: none
- **Category**: tech-debt
- **Planned at**: commit `fbf4130`, v5.11.52, 2026-07-28
- **Closes**: audit findings B8 (diagnostic volume) and B9 (UTF-8 in output)

## Why this matters

`help()` advertises the program as a Unix filter (`cat data | smooth`). On the
auto-lambda path it is not one. Measured at `fbf4130` on the tracked
100-point `test_data.dat`:

```
$ ./smooth -m 2 -l auto test_data.dat | grep -c '^#'
28
$ ./smooth -m 2 -l auto test_data.dat | grep -cP '[^\x00-\x7F]'
21
```

28 comment lines for 100 data rows, 21 of them non-ASCII. Every one of those 21
lines goes to **stdout**, i.e. straight into the file the user redirected the
result to. 21 of the 28 are a *progress trace*: one row per trial $\lambda$ in
the GCV sweep, a quantity nobody keeps alongside a smoothed dataset.

It gets worse when the grid is non-uniform. On a 60-point alternating
0.1/0.4 mesh:

```
$ ./smooth -m 2 -l auto nonuni.dat 2>&1 | grep -c '^#'
62
$ ./smooth -m 2 -l auto nonuni.dat 2>&1 | grep -c 'Trace(H) approximation'
21
```

62 diagnostic lines for 60 data rows — **more diagnostics than data** — and 21
of them are the *same sentence repeated*, because `tikhonov.c:334` sits inside
`compute_gcv_score_robust()`, which the sweep calls once per trial $\lambda$.

This plan does three separable things:

1. Route the per-$\lambda$ trace to **stderr**, keeping on stdout only what the
   existing convention calls "info that should be preserved with the saved
   data". Result: 28 stdout lines → 7 on `test_data.dat`.
2. Hoist the repeated `Trace(H)` note out of the per-$\lambda$ function so it
   prints once instead of 21 times.
3. Replace the 7 non-ASCII glyphs in program output strings with ASCII. Result:
   0 non-ASCII bytes on any stream.

Nothing about the numerics, the chosen $\lambda$, or the data rows changes.

### Why not a `-v` / `-q` flag

Because the repo already tried and deleted it. `revision.h:265-266` records, for
v5.11.43:

```
 *           tikhonov.c: compute_gcv_score_robust() loses its `verbose`
 *           parameter (both call sites passed 1).
```

A verbosity parameter existed, was never once set to anything but "on", and was
removed as dead flexibility in a ponytail audit. Re-adding it as a CLI flag
would re-litigate that decision and add CLI surface, plumbing, and a new state
variable to solve a problem that stream selection already solves for free.
**Do not add a flag.** If you find yourself wanting one, that is a STOP
condition.

### The rule that decides each line

It is already written down, in `.claude/skills/smooth-dev-tasks/SKILL.md:70-76`
and mirrored in `butterworth.c:7-13`:

> - `stdout` as `# ...` — info that should be preserved with the saved data
>   (selected parameters, grid CV, numerical-quality warnings).
> - `stderr` as `Warning: ...` — runtime/operational concerns (memory usage)
>   that do not belong in the data file.
> - `stderr` as `ERROR: ...` — hard failures.

Applying it: the *selected* $\lambda$ and the caveats about its reliability are
preserved-with-data. The 21-row search that found it is not. The convention
already answers the question; the auto-lambda path just never followed it.

## Current state

### The seven non-ASCII glyphs in output strings

Complete list — `grep -nP "[^\x00-\x7F]" *.c *.h | grep -vP ':\s*[*/]'` at
`fbf4130`, minus comment lines:

| Site | Glyph | Text |
|---|---|---|
| `tikhonov.c:345` | λ | `"# λ=%9.3e: J=..."` |
| `tikhonov.c:397` | λ | `"# Refinement around λ=%.6e\n"` |
| `tikhonov.c:417` | λ | `"# WARNING: optimal λ = %.3e lies at the edge..."` |
| `tikhonov.c:419` | λ | `"...consider setting λ manually (-l <value>)."` |
| `tikhonov.c:422` | λ | `"# Optimal λ: %.6e (GCV=%.3e)\n"` |
| `smooth.c:330` | × | `"# Actual cutoff frequency: f_cutoff = %.6lG (= fc × f_Nyquist)\n"` |
| `butterworth.c:512` | — (em dash) | `"Signal may be broadband — consider setting fc manually.\n"` |

The audit entry recorded only the λ occurrences. The `×` in `smooth.c` and the
em dash in `butterworth.c` are the same defect in two other files and are
**in scope** — closing the class, the way plan 003 folded `savgol.h` in.

The rest of the program is already ASCII, including `smooth.c:262` and
`smooth.c:274`, which spell the same parameter `lambda`.

### The routing table

Every diagnostic the auto-lambda path emits, and where it goes after this plan.
`->` means the stream changes.

| Site | Text | Now | After | Why |
|---|---|---|---|---|
| `tikhonov.c:334` | `# Note: Trace(H) approximation less accurate...(ratio=%.2f)` | stdout, **x21** | stdout, **x1**, hoisted | Quality caveat = preserve. Repetition is a defect. |
| `tikhonov.c:345` | `# lambda=...: J=..., RSS=..., tr(H)=..., GCV=...` | stdout | -> **stderr** | Search trace. |
| `tikhonov.c:372` | `# GCV optimization for n=%d points...` | stdout | -> **stderr** | Header of the search trace. |
| `tikhonov.c:373` | `# Grid CV = %.3f (integral-measure...)` | stdout | -> **stderr** | Also on stdout via `-g` and the grid-warning block. |
| `tikhonov.c:376` | `# WARNING: Highly non-uniform grid detected...` | stdout | stdout | Quality caveat on the returned lambda. |
| `tikhonov.c:397` | `# Refinement around lambda=...` | stdout | -> **stderr** | Search trace. |
| `tikhonov.c:417,419` | `# WARNING: optimal lambda ... edge of the search range` | stdout | stdout | Quality caveat on the returned lambda. |
| `tikhonov.c:422` | `# Optimal lambda: %.6e (GCV=%.3e)` | stdout | -> **stderr** | Conclusion of the trace. The value is **already** on stdout from `smooth.c:262`. |
| `smooth.c:262` | `# Automatic lambda selection using GCV: lambda = %.6e` | stdout | stdout | The caller's statement of what was applied. |

Note the duplication this resolves: at `fbf4130` the chosen $\lambda$ is printed
to stdout **three times** — `tikhonov.c:422`, `smooth.c:262`, `smooth.c:274`.
Moving 422 to stderr removes one without deleting or reformatting anything a
user might parse.

Lines that keep the `#` prefix **even on stderr**. That costs nothing and keeps
`2>&1` safe, which matters: `tests/test_parser.c:36` runs the binary with
`2>&1` and treats every `#`-prefixed line as a comment to skip.

### `compute_gcv_score_robust()` — `tikhonov.c:294-351`

The relevant fragment, verbatim:

```c
    /* Grid uniformity for trace-approximation validity (already in grid_info) */
    double ratio = grid_info->ratio_max_min;
    double h_avg = grid_info->h_avg;
    ...
    if (ratio > 2.0) {
        printf("# Note: Trace(H) approximation less accurate for non-uniform grid (ratio=%.2f)\n", ratio);
    }
```

It is `static` and has exactly two callers, `tikhonov.c:386` and `tikhonov.c:403`,
both inside `find_optimal_lambda_gcv()` (verified: `grep -rn
compute_gcv_score_robust *.c *.h tests/*.c`). Hoisting the note into
`find_optimal_lambda_gcv()` is therefore safe and leaves the scoring function
free of output entirely — the state v5.11.43 was reaching for when it removed
the `verbose` parameter.

`ratio` is used only by that `if`. Once the note moves, `ratio` becomes an
orphan local and must go with it; `h_avg` stays, it feeds the eigenvalue sum.

### The header claim that stops being true

`tikhonov.h:77`:

```
 *   - Prints detailed optimization information to stdout
```

After this plan the detailed information goes to stderr and only the summary
warnings stay on stdout. Plan 003 was written because header prose had drifted
from the code; do not create the next instance of it.

## Conventions

- `tikhonov.c` documents *why* in comments — see the block at `tikhonov.c:318-323`.
  Match that density: one short comment at the trace sites saying the stream
  choice is deliberate and pointing at the convention, not a comment per line.
- Tests use Unity with ARRANGE / ACT / ASSERT / CLEANUP comments; `tests/test_parser.c`
  is the model for driving the binary via `popen()`.
- `make test` must stay green, `make test-valgrind` at zero leaks.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build | `make` | exit 0, no compiler warnings |
| Tests | `make test` | `140 Tests 0 Failures 0 Ignored` after this plan |
| Leak check | `make test-valgrind` | exit 0 |
| Clean rebuild | `make clean && make` | exit 0, no warnings |

## Scope

**In scope** (the only files you may modify):

- `tikhonov.c` — the nine `printf` sites listed in the routing table, and the
  hoist of the ratio note
- `tikhonov.h` — the one stale line at `:77`
- `smooth.c` — the `×` glyph at `:330`, nothing else
- `butterworth.c` — the em dash at `:512`, nothing else
- `tests/test_output.c` — new file, two tests
- `tests/test_main.c` — two declarations, two `RUN_TEST()` calls
- `Makefile` — append the new test file to `TEST_SRC`
- `.claude/skills/smooth-dev-tasks/SKILL.md` — one bullet in the convention section
- `CLAUDE.md` — test counts
- `revision.h` — version bump to 5.11.53 plus a history entry

**Out of scope** (do NOT touch):

- **Any numerics.** The sweep range, the 13 points, the `n <= 5000` refinement
  rule, the trace formula, the GCV expression. This plan moves text between file
  descriptors; if a data row changes, you broke something.
- **Any CLI flag.** No `-v`, no `-q`, no `--quiet`. See "Why not a `-v` / `-q`
  flag" above.
- **`smooth.c:262`, `:274`, `:275`, `:280`.** Do not delete or reword the
  caller's lambda reporting. The duplication is resolved by moving
  `tikhonov.c:422`, not by deleting a line whose format a user may parse.
- **`stderr` messages that already exist** (`ERROR:`, `Warning:` in every
  module). They are already on the right stream.
- **`butterworth.c`'s auto-cutoff prints** other than the em dash. They are 3
  lines, not 21, and they name the selected parameter — correctly on stdout.
- **`README.md`.** Verified at `fbf4130`: it contains no sample output block
  showing any of the moved lines (`grep -n "GCV optimization\|Optimal λ\|Refinement around" README.md` → no match).
- **`grid_analysis.c`** and the `-g` report.

## Git workflow

- Branch: `advisor/007-gcv-trace-off-stdout`
- One commit.
- Message style follows `git log` (`fix:` prefix with the version in parens):
  `fix: keep the GCV sweep trace out of the data stream (v5.11.53)`
- Do NOT push or open a PR.

## Steps

### Step 1: Capture the baseline

```bash
mkdir -p /tmp/plan007
make >/dev/null

python3 - <<'EOF' > /tmp/plan007/nonuni.dat
import math
x = 0.0
for i in range(60):
    print("%.6f %.6f" % (x, math.sin(x) + 0.01 * ((i * 37) % 7 - 3)))
    x += 0.1 if i % 2 else 0.4
EOF

# Auto-lambda, streams kept apart
./smooth -m 2 -l auto test_data.dat >/tmp/plan007/before_uni.out 2>/tmp/plan007/before_uni.err
./smooth -m 2 -l auto /tmp/plan007/nonuni.dat >/tmp/plan007/before_nu.out 2>/tmp/plan007/before_nu.err

# Data rows only, for the "numerics did not move" check
for f in test_data.dat /tmp/plan007/nonuni.dat; do
  b=$(basename $f .dat)
  ./smooth -m 2 -l auto  "$f" 2>/dev/null | grep -v '^#' > /tmp/plan007/rows_auto_$b.txt
  ./smooth -m 2 -l 0.1 -d "$f" 2>/dev/null > /tmp/plan007/full_manual_$b.txt
  ./smooth -m 0 -d      "$f" 2>/dev/null > /tmp/plan007/full_poly_$b.txt
  ./smooth -m 3         "$f" 2>/dev/null > /tmp/plan007/full_butter_$b.txt
  ./smooth -g           "$f" 2>/dev/null > /tmp/plan007/full_grid_$b.txt
done

printf 'uni stdout #: %s  non-ascii: %s\n' \
  "$(grep -c '^#' /tmp/plan007/before_uni.out)" \
  "$(grep -cP '[^\x00-\x7F]' /tmp/plan007/before_uni.out)"
printf 'nu  stdout #: %s  Trace(H) repeats: %s\n' \
  "$(grep -c '^#' /tmp/plan007/before_nu.out)" \
  "$(grep -c 'Trace(H) approximation' /tmp/plan007/before_nu.out)"
make test 2>&1 | tail -3
```

**Verify**, and STOP on any mismatch:

- `uni stdout #: 28  non-ascii: 21`
- `nu  stdout #: 62  Trace(H) repeats: 21`
- both `.err` files are **empty** (`wc -c` → 0) — nothing goes to stderr today
- `make test` → `138 Tests 0 Failures 0 Ignored`

`full_*` files include the `#` header lines on purpose: for `-m 0`, `-m 3`,
`-l 0.1` and `-g` this plan must not change **one byte**, headers included.
Only `-l auto` is allowed to change, and there only on stdout comment lines.

### Step 2: De-glyph the seven output strings

Mechanical, no logic change. Replace, in this order:

- `tikhonov.c:345,397,417,419,422` — every `λ` in a string literal becomes
  `lambda`. Watch the two format-column widths: line 345's `"# λ=%9.3e: ..."`
  becomes `"# lambda=%9.3e: ..."`, which is 5 characters wider. That is fine —
  the columns after it stay aligned with each other because every row uses the
  same format string.
- `smooth.c:330` — `fc × f_Nyquist` becomes `fc * f_Nyquist`.
- `butterworth.c:512` — `broadband — consider` becomes `broadband; consider`.

Do not touch the λ, ², ᵀ, · and — characters inside **comments**. They are not
output, `README.md`'s font restrictions do not apply to source comments, and
touching them would bloat the diff (`butterworth.c` alone has ~25).

**Verify**:

```bash
make 2>&1 | grep -i warn; echo "build rc=$?"
grep -nP "[^\x00-\x7F]" *.c *.h | grep -vP ':\s*[*/]'
```

Expected: no compiler warnings, and the second command prints **nothing** —
every remaining non-ASCII byte in the sources is inside a comment.

### Step 3: Hoist the repeated `Trace(H)` note

In `compute_gcv_score_robust()` (`tikhonov.c`), delete:

```c
    if (ratio > 2.0) {
        printf("# Note: Trace(H) approximation less accurate for non-uniform grid (ratio=%.2f)\n", ratio);
    }
```

and the now-unused local:

```c
    double ratio = grid_info->ratio_max_min;
```

Keep `double h_avg = grid_info->h_avg;` — it feeds the eigenvalue sum.
Adjust the comment above them so it still describes what remains.

In `find_optimal_lambda_gcv()`, put the note next to the existing `cv > 0.2`
warning, so both non-uniformity caveats print once, on stdout, adjacent:

```c
    if (grid_info->cv > 0.2) {
        printf("# WARNING: Highly non-uniform grid detected. Trace approximation less accurate.\n");
    }
    if (grid_info->ratio_max_min > 2.0) {
        printf("# Note: Trace(H) approximation less accurate for non-uniform grid (ratio=%.2f)\n",
               grid_info->ratio_max_min);
    }
```

Both stay on **stdout**: they qualify the reliability of the lambda that is
saved with the data. The two predicates differ (`cv` vs `ratio_max_min`) and
both are kept — narrowing them to one is a separate decision, out of scope.

**Verify**:

```bash
make >/dev/null && ./smooth -m 2 -l auto /tmp/plan007/nonuni.dat 2>/dev/null \
  | grep -c 'Trace(H) approximation'
grep -n 'printf\|fprintf' tikhonov.c | sed -n '/29[0-9]:/,/35[0-9]:/p'
```

Expected: `1` (was 21), and `compute_gcv_score_robust()` contains **no**
`printf` or `fprintf` at all.

### Step 4: Route the trace to stderr

Convert exactly these five sites from `printf(...)` to
`fprintf(stderr, ...)`, keeping the `#` prefix and the text otherwise identical:

- the per-lambda row (was `tikhonov.c:345`)
- `# GCV optimization for n=%d points ...` (was `:372`)
- `# Grid CV = %.3f ...` (was `:373`)
- `# Refinement around lambda=%.6e` (was `:397`)
- `# Optimal lambda: %.6e (GCV=%.3e)` (was `:422`)

Leave `:376` and `:417,419` as `printf`.

Add one comment above the trace header explaining the split, in the module's
style — something like:

```c
    /* The per-lambda search trace goes to stderr: it is progress output, not
     * something to preserve alongside the smoothed data. The chosen lambda and
     * the caveats about its reliability stay on stdout (see the diagnostic
     * output convention in the smooth-dev-tasks skill). The '#' prefix is kept
     * so a caller merging the streams with 2>&1 still sees valid comments. */
```

**Verify**:

```bash
make >/dev/null
./smooth -m 2 -l auto test_data.dat        >/tmp/plan007/after_uni.out 2>/tmp/plan007/after_uni.err
./smooth -m 2 -l auto /tmp/plan007/nonuni.dat >/tmp/plan007/after_nu.out 2>/tmp/plan007/after_nu.err
printf 'uni  stdout #: %s (was 28)   stderr #: %s\n' \
  "$(grep -c '^#' /tmp/plan007/after_uni.out)" "$(grep -c '^#' /tmp/plan007/after_uni.err)"
printf 'nu   stdout #: %s (was 62)   stderr #: %s\n' \
  "$(grep -c '^#' /tmp/plan007/after_nu.out)" "$(grep -c '^#' /tmp/plan007/after_nu.err)"
grep -cP '[^\x00-\x7F]' /tmp/plan007/after_uni.out /tmp/plan007/after_uni.err
```

Expected:

- `uni  stdout #: 7   stderr #: 21`
- `nu   stdout #: 17  stderr #: 25`
- both non-ASCII counts `0`

The 7 surviving stdout lines on `test_data.dat` must be exactly:

```
# WARNING: optimal lambda = 1.000e-08 lies at the edge of the search range [1e-08, 1e+00].
#          The true optimum may lie outside this range; consider setting lambda manually (-l <value>).
# Automatic lambda selection using GCV: lambda = 1.000000e-08
# Data smooth - Tikhonov regularization with lambda = 1E-08
# Functional J = 3.55564E-08 (Data: 1.06876E-14 + Regularization: 3.55564E-08)
# Data/Total ratio = 0.000, Regularization/Total ratio = 1.000
#    x          y
```

If `nu` stdout is 16 rather than 17, the ratio note from step 3 is not firing —
go back to step 3. If it is 18+, you moved too little.

### Step 5: Prove no data row moved

```bash
for f in test_data.dat /tmp/plan007/nonuni.dat; do
  b=$(basename $f .dat)
  ./smooth -m 2 -l auto  "$f" 2>/dev/null | grep -v '^#' | diff -q - /tmp/plan007/rows_auto_$b.txt
  ./smooth -m 2 -l 0.1 -d "$f" 2>/dev/null | diff -q - /tmp/plan007/full_manual_$b.txt
  ./smooth -m 0 -d      "$f" 2>/dev/null | diff -q - /tmp/plan007/full_poly_$b.txt
  ./smooth -m 3         "$f" 2>/dev/null | diff -q - /tmp/plan007/full_butter_$b.txt
  ./smooth -g           "$f" 2>/dev/null | diff -q - /tmp/plan007/full_grid_$b.txt
done; echo "all diffs rc=$?"
```

**Verify**: no `differ` line, `rc=0`. Byte-identity for `-m 0`, `-m 3`,
`-l 0.1` and `-g` including their `#` headers; for `-l auto`, identity of the
data rows.

Caveat worth knowing: `-m 3` output *does* contain the changed `×` line — but
only when the auto-cutoff path prints it. If `full_butter_*` differs on exactly
that one line, that is step 2 working as intended, not a regression. Confirm
with `diff` and move on. Any difference in a **numeric** field is a STOP.

### Step 6: Add the regression tests

New file `tests/test_output.c`. It guards the program's *stream contract*, which
is why it does not live in `tests/test_tikhonov.c` (unit tests, no binary) or in
`tests/test_parser.c` (whose helper merges the streams with `2>&1`, so it is
structurally unable to tell the two apart).

Two tests, both driving `./smooth` via `popen()`:

1. `test_output_gcv_trace_not_on_stdout` — run
   `./smooth -m 2 -l auto <fixture> 2>/dev/null`, read stdout, assert:
   - no line contains `"lambda="` followed by a `GCV=` field (the sweep row
     signature); match on the substring `"GCV="` in a line that also starts
     `"# lambda="`
   - no line contains `"GCV optimization for n="`
   - no line contains `"Refinement around"`
   - at least one line contains `"Automatic lambda selection using GCV"` —
     the positive half, so a change that silences *everything* also fails
   - `"Trace(H) approximation"` appears **at most once**

2. `test_output_is_ascii_on_both_streams` — for each of
   `-m 2 -l auto`, `-m 3`, `-m 0 -d`, run with `2>&1` and assert every byte read
   is `< 0x80`. Merging the streams is correct *here*: the claim is that
   nothing the program writes anywhere is non-ASCII.

Use a fixture generated into `/tmp` by the test itself, following
`write_fixture()` in `tests/test_parser.c:27-33`. Do not depend on
`test_data.dat` — it is untracked (this is the exact defect that made plan 001's
first dispatch stop; see `plans/README.md:36-40`). A ~40-point noisy sine
written by the test is enough; assert `n >= 20` data rows came back so a fixture
that fails to parse cannot pass the test vacuously.

Wire it up:

- `Makefile:19` — append `$(TEST_DIR)/test_output.c` to `TEST_SRC`.
  Nothing else: `TEST_MODULES` is for new *source* modules, not test files.
- `tests/test_main.c` — two declarations in a new
  `// Testy pro output stream contract (smooth.c via popen)` block after the
  parser block, and two `RUN_TEST()` calls in the matching place.

**Verify**: `make test` → `140 Tests 0 Failures 0 Ignored`.

### Step 7: Prove the tests have teeth

Both new tests must fail against the pre-change code. In a disposable worktree,
so nothing in your branch is disturbed:

```bash
git worktree add /tmp/plan007/old fbf4130
cp tests/test_output.c /tmp/plan007/old/tests/
cd /tmp/plan007/old
# apply the same three wiring edits (Makefile + test_main.c) by hand
make test 2>&1 | tail -20
```

**Verify**: `test_output_gcv_trace_not_on_stdout` **FAILS** (the sweep is on
stdout at `fbf4130`) and `test_output_is_ascii_on_both_streams` **FAILS** (λ is
on stdout at `fbf4130`). If either passes, the assertion is toothless — rewrite
it before continuing.

Then `cd` back and `git worktree remove /tmp/plan007/old --force`.

### Step 8: Make the header and the skill true again

- `tikhonov.h:77` — replace
  `- Prints detailed optimization information to stdout`
  with two lines that say what is actually true: the per-lambda search trace
  goes to stderr as `# ...`, the selected lambda and its reliability warnings go
  to stdout as `# ...`.
- `.claude/skills/smooth-dev-tasks/SKILL.md:70-76` — add a fourth bullet:

```
- `stderr` as `# ...` — iterative/progress traces (the Tikhonov GCV sweep) that
  would otherwise flood a redirected data file. Keep the `#` prefix so callers
  merging with `2>&1` still see valid comments.
```

- `butterworth.c:7-13` mirrors the same convention block. Add the matching line
  there so the two do not drift apart — that drift is what plan 003 existed to
  fix.

**Verify**: `git diff -- '*.h' '*.md'` touches only these three places, and
`git diff -- butterworth.c` shows only the comment block plus the step-2 em dash.

### Step 9: Leak check and full sweep

```bash
make clean && make 2>&1 | grep -iE 'warn|error'; make test 2>&1 | tail -3
make test-valgrind; echo "valgrind rc=$?"
```

**Verify**: no warnings, `140 Tests 0 Failures 0 Ignored`, `valgrind rc=0`.

### Step 10: Version bump and history entry

- `revision.h` — `#define VERSION "5.11.53"` and a history entry at the top.
  Mark it as a behaviour change, following the v5.11.44 precedent
  (`revision.h:274` uses `BEHAVIOUR CHANGE:` in caps). Say: the GCV sweep trace
  moved from stdout to stderr, so `smooth -m 2 -l auto data > out.dat` no longer
  captures it — `2>>out.dat` recovers the old behaviour; the repeated `Trace(H)`
  note now prints once; all program output is ASCII.
- `CLAUDE.md` — `138 tests total` → `140 tests total`, and append `output (2)`
  to the per-module list.

**Verify**:

```bash
grep -c 'define VERSION "5.11.53"' revision.h    # 1
grep -c "140 tests total" CLAUDE.md              # 1
grep -c "output (2)" CLAUDE.md                   # 1
git status --porcelain
```

The last must list exactly: `tikhonov.c`, `tikhonov.h`, `smooth.c`,
`butterworth.c`, `tests/test_output.c`, `tests/test_main.c`, `Makefile`,
`.claude/skills/smooth-dev-tasks/SKILL.md`, `CLAUDE.md`, `revision.h`
(plus the pre-existing untracked scratch files, which you must not add).

## Test plan

Two tests, described in step 6, with step 7 as the proof they discriminate.

The design point: **the positive assertion is what makes the negative ones
safe.** "stdout contains no sweep rows" is trivially satisfiable by a change
that breaks auto-lambda entirely, so the first test also requires
`"Automatic lambda selection using GCV"` to be present and `n >= 20` data rows
to come back. Without those, the test would bless a program that prints nothing
at all.

Deliberately not tested: the exact *content* of stderr. Asserting that the sweep
rows land on stderr, in order, with a particular format, would pin down the one
thing this plan makes explicitly cheap to change — the trace is now
developer-facing output. What must not regress is stdout.

Also not tested: interleaving order between the two streams. See maintenance
notes.

## Done criteria

ALL must hold:

- [ ] `make clean && make` exits 0 with zero compiler warnings
- [ ] `grep -nP "[^\x00-\x7F]" *.c *.h | grep -vP ':\s*[*/]'` prints nothing
- [ ] `./smooth -m 2 -l auto test_data.dat 2>/dev/null | grep -c '^#'` → `7` (was 28)
- [ ] `./smooth -m 2 -l auto test_data.dat 2>&1 | grep -cP '[^\x00-\x7F]'` → `0` (was 21)
- [ ] the nonuni fixture reports `Trace(H) approximation` exactly once (was 21)
- [ ] step 5 shows byte-identity for `-m 0`, `-m 3`, `-l 0.1`, `-g` (modulo the
      single `×` → `*` line) and identical data rows for `-l auto`
- [ ] `make test` → `140 Tests 0 Failures 0 Ignored`
- [ ] step 7 shows both new tests FAIL against `fbf4130`
- [ ] `make test-valgrind` exits 0
- [ ] `grep -c 'define VERSION "5.11.53"' revision.h` → `1`, with a matching
      history entry marked `BEHAVIOUR CHANGE:`
- [ ] `tikhonov.h` no longer claims the detailed information goes to stdout
- [ ] no new CLI option exists: `git diff smooth.c` shows no change to the
      `getopt` string or the `switch (ch)` block
- [ ] `git status --porcelain` lists only the ten files from step 10

## STOP conditions

Stop and report back (do not improvise) if:

- Step 1's baseline does not reproduce (28/21, 62/21, empty stderr, 138 tests).
  The plan was measured at `fbf4130`; a mismatch means the code moved.
- Any **numeric** field in any output differs before vs after. This plan
  redirects text; it must not perturb a single digit.
- You find yourself wanting a `-v`, `-q`, or `verbose` parameter. That path was
  deleted in v5.11.43 on purpose — report instead.
- You find yourself editing the sweep range, the point count, the refinement
  threshold, or the GCV formula.
- Either new test passes against `fbf4130` in step 7.
- `tests/test_parser.c` starts failing. Its helper merges `2>&1`; if a moved
  line breaks it, the `#` prefix was dropped somewhere — restore it rather than
  changing the parser tests.
- `make test-valgrind` reports any definite or indirect leak.

## Maintenance notes

- **What a reviewer should scrutinise**: step 7 first — without it the two new
  tests prove nothing. Then step 5's byte-identity, because the one genuine risk
  in this plan is an accidental format-string edit while retyping a `printf` as
  an `fprintf`. Diff the string literals character by character; the change to
  each of the five lines should be exactly `printf(` → `fprintf(stderr, ` plus
  `λ` → `lambda`.
- **Stream interleaving is now possible and is not a bug.** stdout is
  block-buffered when redirected to a file or pipe while stderr is unbuffered,
  so on a terminal the trace may appear before stdout lines that were "printed"
  earlier. Nothing in the program depends on relative order, and the `#` prefix
  keeps a merged stream valid. Do **not** add `fflush(stdout)` calls to force an
  order — that is real complexity bought for a cosmetic property nobody asked
  for. If someone eventually does ask, the cheap fix is one `setvbuf` call in
  `main()`, not scattered flushes.
- **Documented behaviour change.** `smooth -m 2 -l auto data > out.dat` no
  longer captures the GCV table. Interactively nothing is lost — stderr still
  goes to the terminal. `2> trace.log` or `2>> out.dat` recovers it. This is why
  the revision entry carries the `BEHAVIOUR CHANGE:` marker.
- **The two non-uniformity caveats are still near-duplicates.** `cv > 0.2` and
  `ratio_max_min > 2.0` produce two lines that say almost the same thing, now
  adjacent on stdout. Merging them into one predicate is a defensible follow-up,
  but it is a *policy* decision about grid thresholds and `CLAUDE.md` requires
  those to be changed consistently across every method — out of scope for a
  stream-routing plan.
- **Rejected alternative — delete the trace entirely.** Tempting (it is 21 lines
  nobody reads) but it is the only visibility into why a given lambda was
  chosen, and `tikhonov.c` has a real history of lambda-selection bugs
  (v5.11.39, v5.11.44). Redirecting costs one word per line and keeps the
  diagnostic; deleting saves nothing measurable and loses it.
- **Rejected alternative — gate the trace behind the existing `-g` flag.** `-g`
  means "grid analysis report"; overloading it with "and also print the GCV
  sweep" makes one flag mean two things and would still put the trace on stdout,
  which is the actual complaint.
- **Third instance of the same class.** Plan 003 fixed prose in two headers that
  had drifted from the code; step 8 fixes a third line in a third place. Nothing
  structural prevents a fourth, because the compiler never reads comments. If it
  happens again, that is the signal to add "grep the headers and the convention
  blocks when you change a public signature or an output stream" to the
  `smooth-dev-tasks` skill as a standing checklist item.
