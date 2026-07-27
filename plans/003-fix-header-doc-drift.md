# Plan 003: Make the header documentation in `tikhonov.h` and `savgol.h` true

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md` — unless a reviewer dispatched you and told you they
> maintain the index.
>
> **Drift check (run first)**: `git diff --stat ce6016e..HEAD -- tikhonov.h savgol.h`
> If either file changed since this plan was written, compare the "Current state"
> excerpts against the live code before proceeding; on a mismatch, treat it as a
> STOP condition.

## Status

- **Priority**: P2
- **Effort**: S
- **Risk**: LOW
- **Depends on**: none
- **Category**: docs
- **Planned at**: commit `ce6016e`, 2026-07-27

## Why this matters

`tikhonov.h` is the first thing anyone reads before calling the module, and its
`Example usage` block **does not compile**. Three independent reasons: two calls
are missing the `grid_info` argument, the variable `result` is declared twice in
one scope, and the stated compile line omits a source file the example needs.
Three further claims in the same header state things the code does not do.

Documentation that is merely absent costs a reader time. Documentation that is
confidently wrong costs them a debugging session — they write the call the header
shows, it fails to compile, and they have to go read `tikhonov.c` to find out what
the real signature is. At that point the header is worse than nothing.

This finding (TK6) has been recorded as open since v5.11.39 and re-verified in the
v5.11.41 audit, which recommended fixing it first precisely because it is trivial
and the header actively lies.

`savgol.h` has the same class of defect and is included here: fixing only the
header the ticket names would leave its sibling equally uncompilable, and the next
audit would rediscover it.

After this plan: every example in both headers compiles and runs, and every factual
claim matches the code.

## Current state

### Files

- `tikhonov.h` — public interface of the Tikhonov module. Comments only; no code
  changes anywhere in this plan.
- `savgol.h` — public interface of the Savitzky-Golay module.

### The seven defects, each verified against the code

**D1 — `tikhonov.h:36`**, the `grid_info` parameter description:

```c
 *   grid_info - Grid analysis results (used for discretization method selection)
```

There is no discretization method selection. It was removed in v5.11.34, which
unified the module on a single integral-measure scheme. `build_band_matrix()`
(`tikhonov.c:34-85`) has one code path and never branches on grid statistics.

**D2 — `tikhonov.h:75`**, the GCV search range:

```c
 *   - Search range: 1e-6 to 1e0
```

The real range is at `tikhonov.c:334-335`:

```c
    const double lambda_min = 1e-8;
    const double lambda_max = 1e0;
```

**D3 — `tikhonov.h:77`**, a warning that does not exist:

```c
 *   - Warns if regularization term dominates excessively
```

`grep -rn "dominat" *.c *.h` returns exactly one line — this claim itself. What
the code actually warns about is different and undocumented: when the chosen
lambda lands on the edge of the search range (`tikhonov.c:391-395`).

**D4/D5/D6 — `tikhonov.h:91-123`**, the example block. Current text:

```c
/* Example usage:
 * 
 * #include "tikhonov.h"
 * 
 * int main() {
 *     double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
 *     double y[] = {1.0, 2.1, 3.9, 6.2, 8.1};
 *     int n = 5;
 *     
 *     // Method 1: Manual lambda
 *     TikhonovResult *result = tikhonov_smooth(x, y, n, 0.1);
 *     
 *     // Method 2: Automatic lambda selection
 *     double optimal_lambda = find_optimal_lambda_gcv(x, y, n);
 *     TikhonovResult *result = tikhonov_smooth(x, y, n, optimal_lambda);
 *     
 *     if (result != NULL) {
 *         printf("# Functional J = %.3e\n", result->total_functional);
 *         
 *         for (int i = 0; i < n; i++) {
 *             printf("%8.4f %8.4f %8.4f\n", 
 *                    x[i], result->y_smooth[i], result->y_deriv[i]);
 *         }
 *         
 *         free_tikhonov_result(result);
 *     }
 *     
 *     return 0;
 * }
 * 
 * Compilation:
 *   gcc -o program program.c tikhonov.c -llapack -lblas -lm
 */
```

- **D4**: lines `:101`, `:104`, `:105` omit `grid_info`. Real signatures
  (`tikhonov.h:55-56` and `:79`):
  ```c
  TikhonovResult* tikhonov_smooth(const double *x, const double *y, int n, double lambda,
                                  const GridAnalysis *grid_info);
  double find_optimal_lambda_gcv(const double *x, const double *y, int n, const GridAnalysis *grid_info);
  ```
- **D5**: `TikhonovResult *result` is declared at `:101` and again at `:105`, in
  the same scope. A redefinition error independent of D4.
- **D6**: the compile line omits `grid_analysis.c`. Once the example calls
  `analyze_grid()`, linking without it fails on an undefined reference. The
  example also never creates or frees a `GridAnalysis`, and never includes
  `<stdio.h>` despite calling `printf`.

**D7 — `savgol.h:73` and `savgol.h:79`**, both identical:

```c
 *   GridAnalysis *grid = analyze_grid(x, 6, 0);
```

The real signature has two parameters (`grid_analysis.h:36`):

```c
GridAnalysis* analyze_grid(const double *x, int n);
```

The third argument was `store_spacings`, removed earlier — `revision.h:98` records
the removal. Note that `savgol_smooth(x, y, 6, 5, 2, grid)` on the adjacent lines
`:74` and `:80` is **correct**; only the `analyze_grid` calls are stale.

### Claims that are TRUE — do not "fix" these

Each was checked against the code. If you change any of them you are introducing a
defect, not removing one.

| Claim | Location | Verified against |
|---|---|---|
| "For small datasets (n < 3), returns conservative default" | `tikhonov.h:76` | `tikhonov.c:342-345` returns `0.01` |
| "Uses grid search with refinement" | `tikhonov.h:73` | `tikhonov.c:355` (13 points) + `:371-385` (8 refinement steps) |
| "Prints detailed optimization information to stdout" | `tikhonov.h:74` | `tikhonov.c:320-321` |
| "Uses natural boundary conditions (second derivative = 0 at ends)" | `tikhonov.h:50` | `tikhonov.c:56-59` |
| "Efficient pentadiagonal band matrix (bandwidth 2) for O(n) memory" | `tikhonov.h:51` | `tikhonov.c:208-211`, `kd = 2` |
| "Correct discretization for both uniform and non-uniform grids" | `tikhonov.h:52` | `tikhonov.c:60-84` |
| `n` "must be >= 1", `lambda` ">= 0" | `tikhonov.h:32-33` | `tikhonov.c:163` |
| `savgol_smooth(x, y, 6, 5, 2, grid)` | `savgol.h:74,80` | `savgol.h:83-84` |
| `tikhonov_smooth(x, y, n, lambda, grid_info)` and `find_optimal_lambda_gcv(x, y, n, grid_info)` in the cross-reference | `savgol.h:98,101` | correct — `savgol.h` documents Tikhonov's API better than `tikhonov.h` does |
| `polyfit_smooth(x, y, n, window_size, poly_degree)` | `savgol.h:104` | `polyfit.h` |
| CV bands 0.01 / 0.05 | `savgol.h:127-129` | `savgol.c:24`, `savgol.c:217`, `savgol.c:244` |

### Repo conventions

- Header comment style is `/* ... */` with ` * ` continuation. Match the
  surrounding blocks exactly; do not reflow unrelated lines.
- Doc-only commits do **not** bump the version in this repo. Precedent:
  `902d74e`, `a3be811`, `732e930`, `d6fb01d` — all `doc:` prefixed, no version
  change.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build | `make` | exit 0, no compiler warnings |
| Tests | `make test` | `123 Tests 0 Failures 0 Ignored` / `OK` |
| Clean rebuild | `make clean && make` | exit 0, no warnings |
| Example compile | see step 5 | exit 0, no warnings |

Build uses `clang` with `-Wall -Wextra -pedantic -O2`, links `-llapack -lblas -lm`.
A clean build on `ce6016e` produces zero warnings.

## Scope

**In scope** (the only files you may modify):

- `tikhonov.h`
- `savgol.h`

**Out of scope** (do NOT touch):

- Every `.c` file. This plan changes comments only. If you find yourself editing
  code to make a doc claim true, STOP — that is a separate finding.
- `revision.h` — doc-only commits carry no version bump here.
- `README.md` — its Tikhonov section is separate and was reconciled in `732e930`.
- `doc/audit-tikhonov.md` — retiring the TK6 entry there is tracked separately in
  `plans/README.md`.
- `grid_analysis.h`, `polyfit.h`, `butterworth.h`, `parser.h`, `timestamp.h` — all
  checked; none contains an example block or a stale signature reference.

## Git workflow

- Branch: `advisor/003-fix-header-doc-drift`
- One commit.
- Message style follows `git log`; doc commits use the `doc:` prefix with no
  version suffix: `doc: fix stale API claims and examples in tikhonov.h and savgol.h`
- Do NOT push or open a PR.

## Steps

### Step 1: Capture the byte-identity baseline

A comments-only change must not alter program behaviour at all. Capture this
first so step 7 can prove it.

```bash
mkdir -p /tmp/plan003
printf '0 1.0\n1 2.1\n2 3.9\n3 6.2\n4 8.1\n5 9.8\n6 12.2\n' > /tmp/plan003/f.dat
make >/dev/null
./smooth -m2 -l 0.1 /tmp/plan003/f.dat  > /tmp/plan003/before_tik.txt 2>&1
./smooth -m1 -n5 -p2 /tmp/plan003/f.dat > /tmp/plan003/before_sg.txt  2>&1
wc -l /tmp/plan003/before_*.txt
```

**Verify**: both files non-empty.

### Step 2: Fix the three false claims in `tikhonov.h`

- `:36` — replace `(used for discretization method selection)` with a description
  of the real role: grid analysis results, required (non-NULL); the module uses a
  single integral-measure discretization and does not select between schemes.
- `:75` — `1e-6` becomes `1e-8`. Consider naming `tikhonov.c` as the source of
  truth so the next reader can check it in one step.
- `:77` — delete the "regularization term dominates" line. Replace it with what
  the code does warn about: the optimal lambda landing on the edge of the search
  range (`tikhonov.c:391-395`).
- `:73` — keep "Uses grid search with refinement" and add that refinement runs
  only for `n <= 5000` (`tikhonov.c:371`).

**Verify**:

```bash
grep -c "1e-6" tikhonov.h        # → 0
grep -rc "dominat" tikhonov.h    # → 0
grep -c "discretization method selection" tikhonov.h  # → 0
```

### Step 3: Rewrite the `tikhonov.h` example block

Replace lines `91-123` with the following. This exact program was compiled and run
before this plan was written — it is known-good, not a sketch. Keep the ` * `
comment prefixes.

```c
/* Example usage:
 *
 * #include <stdio.h>
 * #include "tikhonov.h"
 * #include "grid_analysis.h"
 *
 * int main(void) {
 *     double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
 *     double y[] = {1.0, 2.1, 3.9, 6.2, 8.1};
 *     int n = 5;
 *
 *     // Grid analysis is required; every method takes it as input.
 *     GridAnalysis *grid = analyze_grid(x, n);
 *     if (grid == NULL) return 1;
 *
 *     // Manual lambda:   tikhonov_smooth(x, y, n, 0.1, grid)
 *     // Automatic lambda (GCV):
 *     double lambda = find_optimal_lambda_gcv(x, y, n, grid);
 *     TikhonovResult *result = tikhonov_smooth(x, y, n, lambda, grid);
 *
 *     if (result != NULL) {
 *         printf("# Functional J = %.3e\n", result->total_functional);
 *
 *         for (int i = 0; i < n; i++) {
 *             printf("%8.4f %8.4f %8.4f\n",
 *                    x[i], result->y_smooth[i], result->y_deriv[i]);
 *         }
 *
 *         free_tikhonov_result(result);
 *     }
 *
 *     free_grid_analysis(grid);
 *     return 0;
 * }
 *
 * Compilation:
 *   gcc -o program program.c tikhonov.c grid_analysis.c -llapack -lblas -lm
 */
```

Note what changed and why, so you do not "simplify" it back:

- `<stdio.h>` and `grid_analysis.h` added — `printf` and `analyze_grid` need them.
- The two methods are collapsed into one path. The old block declared `result`
  twice; showing manual lambda as a comment keeps both documented without the
  redefinition.
- `analyze_grid()` / `free_grid_analysis()` added — the API cannot be called
  without a `GridAnalysis`.
- `int main(void)` rather than `int main()` — the build uses `-Wpedantic`.
- `grid_analysis.c` added to the compile line.

**Verify**: `grep -c "grid_analysis.c" tikhonov.h` → `1`.

### Step 4: Fix the two stale calls in `savgol.h`

Lines `:73` and `:79`, identical, one token each:

```c
 *   GridAnalysis *grid = analyze_grid(x, 6, 0);
```

becomes

```c
 *   GridAnalysis *grid = analyze_grid(x, 6);
```

Do not touch lines `:74`, `:80`, or the cross-reference block at `:95-106` — all
verified correct.

**Verify**:

```bash
grep -c "analyze_grid(x, 6, 0)" savgol.h   # → 0
grep -c "analyze_grid(x, 6)" savgol.h      # → 2
```

### Step 5: Prove the examples compile

This is the criterion that matters. Extract each example into a real translation
unit and build it.

```bash
mkdir -p /tmp/plan003
# tikhonov.h example -> strip the leading " * " comment prefixes
sed -n '/^\/\* Example usage:/,/^ \*\//p' tikhonov.h \
  | sed -e 's|^/\* Example usage:||' -e 's|^ \*/||' -e 's|^ \* \?||' \
  | sed -n '/#include/,/^}/p' > /tmp/plan003/ex_tik.c

clang -Wall -Wextra -pedantic -I. -o /tmp/plan003/ex_tik \
      /tmp/plan003/ex_tik.c tikhonov.c grid_analysis.c -llapack -lblas -lm
/tmp/plan003/ex_tik
```

**Verify**: the compile emits **no output at all** (no warnings) and exits 0; the
run prints a `# Functional J = ...` line followed by five data rows and exits 0.

Then the savgol snippets. They are fragments, not a program, so wrap them:

```bash
cat > /tmp/plan003/ex_sg.c <<'EOF'
#include <stdio.h>
#include "savgol.h"
#include "grid_analysis.h"
int main(void) {
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double y[] = {1.0, 2.1, 2.9, 4.2, 5.1, 6.0};
    GridAnalysis *grid = analyze_grid(x, 6);
    SavgolResult *result = savgol_smooth(x, y, 6, 5, 2, grid);
    printf("%s\n", result ? "ok" : "NULL");
    free_savgol_result(result);
    free_grid_analysis(grid);
    return 0;
}
EOF
clang -Wall -Wextra -pedantic -I. -o /tmp/plan003/ex_sg \
      /tmp/plan003/ex_sg.c savgol.c grid_analysis.c -llapack -lblas -lm
/tmp/plan003/ex_sg
```

**Verify**: no compiler output, exit 0; the run prints `ok`.

If the `sed` extraction in the first block produces something that does not look
like the example (check with `cat /tmp/plan003/ex_tik.c`), do not fight the
regex — paste the example into the file by hand and compile that. The point is to
compile the header's text, not to build a clever extractor.

### Step 6: Confirm no code changed

```bash
git diff --stat
git diff -- '*.c'
```

**Verify**: the first lists only `tikhonov.h` and `savgol.h`. The second produces
**no output**. If any `.c` file appears, STOP.

### Step 7: Build, test, and prove behaviour is unchanged

```bash
make clean && make
make test 2>&1 | tail -3
./smooth -m2 -l 0.1 /tmp/plan003/f.dat  > /tmp/plan003/after_tik.txt 2>&1
./smooth -m1 -n5 -p2 /tmp/plan003/f.dat > /tmp/plan003/after_sg.txt  2>&1
diff /tmp/plan003/before_tik.txt /tmp/plan003/after_tik.txt
diff /tmp/plan003/before_sg.txt  /tmp/plan003/after_sg.txt
```

**Verify**: clean build with zero warnings; `123 Tests 0 Failures 0 Ignored`; both
`diff` commands produce no output. A comments-only change that alters output means
something went badly wrong — STOP.

## Test plan

No new unit tests. This plan changes only comments, and the repo has no
documentation-testing harness; adding one for two examples would cost more than it
saves.

The verification that carries the weight is step 5: the header text is compiled as
real code. That is a stronger check than any assertion about the comment's
contents, and it is what a reader would do anyway.

Deliberately not done: a permanent `make check-docs` target that extracts and
compiles header examples on every build. Two examples do not justify the
extraction machinery, and the `sed` pipeline above is fragile enough that it would
become its own maintenance burden. Revisit only if a third example appears.

## Done criteria

ALL must hold:

- [ ] `grep -c "1e-6" tikhonov.h` → `0`
- [ ] `grep -rc "dominat" tikhonov.h` → `0`
- [ ] `grep -c "discretization method selection" tikhonov.h` → `0`
- [ ] `grep -c "grid_analysis.c" tikhonov.h` → `1`
- [ ] `grep -c "analyze_grid(x, 6, 0)" savgol.h` → `0` and `grep -c "analyze_grid(x, 6)" savgol.h` → `2`
- [ ] The `tikhonov.h` example compiles with zero warnings and runs (step 5)
- [ ] The `savgol.h` snippet compiles with zero warnings and prints `ok` (step 5)
- [ ] `git diff -- '*.c'` is empty
- [ ] `make clean && make` exits 0 with zero compiler warnings
- [ ] `make test` → `123 Tests 0 Failures 0 Ignored`
- [ ] Both `diff` commands in step 7 produce no output
- [ ] `git status --porcelain` lists only `tikhonov.h`, `savgol.h`, and `plans/README.md`

## STOP conditions

Stop and report back (do not improvise) if:

- `tikhonov.h:91-123` or `savgol.h:69-82` does not match the "Current state"
  excerpts — the files drifted since this plan was written.
- Any claim in the "Claims that are TRUE" table turns out to be false. That is a
  finding against the *code*, not the docs, and belongs in its own plan. Report it
  and stop; do not change the code to match the comment, and do not change the
  comment to match a behaviour you have not verified.
- Any `.c` file shows up in `git diff`.
- The example still fails to compile after your edit and you cannot see why in two
  attempts.
- `make test` reports anything other than 123 passing.

## Maintenance notes

- **What a reviewer should scrutinise**: that step 5 was actually run and the
  compiler was actually silent. "The example looks right" is exactly the review
  that let this finding survive five audits.
- **Why the example collapses two methods into one**: the old block demonstrated
  manual and automatic lambda by declaring `result` twice. Any future edit that
  re-adds a second `TikhonovResult *result` in the same scope reintroduces the
  compile error.
- **Root cause, unaddressed**: nothing prevents this drift from recurring. The
  headers are prose; the compiler never reads them. If signatures change again,
  these examples rot again. The cheap partial mitigation is the rule "if you
  change a public signature, grep the headers for its name" — worth adding to the
  `smooth-dev-tasks` skill if this happens a third time.
- **Related, out of scope**: `doc/audit-tikhonov.md` still lists TK6 as OPEN, and
  also lists TK3 as OPEN even though the L-curve code it describes was deleted in
  `tikhonov.c` V5.4. Both are tracked in `plans/README.md`.
