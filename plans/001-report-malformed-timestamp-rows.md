# Plan 001: Report malformed timestamp rows instead of dropping them silently

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md`.
>
> **Drift check (run first)**: `git diff --stat 86764a9..HEAD -- parser.c revision.h`
> If either file changed since this plan was written, compare the "Current state"
> excerpts against the live code before proceeding; on a mismatch, treat it as a
> STOP condition.

## Status

- **Priority**: P1
- **Effort**: S
- **Risk**: LOW
- **Depends on**: none
- **Category**: bug
- **Planned at**: commit `86764a9`, 2026-07-27

## Why this matters

In timestamp mode (`-T`), `parser.c` silently discards any input row where the
timestamp column exists but the second whitespace token is missing — a bare
`2026-01-01`, a label like `BADROW`, or a truncated final line. No counter is
incremented and no message is printed. This is the only `continue` branch in the
whole file that drops a data row without accounting for it.

The failure is worse than plain data loss, because dropping rows changes the x
grid. Reproduced on `86764a9` with a 7-row file containing 2 malformed rows: the
program emits 5 data rows and the *only* thing the user sees is

```
# WARNING: Significant spacing variation detected: CV = 0.38
```

The grid became non-uniform *because* rows were dropped, so the program blames
the user's sampling instead of admitting it discarded input. A user who trusts
that warning will go resample data that was fine.

This is the same class of defect the project already fixed as A1 in v5.11.42
(`y` silently shifting against `x` when a timestamp was dropped). The remaining
sibling path was missed because timestamp-mode parsing has no end-to-end tests —
that gap is plan 002, which also carries the regression test for this fix.

After this plan: every dropped row is counted and reported, using the same
mechanism and output stream as the existing `skipped_nonnumeric` counter. **Which
rows get dropped does not change** — only whether the user is told.

## Current state

Files in scope, and their role:

- `parser.c` — the input parser. `parse_input()` is the only function; the
  timestamp-mode branch runs from line 91 to line 191.
- `revision.h` — version string plus the version-history comment block at the top.

### The three sites you will touch

**Site 1** — counter declarations, `parser.c:22-26`:

```c
  char line[MAX_LINE];
  int line_number = 0;
  int skipped_nonnumeric = 0;
  int n = 0;
  int abuf = 0;
```

**Site 2** — the silent drop, `parser.c:132-143`:

```c
      /* Detect timestamp format and assemble timestamp_str */
      int ts_token_count;
      if (strchr(tokens[ts_tok_start], 'T') != NULL) {
        ts_token_count = 1;
        snprintf(timestamp_str, sizeof(timestamp_str), "%s", tokens[ts_tok_start]);
      } else if (ts_tok_start + 1 < ntok) {
        ts_token_count = 2;
        snprintf(timestamp_str, sizeof(timestamp_str), "%s %s",
                 tokens[ts_tok_start], tokens[ts_tok_start + 1]);
      } else {
        continue;  /* malformed: space-format expects two tokens, only one present */
      }
```

The `else` on line 141 is the bug. Compare with the correctly-accounted sibling
30 lines below, `parser.c:161-164`, which is the pattern to copy:

```c
      if (endptr != y_tok_end || errno != 0 || isnan(y_value) || isinf(y_value)) {
        skipped_nonnumeric++;
        continue;  /* token is not a fully numeric finite value */
      }
```

**Site 3** — the reporting block, `parser.c:276-284`:

```c
  if (skipped_nonnumeric > 0) {
    if (timestamp_mode) {
      printf("# Skipped %d data row(s) with non-numeric or NaN/Inf value in column %d (y)\n",
             skipped_nonnumeric, y_column);
    } else {
      printf("# Skipped %d data row(s) with non-numeric or NaN/Inf value in column %d (x) or %d (y)\n",
             skipped_nonnumeric, x_column, y_column);
    }
  }
```

### Conventions this repo enforces — follow them

**Diagnostic output convention.** Documented in `CLAUDE.md` and spelled out in
the header comment of `butterworth.c:7-13`:

```
 *   stdout "# ..."        — info that should be preserved with the saved data
 *                           (affects interpretation of the result: selected fc,
 *                           grid CV, numerical-quality warnings).
 *   stderr "Warning: ..." — runtime/operational warnings irrelevant to the
 *                           filtered output (memory usage, etc.).
 *   stderr "ERROR: ..."   — hard failures; function returns NULL/non-zero.
```

A dropped row changes how the result must be interpreted (it moved CV from
uniform to 0.38 in the reproduction above), so the new message goes to
**stdout with a `#` prefix** — same as its sibling `skipped_nonnumeric`.

**Byte-identical output when nothing is wrong.** This project verifies releases
by diffing program output byte-for-byte against the previous version (see the
v5.11.44–v5.11.46 entries in `revision.h`). Your change must not alter a single
byte of output for input that has no malformed rows. That is why step 3 adds a
*separate* `if` block rather than extending the existing message.

**Version history style.** `revision.h` opens with a comment block, newest entry
first. Each entry states what was wrong, why it mattered, what the fix does, and
what was verified. Read the `v5.11.47` entry at the top of `revision.h` as your
model — match its density and tone.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build | `make` | exit 0, no compiler warnings |
| Unit tests | `make test` | `116 Tests 0 Failures 0 Ignored` / `OK` |
| Leak check | `make test-valgrind` | exit 0 |
| Clean rebuild | `make clean && make` | exit 0, no warnings |

The build uses `clang` with `-Wall -Wextra -pedantic -O2` and requires
`-llapack -lblas`. Both are already present on this machine; a clean build on
`86764a9` produces **zero** warnings. If your change produces any warning, fix it
before proceeding.

## Scope

**In scope** (the only files you may modify):

- `parser.c`
- `revision.h`
- `plans/README.md` (status row only)

**Out of scope** (do NOT touch, even though they look related):

- `timestamp.c` / `timestamp.h` — `convert_timestamps_to_relative()` already
  compacts the parallel `y` array correctly (fixed as A1 in v5.11.42) and reports
  its own invalid-timestamp count via `ts_ctx->errors_encountered`. It is not
  involved in this bug.
- The normal (non-`-T`) parsing branch, `parser.c:193-273` — it has no equivalent
  gap; every drop there is already counted.
- Any change to **which** rows are dropped. Do not convert the `continue` into
  `goto fail`. Files that parse today must still parse after this change; a
  truncated final line is common in real data and must stay tolerated.
- `tests/` — the regression test for this fix belongs to plan 002.
- `README.md` — no user-facing documentation change is needed; the new line is
  self-explanatory and matches an already-documented sibling.

## Git workflow

- Branch: `advisor/001-report-malformed-timestamp-rows`
- One commit for the whole change.
- Message style follows `git log`: a `fix:` prefix and the version in parens,
  e.g. `fix: prefix every line of multi-line grid warnings (v5.11.47)`.
  Use: `fix: report malformed timestamp rows instead of dropping them (v5.11.48)`
- Do NOT push or open a PR.

## Steps

### Step 1: Capture the before-baseline

Create the reproduction fixture and record current behaviour. You will diff
against this in step 5.

```bash
mkdir -p /tmp/plan001
cat > /tmp/plan001/malformed.dat <<'EOF'
2026-01-01 00:00:00 1
2026-01-01 00:00:01 2
BADROW
2026-01-01 00:00:03 4
alsobad
2026-01-01 00:00:05 6
2026-01-01 00:00:06 7
EOF
make >/dev/null
./smooth -T -m0 -n3 -p1 /tmp/plan001/malformed.dat 2>&1 | grep -ci skip
./smooth -T -m0 -n3 -p1 /tmp/plan001/malformed.dat 2>&1 | grep -c '^2026'
```

**Verify**: first command prints `0` (no diagnostic today), second prints `5`
(7 rows in, 5 out). If the first prints anything other than `0`, the bug is
already fixed — STOP and report.

Also capture a clean-input baseline for the byte-identity check. `test_data.txt`
is tracked in git (259 lines, uniform x grid); the timestamp fixture is generated
here because the repo has no tracked timestamp data file.

```bash
for k in $(seq 0 9); do printf '2026-01-01 00:00:%02d %d.0\n' $k $((10+k)); done \
  > /tmp/plan001/clean_ts.dat

./smooth -T -m0 -n3 -p1 /tmp/plan001/clean_ts.dat > /tmp/plan001/before_T.txt 2>&1
./smooth -m2 -l 0.1 test_data.txt                 > /tmp/plan001/before_tik.txt 2>&1
./smooth -m0 -n5 -p2 -d test_data.txt             > /tmp/plan001/before_poly.txt 2>&1
```

**Verify**: all three files are non-empty (`wc -l /tmp/plan001/before_*.txt`), and
the clean timestamp run produces no skip message:

```bash
grep -ci skip /tmp/plan001/before_T.txt
```

→ `0`. If it is not 0, the "clean" fixture is not clean — STOP and report.

### Step 2: Add the counter

In `parser.c`, add one declaration immediately after the existing
`int skipped_nonnumeric = 0;` (site 1, line 24):

```c
  int skipped_malformed_ts = 0;
```

**Verify**: `grep -n "skipped_malformed_ts" parser.c` → exactly 1 match.

### Step 3: Count the dropped row

In `parser.c`, at site 2 (the `else` branch at line 141), add the increment.
Keep the existing comment.

```c
      } else {
        skipped_malformed_ts++;
        continue;  /* malformed: space-format expects two tokens, only one present */
      }
```

**Verify**: `grep -n "skipped_malformed_ts" parser.c` → exactly 2 matches;
`make` → exit 0 with no warnings.

### Step 4: Report the count

In `parser.c`, immediately **after** the closing brace of the existing
`if (skipped_nonnumeric > 0) { ... }` block (site 3, after line 284), add a
separate block:

```c
  if (skipped_malformed_ts > 0) {
    printf("# Skipped %d data row(s) with malformed timestamp in column %d\n",
           skipped_malformed_ts, x_column);
  }
```

Do **not** merge this into the existing `if`, and do **not** guard it with
`if (timestamp_mode)` — the counter can only be non-zero inside the timestamp
branch, so an unguarded `if` is correct and keeps normal-mode output untouched.

**Verify**:

```bash
make >/dev/null && ./smooth -T -m0 -n3 -p1 /tmp/plan001/malformed.dat 2>&1 | grep -i skip
```

→ exactly one line, reading:

```
# Skipped 2 data row(s) with malformed timestamp in column 1
```

And the row count is unchanged:

```bash
./smooth -T -m0 -n3 -p1 /tmp/plan001/malformed.dat 2>&1 | grep -c '^2026'
```

→ `5`. If this is no longer 5, you changed which rows are dropped — STOP.

### Step 5: Confirm byte-identical output on clean input

```bash
./smooth -T -m0 -n3 -p1 /tmp/plan001/clean_ts.dat > /tmp/plan001/after_T.txt 2>&1
./smooth -m2 -l 0.1 test_data.txt                 > /tmp/plan001/after_tik.txt 2>&1
./smooth -m0 -n5 -p2 -d test_data.txt             > /tmp/plan001/after_poly.txt 2>&1
diff /tmp/plan001/before_T.txt    /tmp/plan001/after_T.txt
diff /tmp/plan001/before_tik.txt  /tmp/plan001/after_tik.txt
diff /tmp/plan001/before_poly.txt /tmp/plan001/after_poly.txt
```

**Verify**: all three `diff` commands produce **no output** and exit 0. Any
difference means the new `if` fired on clean input — STOP and report.

### Step 6: Run the full suite and leak check

```bash
make test
make test-valgrind
```

**Verify**: `make test` prints `116 Tests 0 Failures 0 Ignored` then `OK`.
`make test-valgrind` exits 0. The test count stays 116 — the regression test is
plan 002's job. If any previously-passing test now fails, STOP.

### Step 7: Bump the version and write the history entry

In `revision.h`:

1. Line 413: `#define VERSION "5.11.47"` → `#define VERSION "5.11.48"`
2. Line 414: set `#define REVDATE` to today's date in `YYYY-MM-DD` form.
3. Insert a new entry at the top of the version-history comment block, above the
   current `v5.11.47` entry, and change the existing `v5.11.47 (current):` marker
   to plain `v5.11.47:`. The new entry must cover:
   - what was wrong: the `else` branch at `parser.c:141` dropped a row in `-T`
     mode without counting it — the only unaccounted drop in the file;
   - which inputs hit it: timestamp column present but the second whitespace
     token missing (bare date, label, truncated final line);
   - why it was invisible and misleading: no message at all, and the resulting
     gap in the x grid surfaced as a spurious `Significant spacing variation
     detected: CV = 0.38` warning that blamed the user's sampling;
   - the fix: new `skipped_malformed_ts` counter reported on stdout as
     `# Skipped N data row(s) with malformed timestamp in column M`, mirroring
     `skipped_nonnumeric`; which rows are dropped is unchanged;
   - verification: three program outputs byte-identical on clean input, 116 tests
     pass, zero valgrind leaks;
   - that the end-to-end regression guard lands with the timestamp parser tests.

**Verify**:

```bash
grep -n 'define VERSION' revision.h
make clean && make && ./smooth -h 2>&1 | grep -i version
```

→ both report `5.11.48`, and the clean build produces no warnings.

### Step 8: Update the plan index

Set this plan's status to `DONE` in the table in `plans/README.md`.

**Verify**: `grep -n "001" plans/README.md` shows `DONE`.

## Test plan

No new tests in this plan — the end-to-end regression test
(`test_parser_ts_malformed_row_is_reported`) is step 2 of plan 002, which depends
on this one. Splitting them this way keeps `make test` green at every commit,
which `CLAUDE.md` requires.

Your verification here is therefore: the manual reproduction in steps 1 and 4,
the byte-identity diff in step 5, and the unchanged 116-test suite in step 6.

## Done criteria

ALL must hold:

- [ ] `make clean && make` exits 0 with zero compiler warnings
- [ ] `./smooth -T -m0 -n3 -p1 /tmp/plan001/malformed.dat 2>&1 | grep -c '^2026'` → `5`
- [ ] The same run prints exactly `# Skipped 2 data row(s) with malformed timestamp in column 1`
- [ ] All three `diff` commands in step 5 produce no output
- [ ] `make test` → `116 Tests 0 Failures 0 Ignored`
- [ ] `make test-valgrind` exits 0
- [ ] `grep -c 'define VERSION "5.11.48"' revision.h` → `1`, and a matching history entry exists
- [ ] `git status --porcelain` lists only `parser.c`, `revision.h`, `plans/README.md`

## STOP conditions

Stop and report back (do not improvise) if:

- The code at `parser.c:132-143` or `parser.c:276-284` does not match the
  "Current state" excerpts — the file has drifted since this plan was written.
- Step 1's baseline shows a skip message already present (bug already fixed).
- The malformed-row output count changes from 5 to anything else — you altered
  drop behaviour, which is explicitly out of scope.
- Any `diff` in step 5 shows a difference on clean input.
- Any of the 116 existing tests fails, or `make test-valgrind` reports a
  definite/indirect leak.
- The fix appears to require touching `timestamp.c` or the normal-mode branch.

## Maintenance notes

- **What a reviewer should scrutinise**: that step 4 added a *separate* `if`
  rather than folding the message into the `skipped_nonnumeric` block. Folding
  them would change output for files that have only non-numeric skips, breaking
  the byte-identity property this project relies on for release verification.
- **Future interaction**: if anyone later adds a `-q`/quiet flag (audit finding
  B8 — the `-l auto` GCV log floods stdout), this message is data-relevant and
  must survive at the default verbosity. It reports lost input, not progress.
- **Deliberately deferred**: making the drop a hard error. It would be more
  consistent with the neighbouring insufficient-column branches, which use
  `goto fail`, but it would break files that parse today — a truncated final line
  is normal in real data. Chosen deliberately; do not "fix" this later without a
  decision.
- **Related open finding**: `parser.c` is still absent from `TEST_MODULES` in the
  `Makefile`, so it is only ever exercised through the `./smooth` binary. Plan
  002 works within that constraint rather than changing it.
