# Plan 002: Add end-to-end tests for the timestamp-mode input parser

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md`.
>
> **Drift check (run first)**: `git diff --stat 86764a9..HEAD -- tests/test_parser.c tests/test_main.c CLAUDE.md parser.c`
> Changes to `parser.c` from plan 001 are expected. If `tests/test_parser.c` or
> `tests/test_main.c` changed, compare the "Current state" excerpts against the
> live code before proceeding; on a mismatch, treat it as a STOP condition.

## Status

- **Priority**: P1
- **Effort**: S–M
- **Risk**: LOW
- **Depends on**: `plans/001-report-malformed-timestamp-rows.md` (must be DONE)
- **Category**: tests
- **Planned at**: commit `86764a9`, 2026-07-27

## Why this matters

The timestamp branch of `parse_input()` (`parser.c:91-191`, roughly 100 lines) has
zero test coverage. All seven existing parser tests run in normal mode —
`grep -rn '\-T' tests/*.c` returns nothing. `parser.c` is also absent from
`TEST_MODULES` in the `Makefile`, so it is never linked into the unit-test binary
and is only ever exercised indirectly through `./smooth`.

That gap is not theoretical: it is exactly why the silent data-loss bug fixed by
plan 001 survived several full-codebase audits. The untested code includes the
trickiest logic in the parser — the mapping from user-facing *logical* columns to
whitespace *tokens*, which must compensate for a timestamp occupying one token
(`2026-01-01T00:00:00`) or two (`2026-01-01 00:00:00`).

The tests in `tests/test_timestamp.c` do not cover this. They test
`parse_timestamp()` and `convert_timestamps_to_relative()` directly — the
functions the parser *calls*, not the parser logic around them.

After this plan: seven new end-to-end tests pin the timestamp path, including a
regression test for plan 001's bug and coverage of all three branches of the
logical-column mapping.

## Current state

Files in scope:

- `tests/test_parser.c` — the seven existing end-to-end parser tests. Drives
  `./smooth` via `popen()` against fixtures written under `/tmp`.
- `tests/test_main.c` — the Unity runner. Declares every test function at the top
  and calls `RUN_TEST()` for each in `main()`.
- `CLAUDE.md` — records the test count per module.

`tests/test_parser.c` is already listed in `TEST_SRC` (`Makefile:17`), so **no
Makefile change is needed**.

### The existing harness in `tests/test_parser.c`

Header comment, `tests/test_parser.c:1-17`:

```c
/* test_parser.c - End-to-end tests for the smooth input parser.
 *
 * Drives the ./smooth binary via popen() on small generated fixtures
 * under /tmp, parses data rows from its stdout, and verifies column
 * counts, values, and the "# Skipped N ..." summary message.
 *
 * Requires that ./smooth is built (Makefile target `test` adds the
 * dependency).
 */
...
#define SMOOTH_BIN "./smooth"
```

Result struct and fixture writer, `tests/test_parser.c:19-32`:

```c
typedef struct {
    int data_rows;
    int has_skip_msg;
    int skip_count;
    double first_x, first_y;
    double last_x,  last_y;
} SmoothRun;

static void write_fixture(const char *path, const char *content) {
    FILE *f = fopen(path, "w");
    TEST_ASSERT_NOT_NULL_MESSAGE(f, "could not open fixture file for writing");
    fputs(content, f);
    fclose(f);
}
```

The existing runner, `tests/test_parser.c:34-67`:

```c
static SmoothRun run_smooth(const char *args, const char *fixture) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), SMOOTH_BIN " %s %s 2>&1", args, fixture);
    FILE *p = popen(cmd, "r");
    TEST_ASSERT_NOT_NULL_MESSAGE(p, "popen failed");

    SmoothRun r = {0};
    r.skip_count = -1;

    char ln[1024];
    while (fgets(ln, sizeof(ln), p)) {
        if (ln[0] == '#') {
            int n;
            if (sscanf(ln, "# Skipped %d data row", &n) == 1) {
                r.has_skip_msg = 1;
                r.skip_count = n;
            }
            continue;
        }
        ...
        double xv, yv;
        if (sscanf(ln, "%lf %lf", &xv, &yv) == 2) {
            ...
        }
    }
    pclose(p);
    return r;
}
```

**`run_smooth()` cannot be reused for `-T` output, and must not be modified.**
Two independent reasons:

1. Its row parser is `sscanf(ln, "%lf %lf", &xv, &yv)`. Timestamp output rows look
   like `2026-01-01 00:00:00           10`. `%lf` on `"2026-01-01"` consumes
   `2026`, then `%lf` on `"-01-01"` yields `-1` — it returns 2 and silently
   records garbage instead of failing.
2. Its skip matcher is `sscanf(ln, "# Skipped %d data row", &n)`. Plan 001's new
   message — `# Skipped 2 data row(s) with malformed timestamp in column 1` —
   also matches that prefix, so the two counters would be conflated.

Changing `run_smooth()` would put the seven existing green tests at risk. Add a
sibling helper instead.

### Test registration in `tests/test_main.c`

Declarations, `tests/test_main.c:158-164`:

```c
void test_parser_iso_t_timestamp_one_column(void);
void test_parser_iso_space_timestamp_two_columns(void);
void test_parser_partial_numeric_token_is_placeholder(void);
void test_parser_nan_in_y_skips_row(void);
void test_parser_inf_in_x_skips_row(void);
void test_parser_label_outside_xy_is_harmless(void);
void test_parser_comments_are_stripped(void);
```

Runner block, `tests/test_main.c:424-435`:

```c
    printf("========================================\n");
    printf("Running tests for input parser (end-to-end via ./smooth)\n");
    printf("========================================\n\n");

    RUN_TEST(test_parser_iso_t_timestamp_one_column);
    RUN_TEST(test_parser_iso_space_timestamp_two_columns);
    RUN_TEST(test_parser_partial_numeric_token_is_placeholder);
    RUN_TEST(test_parser_nan_in_y_skips_row);
    RUN_TEST(test_parser_inf_in_x_skips_row);
    RUN_TEST(test_parser_label_outside_xy_is_harmless);
    RUN_TEST(test_parser_comments_are_stripped);
```

### The logical-column model you are testing

Documented at `parser.c:92-95` and implemented at `parser.c:146-150`:

```c
      /* Map logical y_column to whitespace-token index. Logical columns before
       * the timestamp are unaffected by its width; columns after shift by
       * ts_token_count - 1. y_column != x_column is enforced at -k parse time. */
      int y_token_idx = (y_column < x_column)
                        ? y_column - 1
                        : y_column - 1 + (ts_token_count - 1);
```

The timestamp counts as **one logical column** regardless of whether it occupies
one or two whitespace tokens. Concretely, for a row with a label between the
timestamp and the value:

| Input row | Correct flag |
|---|---|
| `2026-01-01T00:00:00 expA 10.0` | `-T -k 1:3` |
| `2026-01-01 00:00:00 expA 10.0` | `-T -k 1:3` — **the same flag** |

This is a trap worth stating plainly, because the token counts differ (3 vs 4):
`-T -k 1:4` on the space-separated form fails with
`ERROR: Line 1 has insufficient columns for y column 4`. Test 4 below exists
specifically to pin this equivalence.

Note that the two *existing* tests `test_parser_iso_t_timestamp_one_column` and
`test_parser_iso_space_timestamp_two_columns` use `-k 2:5` / `-k 3:6` — those run
in **normal mode** (no `-T`), where each timestamp token is an ordinary
non-numeric placeholder column. Different model; do not copy their flags.

### Conventions to match

- Test naming: `test_parser_<what>`, descriptive, snake_case.
- Each test: `write_fixture()` → run → asserts → `remove(path)`.
- Fixtures live at hardcoded `/tmp/test_parser_*.dat` paths. Keep that convention
  for consistency with the existing seven; do not switch to `mkstemp`/`TMPDIR` in
  this plan.
- Doubles are compared with `TEST_ASSERT_DOUBLE_WITHIN`, never `==`
  (`tests/README.md:182-183`).
- A short comment above each test explaining the non-obvious case, matching the
  style already used in `tests/test_parser.c:69-71`.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build binary + tests | `make test` | `123 Tests 0 Failures 0 Ignored` / `OK` after this plan |
| Leak check | `make test-valgrind` | exit 0 |
| Clean rebuild | `make clean && make test` | exit 0, no warnings |

Test sources compile with `-Wall -Wextra -pedantic -I. -Itests
-DUNITY_INCLUDE_DOUBLE`. Any new warning must be fixed before you finish.

## Scope

**In scope** (the only files you may modify):

- `tests/test_parser.c` — add one helper and seven tests
- `tests/test_main.c` — seven declarations and seven `RUN_TEST()` calls
- `CLAUDE.md` — test counts
- `plans/README.md` (status row only)

**Out of scope** (do NOT touch):

- `parser.c` and every other source file. This plan adds tests only. If a test
  reveals a *second* bug, record it and STOP — do not fix it here.
- `Makefile` — `tests/test_parser.c` is already in `TEST_SRC` (line 17). In
  particular, do **not** add `parser.o` to `TEST_MODULES`; linking the parser for
  in-process unit tests is a separate, larger decision.
- `run_smooth()` and the seven existing tests in `tests/test_parser.c` — see the
  two reasons above.
- `revision.h` — test-only changes do not bump the version in this repo. See
  commits `902d74e` and `a3be811` for precedent.
- `tests/README.md` — its module table is already stale for several modules;
  fixing it is a separate cleanup.

## Git workflow

- Branch: `advisor/002-timestamp-parser-e2e-tests`
- One commit.
- Message style follows `git log`; test-only commits use no version suffix:
  `test: cover the timestamp-mode parser path end-to-end`
- Do NOT push or open a PR.

## Steps

### Step 1: Confirm plan 001 has landed

```bash
grep -c "malformed timestamp" parser.c
```

**Verify**: prints `1` (the message string plan 001 added). If it prints `0`,
plan 001 has not landed — STOP and report; test 6 below cannot pass.

Record the baseline count:

```bash
make test 2>&1 | tail -3
```

**Verify**: `116 Tests 0 Failures 0 Ignored` / `OK`.

### Step 2: Add the `run_smooth_ts()` helper

In `tests/test_parser.c`, insert this immediately after `run_smooth()` (which
ends at line 67). Leave `run_smooth()` untouched.

```c
/* Timestamp-mode counterpart of run_smooth().
 *
 * A -T output row is "<timestamp tokens...> <y>", so the x field cannot be read
 * with sscanf("%lf %lf") — "2026-01-01" would parse as 2026 followed by -1.
 * Instead the last whitespace token is parsed as y, and a row is only accepted
 * when it starts with a digit (timestamps do; "Warning:"/"ERROR:" lines do not).
 *
 * The two "# Skipped ..." messages share a prefix, so they are told apart by
 * substring, not by sscanf format.
 */
typedef struct {
    int    data_rows;
    int    skip_nonnumeric;   /* -1 if the message was absent */
    int    skip_malformed;    /* -1 if the message was absent */
    double first_y, last_y;
} SmoothRunTs;

static SmoothRunTs run_smooth_ts(const char *args, const char *fixture) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), SMOOTH_BIN " %s %s 2>&1", args, fixture);
    FILE *p = popen(cmd, "r");
    TEST_ASSERT_NOT_NULL_MESSAGE(p, "popen failed");

    SmoothRunTs r = {0};
    r.skip_nonnumeric = -1;
    r.skip_malformed  = -1;

    char ln[1024];
    while (fgets(ln, sizeof(ln), p)) {
        if (ln[0] == '#') {
            int k;
            if (sscanf(ln, "# Skipped %d", &k) == 1) {
                if (strstr(ln, "malformed timestamp") != NULL) {
                    r.skip_malformed = k;
                } else {
                    r.skip_nonnumeric = k;
                }
            }
            continue;
        }

        char *q = ln;
        while (*q == ' ' || *q == '\t') q++;
        if (*q < '0' || *q > '9') continue;  /* not a data row */

        /* Find the last whitespace-separated token and require it to be a
         * complete number. */
        size_t len = strlen(ln);
        while (len > 0 && (ln[len-1] == '\n' || ln[len-1] == '\r' ||
                           ln[len-1] == ' '  || ln[len-1] == '\t')) {
            ln[--len] = '\0';
        }
        if (len == 0) continue;
        char *tok = ln + len;
        while (tok > ln && tok[-1] != ' ' && tok[-1] != '\t') tok--;

        char *endptr;
        double yv = strtod(tok, &endptr);
        if (*endptr != '\0') continue;  /* trailing junk: not a data row */

        if (r.data_rows == 0) r.first_y = yv;
        r.last_y = yv;
        r.data_rows++;
    }
    pclose(p);
    return r;
}
```

**Verify**: `make test 2>&1 | tail -3` still reports
`116 Tests 0 Failures 0 Ignored` and the build emits no warnings. (The helper is
unused so far; `static` + unused would warn, so if `-Wunused-function` fires,
proceed straight to step 3 and re-check there.)

### Step 3: Add the seven tests

Append these to `tests/test_parser.c`, after the last existing test
(`test_parser_label_outside_xy_is_harmless`, ends line 203). All expected values
below were measured by running `./smooth` — they are facts, not guesses.

The `y` values `10 11 12 13 14` are exactly reproduced by `-m0 -n3 -p1` because
the data is exactly linear, so a degree-1 fit is exact. Do not change the
smoothing flags without re-deriving the expectations.

**Test 1 — T-separated timestamp, default columns** (`ts_token_count = 1`):

- fixture `/tmp/test_parser_ts_t.dat`, five rows
  `2026-01-01T00:00:0<k> 1<k>.0` for k = 0..4
- args `-T -m0 -n3 -p1`
- assert: `data_rows == 5`, `first_y ≈ 10.0`, `last_y ≈ 14.0`,
  `skip_nonnumeric == -1`, `skip_malformed == -1`

**Test 2 — space-separated timestamp, default columns** (`ts_token_count = 2`):

- fixture `/tmp/test_parser_ts_space.dat`, same data with
  `2026-01-01 00:00:0<k> 1<k>.0`
- args `-T -m0 -n3 -p1`
- assert: identical to test 1

**Test 3 — `-k` maps y past a T-separated timestamp** (`y_column > x_column`,
shift 0):

- fixture `/tmp/test_parser_ts_k_t.dat`, rows `2026-01-01T00:00:0<k> expA 1<k>.0`
- args `-T -k 1:3 -m0 -n3 -p1`
- assert: `data_rows == 5`, `first_y ≈ 10.0`, `last_y ≈ 14.0`, no skips

**Test 4 — the same `-k 1:3` works for the space-separated form**
(`y_column > x_column`, shift +1). This is the test that pins the logical-column
abstraction; give it a comment saying so, and noting that `-k 1:4` would fail.

- fixture `/tmp/test_parser_ts_k_space.dat`, rows
  `2026-01-01 00:00:0<k> expA 1<k>.0`
- args `-T -k 1:3 -m0 -n3 -p1` — **the same flag as test 3**
- assert: `data_rows == 5`, `first_y ≈ 10.0`, `last_y ≈ 14.0`, no skips

**Test 5 — y before the timestamp** (`y_column < x_column` branch):

- fixture `/tmp/test_parser_ts_y_first.dat`, rows `1<k>.0 2026-01-01 00:00:0<k>`
- args `-T -k 2:1 -m0 -n3 -p1`
- assert: `data_rows == 5`, `first_y ≈ 10.0`, `last_y ≈ 14.0`, no skips

**Test 6 — regression for plan 001: malformed rows are reported.** Comment it as
the regression guard for the silent-drop bug.

- fixture `/tmp/test_parser_ts_malformed.dat`:

```
2026-01-01 00:00:00 1
2026-01-01 00:00:01 2
BADROW
2026-01-01 00:00:03 4
alsobad
2026-01-01 00:00:05 6
2026-01-01 00:00:06 7
```

- args `-T -m0 -n3 -p1`
- assert: `data_rows == 5` (7 in, 2 dropped) **and** `skip_malformed == 2`.
  Assert both — the count alone would pass even if the message vanished, and the
  message alone would pass if drop behaviour changed.

**Test 7 — non-numeric y in timestamp mode is reported** (the
`skipped_nonnumeric` path inside the `-T` branch, currently untested):

- fixture `/tmp/test_parser_ts_nan.dat`, six rows
  `2026-01-01 00:00:0<k> 1<k>.0` for k = 0..5, with the k = 2 row's value
  replaced by `NaN`
- args `-T -m0 -n3 -p1`
- assert: `data_rows == 5`, `skip_nonnumeric == 1`, `skip_malformed == -1`

Every test ends with `remove(path);`, matching the existing seven.

**Verify**: `make test 2>&1 | tail -3` — still `116 Tests` at this point, because
the new tests are not registered yet. The build must be warning-free.

### Step 4: Register the tests

In `tests/test_main.c`, add the seven declarations after line 164 (the existing
parser declarations block), and the seven `RUN_TEST()` calls after line 435 (the
existing parser runner block), in the same order as step 3.

**Verify**:

```bash
grep -c "RUN_TEST(test_parser" tests/test_main.c
make test 2>&1 | tail -3
```

→ first prints `14`; second prints `123 Tests 0 Failures 0 Ignored` then `OK`.

If any of the seven fails, read its output before changing anything: tests 1–5
and 7 failing points at a real parser defect (STOP and report — fixing `parser.c`
is out of scope); test 6 failing points at plan 001 not having landed correctly.

### Step 5: Leak check

```bash
make test-valgrind
```

**Verify**: exits 0. The new tests allocate nothing themselves, but `popen`/
`pclose` pairing is worth confirming.

### Step 6: Update the recorded test counts

In `CLAUDE.md`, lines 117-118, currently:

```
- 116 tests total: grid_analysis (7), polyfit (21), savgol (16), tikhonov (25),
  butterworth (22), timestamp (18), parser (7). Source of truth is `tests/test_main.c`.
```

Change `116 tests total` → `123 tests total` and `parser (7)` → `parser (14)`.
Leave every other module count alone.

**Verify**: `grep -n "123 tests total" CLAUDE.md` → 1 match;
`grep -n "parser (14)" CLAUDE.md` → 1 match.

### Step 7: Update the plan index

Set this plan's status to `DONE` in `plans/README.md`.

**Verify**: `grep -n "002" plans/README.md` shows `DONE`.

## Test plan

This plan *is* the test plan. Summary of what the seven new tests pin:

| Test | Parser path covered |
|---|---|
| 1 | `parser.c:134-136` — `ts_token_count = 1` |
| 2 | `parser.c:137-140` — `ts_token_count = 2` |
| 3 | `parser.c:148-150` — `y_column > x_column`, shift 0 |
| 4 | `parser.c:148-150` — `y_column > x_column`, shift +1 |
| 5 | `parser.c:148-150` — `y_column < x_column` |
| 6 | `parser.c:141-143` — malformed drop is counted (**regression, plan 001**) |
| 7 | `parser.c:159-164` + `276-284` — non-numeric y in `-T` mode |

Structural pattern to follow: the existing tests in the same file, e.g.
`test_parser_nan_in_y_skips_row` (`tests/test_parser.c:129-145`).

Still uncovered after this plan, deliberately: the `MAX_LINE` truncation path
(`parser.c:54-89`), the `MAX_COLS` overflow paths, and the realloc growth beyond
`BUF = 512` rows. All three need multi-kilobyte fixtures; worth a follow-up, not
worth inflating this one.

## Done criteria

ALL must hold:

- [ ] `make clean && make test` exits 0 with zero compiler warnings
- [ ] `make test` prints `123 Tests 0 Failures 0 Ignored` then `OK`
- [ ] `grep -c "RUN_TEST(test_parser" tests/test_main.c` → `14`
- [ ] `make test-valgrind` exits 0
- [ ] `grep -c "123 tests total" CLAUDE.md` → `1` and `grep -c "parser (14)" CLAUDE.md` → `1`
- [ ] `git diff --stat -- parser.c` is empty (no source change in this plan)
- [ ] `git status --porcelain` lists only `tests/test_parser.c`, `tests/test_main.c`, `CLAUDE.md`, `plans/README.md`

## STOP conditions

Stop and report back (do not improvise) if:

- Step 1 shows plan 001 has not landed (`grep -c "malformed timestamp" parser.c` → 0).
- `tests/test_parser.c:34-67` or `tests/test_main.c:158-164` does not match the
  excerpts above — the files drifted since this plan was written.
- Any of the seven existing parser tests, or any other of the 116, starts failing.
- Any of tests 1, 2, 3, 4, 5, or 7 fails. That means the parser has a defect
  beyond the one plan 001 fixed. Report the failing test and its output; do not
  edit `parser.c`.
- You conclude that `run_smooth()` must be modified to make a test pass.
- A test needs `parser.o` linked into `TEST_MODULES` to work.

## Maintenance notes

- **What a reviewer should scrutinise**: that `run_smooth()` and the seven
  original tests are byte-identical in the diff, and that test 6 asserts *both*
  `data_rows == 5` and `skip_malformed == 2`. Dropping either assertion makes the
  regression guard hollow.
- **Fragility to watch**: `run_smooth_ts()` classifies a line as data by "starts
  with a digit, last token parses cleanly as a double". If a future diagnostic is
  ever printed to stdout without a `#` prefix and begins with a digit, these tests
  will miscount. The `#` prefix is the project's documented convention
  (`butterworth.c:7-13`), so this is a safe bet, but it is the assumption to
  revisit if these tests ever start behaving oddly.
- **If the skip messages are reworded**, test 6 and test 7 break by design — they
  match on the substring `"malformed timestamp"` and on the absence of it. That
  coupling is intentional: the message is the user-visible contract.
- **Follow-up worth considering, deliberately not done here**: adding `parser.o`
  to `TEST_MODULES` so `parse_input()` can be unit-tested in process against a
  `tmpfile()` stream. That gives stronger assertions directly on `ParseResult`
  (`n`, `x`, `y`, `ts_ctx`) without scraping stdout, at the cost of making the
  `#` message assertions awkward. The end-to-end approach was chosen for this
  plan because it tests what the user actually sees.
- **Known pre-existing wart, out of scope**: fixtures use hardcoded `/tmp` paths
  rather than `TMPDIR` or `mkstemp`, so two concurrent test runs on one machine
  can collide. Inherited from the existing seven tests; changing it should be one
  focused commit covering all fourteen.
