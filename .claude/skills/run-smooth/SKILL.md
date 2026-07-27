---
name: run-smooth
description: Build, run, and drive the `smooth` CLI - a C data-smoothing filter (polyfit, Savitzky-Golay, Tikhonov, Butterworth). Use when asked to run, build, start, smoke-test, or verify smooth, to check that a smoothing-method change actually works on real data, or to run its unit tests and valgrind leak check.
---

# Running `smooth`

`smooth` is a batch Unix filter: a table of `x y` in, a smoothed table out. No
GUI, no server, no interactive surface - so the useful handle is **numeric**,
not a screenshot. `driver.py` proves each method actually reduces noise against
a known clean signal, which is the layer nearly every change here touches.

All paths below are relative to the repo root (`/home/tch/Lang/c/smooth`).

## Prerequisites

Already present in this container. On a bare Debian/Ubuntu box:

```bash
sudo apt-get install -y clang make liblapack-dev libblas-dev valgrind
```

`python3` (stdlib only) is needed for the driver. **numpy is not installed and
is not required** - the driver deliberately uses stdlib `math` only.

## Build

```bash
make
```

`make` alone is correct and sufficient. **Always `make` before `make test`** -
see Gotchas, the reverse order silently produces an unoptimized binary.

The Makefile passes `-I$(HOME)/include -L$(HOME)/lib` by default. Those paths
not existing is harmless (verified with `LIBDIR=/nonexistent`); system LAPACK
is found regardless. Override with `LIBDIR=/path INCDIR=/path make` if needed.

## Run (agent path)

```bash
python3 .claude/skills/run-smooth/driver.py smoke
```

Takes ~2s, exits 0 on success, prints `ALL CHECKS PASSED`. It covers:

- all four methods on a **uniform** grid - asserts RMSE against the clean
  sine *drops* (this is the check that catches broken math that still exits 0)
- polyfit + Tikhonov on a **non-uniform** grid
- Savgol and Butterworth **rejecting** the non-uniform grid with exit 1 and a
  grid-specific stderr message (the "reject loudly" rule in CLAUDE.md)
- `-d` derivatives validated against `cos(x)`
- `-g` grid analysis, `-T` timestamp mode (both separators), `-k` column
  selection, stdin filtering, and four error paths

Expected output:

```
uniform grid — all four methods must smooth:
  [OK] polyfit  -m 0              RMSE 0.01489 -> 0.01018
  [OK] savgol   -m 1              RMSE 0.01489 -> 0.01020
  [OK] tikhonov -m 2              RMSE 0.01489 -> 0.00655
  [OK] butterworth -m 3           RMSE 0.01489 -> 0.00605
...
ALL CHECKS PASSED
```

One line is **expected to report a real, still-open defect** (see Gotchas):

```
stdout table integrity (known defect — see SKILL.md Gotchas):
  [KNOWN BUG] 1 unprefixed diagnostic line(s) on stdout, e.g. 'Adaptive methods may improve results.'
```

The driver does not fail on it. If it ever prints `[OK] every diagnostic line
is '#'-prefixed`, the defect was fixed - delete that check.

### Driving it by hand

Every checked-in `.dat` in this repo is non-uniform (CV = 0.177), so **Savgol
and Butterworth reject all of them**. To exercise those two, generate a
uniform fixture first:

```bash
python3 .claude/skills/run-smooth/driver.py fixture /tmp/uniform.dat
./smooth -m 1 -n 7 -p 2 -d /tmp/uniform.dat     # Savitzky-Golay + derivative
./smooth -m 3 -f auto /tmp/uniform.dat          # Butterworth, auto cutoff
```

Methods are `-m 0` polyfit (default), `1` savgol, `2` tikhonov, `3`
butterworth. Grid-tolerant methods work on the repo's own data:

```bash
./smooth -m 0 -n 7 -p 3 test_data.dat           # polyfit
./smooth -m 2 -l auto -d test_data.dat          # Tikhonov, GCV auto-lambda
./smooth -g test_data.dat                       # grid analysis only
./smooth -T -m 2 -l 0.1 -d example.dat          # timestamped input
cat test_data.dat | ./smooth -m 0 -n 5          # stdin filter
./smooth -h                                     # full option list, exit 0
```

## Test

```bash
make test           # 116 Unity tests, ~1s
make test-valgrind  # same, under valgrind; exits 1 on any definite/indirect leak
```

Both pass clean as of this writing: `116 Tests 0 Failures 0 Ignored`, and
`in use at exit: 0 bytes in 0 blocks`. The binary itself is leak-free too -
`valgrind -q --leak-check=full ./smooth -m 2 -l auto /tmp/uniform.dat` exits 0.

## Gotchas

- **`make test` before `make` gives you an `-O0` binary, silently.**
  `CFLAGS` is a target-specific variable set only on `build` and `debug`.
  `test` depends on `$(PROGRAM)` directly, so from a clean tree it compiles
  every `.o` with **empty** CFLAGS - no `-O2`, no `-Wall -Wextra -pedantic`.
  Worse, a subsequent plain `make` sees those objects as up to date and
  **only relinks them**, so you get an unoptimized `./smooth` that looks
  freshly built and whose warnings never appeared. Always `make` first, or
  `make clean && make && make test`. Same trap for any bare object target
  (`make smooth.o`).
- **Diagnostic lines leak onto stdout unprefixed.** `grid_analysis.c:130-158`
  builds a multi-line `warning_msg`, but only its *first* line gets the
  `# WARNING: ` prefix. Continuation lines land on stdout bare, e.g.
  `Adaptive methods may improve results.` - which corrupts the table for any
  numeric consumer. Triggers whenever CV crosses a warning band (all three
  bands, roughly CV > 0.2). Do not assume `grep -v '^#'` yields clean numbers;
  filter on "parses as floats" instead, as `driver.py:columns()` does.
- **No repo data file is uniform.** All the committed `.dat` files sit at
  CV = 0.177, above Savgol's 0.05 threshold and above Butterworth's 0.15.
  Both reject. Use `driver.py fixture` to get a testable uniform grid.
- **`-T` output field count depends on the separator.** The timestamp is
  echoed verbatim, so `2025-09-25 14:06:06.390` (space form) occupies *two*
  whitespace fields and a `-d` row splits into 4, while the `T` form splits
  into 3. Split from the right, or match on the timestamp prefix.
- **Piping into `head` shows exit 141, not a failure.** That is SIGPIPE from
  `head` closing the pipe. Check `${PIPESTATUS[0]}`, or redirect to a file.
  A genuine rejection is exit 1 with a message on stderr.
- **Errors go to stderr with exit 1; everything else, including all
  diagnostics, goes to stdout** (modulo the leak above). Tikhonov `-l auto`
  prints ~28 lines of GCV trace as `#` comments before the data.
- Method preconditions are enforced, not degraded: Savgol rejects CV > 0.05,
  Butterworth rejects CV > 0.15. That is intended behaviour, not a bug.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ERROR: Savitzky-Golay method not suitable for non-uniform grid!`, exit 1 | Expected on every repo `.dat`. Generate a uniform fixture: `driver.py fixture /tmp/uniform.dat` |
| `ERROR: Butterworth filter not suitable for non-uniform grid (CV=0.1772, threshold=0.1500)!` | Same - needs a uniform grid |
| `ValueError: could not convert string to float: 'Adaptive'` while parsing output | The unprefixed-diagnostic defect above. Skip non-numeric lines |
| `Need more data (n < 5)!`, exit 1 | Input has fewer than 5 usable rows (empty file, or all rows filtered as comments/NaN) |
| `Unknown method number: 9`, exit 1 | `-m` takes 0-3 only |
| `Lambda must be non-negative!`, exit 1 | `-l` needs `>= 0` or the literal `auto` |
| `Failed to open file 'X' (2: No such file or directory)`, exit 1 | Bad path; `-` or no argument reads stdin |
| `# WARNING: Data length (n=200) too short to absorb the filter transient for fc=0.0100` | Butterworth with a very small `-f`; padding hit its `n-1` cap. Use a larger fc or more data |
| Build cannot find `-llapack` | `sudo apt-get install -y liblapack-dev libblas-dev` |
| `./smooth` exists but behaves like an old build | Stale objects from the CFLAGS trap. `make clean && make` |

## Harness

`.claude/skills/run-smooth/driver.py` - stdlib-only Python, no framework.
`smoke` (default) runs every check; `fixture PATH` writes a uniform-grid noisy
sine. Fixtures are deterministic (`random.Random(7)`) and land in a temp dir,
so the driver never writes into the repo.
