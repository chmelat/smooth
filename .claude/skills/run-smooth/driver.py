#!/usr/bin/env python3
"""Drive the `smooth` CLI end-to-end.

`smooth` is a batch filter: stdin/file in, smoothed table out. There is no
interactive surface, so the handle a future agent needs is *numeric* — proof
that each method actually smooths, not just that it exits 0.

Subcommands:
  smoke      build-free end-to-end run of all four methods + modes (default)
  fixture P  write a uniform-grid noisy-sine fixture to path P

Only the Python stdlib is used (no numpy in this container).

Usage:
  python3 .claude/skills/run-smooth/driver.py smoke
  python3 .claude/skills/run-smooth/driver.py fixture /tmp/uniform.dat
"""
import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SMOOTH = ROOT / "smooth"

N = 200
H = 0.05
NOISE = 0.05


def signal(x):
    return math.sin(x)


def make_data(path, uniform=True, seed=7):
    """Noisy sine on a uniform (or jittered) grid. Returns the clean y values."""
    rng = random.Random(seed)
    xs, clean, rows = [], [], []
    for i in range(N):
        x = i * H
        if not uniform:
            # +-25% jitter -> CV well above every method's uniformity threshold
            x += H * 0.5 * (rng.random() - 0.5)
        y0 = signal(x)
        xs.append(x)
        clean.append(y0)
        rows.append(f"{x:.6f} {y0 + NOISE * (rng.random() - 0.5):.6f}")
    Path(path).write_text("\n".join(rows) + "\n")
    return xs, clean


def run(args, expect=0):
    p = subprocess.run([str(SMOOTH)] + args, capture_output=True, text=True)
    if p.returncode != expect:
        raise AssertionError(
            f"smooth {' '.join(args)}: exit {p.returncode}, expected {expect}\n"
            f"--- stderr ---\n{p.stderr}"
        )
    return p


def columns(stdout):
    """Data rows only. Diagnostics are *mostly* `#` comments — but multi-line
    grid warnings leak their continuation lines unprefixed (see leaked_lines),
    so anything non-numeric is skipped rather than trusted."""
    out = []
    for line in stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        try:
            out.append([float(f) for f in line.split()])
        except ValueError:
            continue
    return out


def leaked_lines(stdout):
    """Non-comment, non-numeric lines on stdout — each one corrupts the table
    for any downstream numeric consumer."""
    bad = []
    for line in stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        try:
            [float(f) for f in line.split()]
        except ValueError:
            bad.append(line)
    return bad


def rmse(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)) / len(a))


def check_smooths(name, args, data, clean):
    """The real assertion: output must be closer to the clean signal than the
    noisy input was. Catches broken math that still exits 0."""
    rows = run(args + [data]).stdout
    got = columns(rows)
    assert len(got) == N, f"{name}: got {len(got)} rows, expected {N}"
    ys = [r[1] for r in got]
    noisy = [r[1] for r in columns(Path(data).read_text())]
    before, after = rmse(noisy, clean), rmse(ys, clean)
    assert after < before, f"{name}: RMSE got worse ({before:.5f} -> {after:.5f})"
    print(f"  [OK] {name:<26} RMSE {before:.5f} -> {after:.5f}")
    return got


def main():
    argv = sys.argv[1:] or ["smoke"]
    if not SMOOTH.exists():
        sys.exit(f"{SMOOTH} not built — run `make` in {ROOT} first")

    tmp = Path(tempfile.mkdtemp(prefix="smooth-driver-"))
    if argv[0] == "fixture":
        dest = argv[1] if len(argv) > 1 else "uniform.dat"
        make_data(dest, uniform=True)
        print(f"wrote uniform-grid fixture: {dest}")
        return

    uni, non = str(tmp / "uniform.dat"), str(tmp / "nonuniform.dat")
    _, clean = make_data(uni, uniform=True)
    _, clean_nu = make_data(non, uniform=False)

    print("uniform grid — all four methods must smooth:")
    check_smooths("polyfit  -m 0", ["-m", "0", "-n", "7", "-p", "3"], uni, clean)
    check_smooths("savgol   -m 1", ["-m", "1", "-n", "7", "-p", "2"], uni, clean)
    check_smooths("tikhonov -m 2", ["-m", "2", "-l", "auto"], uni, clean)
    check_smooths("butterworth -m 3", ["-m", "3", "-f", "auto"], uni, clean)

    print("non-uniform grid — grid-tolerant methods:")
    check_smooths("polyfit  -m 0", ["-m", "0", "-n", "7", "-p", "3"], non, clean_nu)
    check_smooths("tikhonov -m 2", ["-m", "2", "-l", "auto"], non, clean_nu)

    print("non-uniform grid — preconditions must be rejected LOUDLY (exit 1):")
    for name, args in [("savgol", ["-m", "1", "-n", "7"]), ("butterworth", ["-m", "3"])]:
        p = run(args + [non], expect=1)
        assert "grid" in p.stderr.lower(), f"{name}: rejection did not mention the grid"
        print(f"  [OK] {name} rejected: {p.stderr.strip().splitlines()[-1][:60]}")

    print("stdout table integrity — every diagnostic line must be `#`-prefixed:")
    # Regression guard for the multi-line warning_msg leak fixed in v5.11.47:
    # continuation lines used to land on stdout bare and corrupt the table.
    for label, args in [
        ("CV>0.2  (noticeable, 2-line msg)", ["-m", "0", "-n", "7", non]),
        ("CV>0.2  (-g verbose)", ["-g", non]),
    ]:
        leaks = leaked_lines(run(args).stdout)
        assert not leaks, f"{label}: {len(leaks)} unprefixed line(s), e.g. {leaks[0]!r}"
        print(f"  [OK] {label}: no unprefixed lines")

    print("derivatives (-d) add a third column:")
    rows = columns(run(["-m", "2", "-l", "auto", "-d", uni]).stdout)
    assert all(len(r) == 3 for r in rows), "-d did not produce 3 columns"
    # d/dx sin(x) = cos(x); a correct derivative tracks it
    derr = rmse([r[2] for r in rows], [math.cos(r[0]) for r in rows])
    assert derr < 0.2, f"derivative far from cos(x): RMSE {derr:.3f}"
    print(f"  [OK] dy/dx vs cos(x) RMSE {derr:.5f}")

    print("grid analysis (-g):")
    g = run(["-g", uni]).stdout
    assert "UNIFORM" in g and "CV = 0.000" in g, g
    print("  [OK] -g reports CV = 0.000, UNIFORM")

    print("timestamp mode (-T), both separators:")
    for sep in (" ", "T"):
        ts = tmp / f"ts{sep.strip() or 'sp'}.dat"
        ts.write_text(
            "".join(
                f"2025-09-25{sep}14:06:{i:02d}.000 {signal(i * 0.1):.5f}\n"
                for i in range(20)
            )
        )
        out = run(["-T", "-m", "2", "-l", "0.1", "-d", str(ts)]).stdout
        rows = [l for l in out.splitlines() if l.startswith("2025-")]
        assert len(rows) == 20, f"sep {sep!r}: {len(rows)} rows"
        # The timestamp is echoed verbatim, so a space-separated one costs TWO
        # whitespace fields: 'DATE TIME y dy/dt' vs 'DATEtTIME y dy/dt'.
        want = 4 if sep == " " else 3
        assert len(rows[0].split()) == want, f"sep {sep!r}: {rows[0]!r}"
        assert rows[0].startswith(f"2025-09-25{sep}14:06:00.000"), rows[0]
        print(f"  [OK] separator {sep!r}: 20 rows, {want} fields, verbatim timestamps")

    print("column selection (-k):")
    wide = tmp / "wide.dat"
    src = columns(Path(uni).read_text())
    wide.write_text("".join(f"{x} 999 {y} 999\n" for x, y in src))
    a = columns(run(["-k", "1:3", "-m", "0", str(wide)]).stdout)
    b = columns(run(["-m", "0", uni]).stdout)
    assert a == b, "-k 1:3 on a padded file != plain 2-column file"
    print("  [OK] -k 1:3 selects x=col1, y=col3")

    print("stdin filter:")
    p = subprocess.run(
        [str(SMOOTH), "-m", "0", "-n", "5"],
        input=Path(uni).read_text(), capture_output=True, text=True,
    )
    assert p.returncode == 0 and len(columns(p.stdout)) == N, p.stderr
    print("  [OK] reads stdin when no file argument is given")

    print("error paths:")
    for label, args in [
        ("missing file", ["-m", "0", str(tmp / "nope.dat")]),
        ("bad method", ["-m", "9", uni]),
        ("negative lambda", ["-m", "2", "-l", "-5", uni]),
        ("too few points", ["-m", "0", str(tmp / "tiny.dat")]),
    ]:
        (tmp / "tiny.dat").write_text("0 1\n1 2\n")
        p = run(args, expect=1)
        assert p.stderr.strip(), f"{label}: exited 1 with an empty stderr"
        print(f"  [OK] {label:<16} exit 1: {p.stderr.strip().splitlines()[0][:50]}")

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
