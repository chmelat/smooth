---
name: smooth-dev-tasks
description: Step-by-step recipes for developing on the `smooth` codebase — adding a new smoothing method, modifying shared grid analysis, the result-struct memory-management pattern, the diagnostic output convention, and adding a Unity test. Load when implementing any of these.
---

# smooth development tasks

Recipes for common changes to the `smooth` codebase. Design rationale, LAPACK
choices, grid-uniformity policy, and hard rules live in the root `CLAUDE.md`;
these are the how-to steps.

## Adding a test

1. Add the test function in the relevant `tests/test_<module>.c`. Use
   `TEST_ASSERT_DOUBLE_WITHIN()` for floating-point comparisons.
2. Declare the function and `RUN_TEST(...)` it in `tests/test_main.c`.
3. If the test exercises a new module, append the source to `TEST_SRC` and the
   object to `TEST_MODULES` in the `Makefile`.
4. Run `make test` and `make test-valgrind` before committing.

Test best practices:

- One logical assertion per test.
- Cover edge cases: `n=1`, `n=2`, NULL inputs, empty arrays.
- Always free in the test body — the leak budget is fixed.

## Adding a new smoothing method

1. Create `newmethod.c` / `newmethod.h` following the existing pattern.
2. Define `NewmethodResult` with at minimum `y_smooth`, `y_deriv`, `n`.
3. Implement `newmethod_smooth(...)` accepting `const GridAnalysis*` —
   policy decisions about grid uniformity belong here.
4. Implement `free_newmethod_result()`.
5. In `smooth.c`: add `METHOD_NEWMETHOD` constant, dispatch case, output path.
6. Update `Makefile` (`SRC`, `HEAD`).
7. Bump `revision.h` (version + date + history line).
8. Add `tests/test_newmethod.c` and wire it into `tests/test_main.c` and the
   `Makefile`.

## Modifying grid analysis

`GridAnalysis` is shared state — every change ripples to all methods.

1. Update the struct in `grid_analysis.h`.
2. Populate the new field in `analyze_grid()`.
3. Update every method that should react to it.
4. Update `tests/test_grid_analysis.c`.

## Memory management

```c
NewmethodResult* newmethod_smooth(...) {
    NewmethodResult *r = malloc(sizeof(*r));
    r->y_smooth = malloc(n * sizeof(double));
    r->y_deriv  = malloc(n * sizeof(double));
    return r;
}

void free_newmethod_result(NewmethodResult *r) {
    if (r) {
        free(r->y_smooth);
        free(r->y_deriv);
        free(r);
    }
}
```

Verify with `make test-valgrind` before committing.

## Diagnostic output convention (butterworth, applies to new modules)

- `stdout` as `# ...` — info that should be preserved with the saved data
  (selected parameters, grid CV, numerical-quality warnings).
- `stderr` as `Warning: ...` — runtime/operational concerns (memory usage)
  that do not belong in the data file.
- `stderr` as `ERROR: ...` — hard failures; function returns NULL/non-zero.
