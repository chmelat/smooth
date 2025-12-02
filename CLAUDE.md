# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`smooth` is a scientific data smoothing program implementing four mathematical methods for noise reduction and derivative computation in experimental data. The codebase is ~3,600 lines of modular C code with LAPACK dependencies.

**Current version:** 5.7.1 (2025-11-23)

## Documentation Guidelines

### README.md Character Restrictions

**IMPORTANT:** The README.md file uses Unicode mathematical symbols but **restricts non-mathematical special characters** for PDF conversion compatibility.

**Allowed characters:**
- **Mathematical symbols:** All Greek letters (λ, ω, π, Σ, ε, μ, θ, etc.)
- **Mathematical operators:** ², ³, ⁴, ₀, ₁, ₂, ∫, ∂, √, ±, ×, ≈, ≥, ≤, ·
- **Standard ASCII:** Letters, numbers, punctuation

**Prohibited characters (use ASCII alternatives):**
- **Box-drawing:** Use `|-`, `-`, `|`, `+-`, `=` instead of ├ ─ │ └ ━
- **Checkmarks/warnings:** Use `[OK]`, `[X]`, `[WARNING]` instead of ✓ ✗ ⚠️
- **Arrows:** Use `->`, `==>` instead of → ⟹
- **Ellipsis:** Use `...` instead of ⋮ ⋱

**Reason:** DejaVu fonts used for PDF generation support mathematical Unicode but have limited support for decorative box-drawing characters. This restriction ensures README.md converts cleanly to PDF while preserving all mathematical notation.

**When editing README.md:** Do not introduce box-drawing characters, special arrows, or decorative Unicode. Use ASCII art alternatives.

## Build System

```bash
# Build production version
make

# Build with debug symbols
make debug

# Run unit tests
make test

# Run tests with Valgrind (memory leak detection)
make test-valgrind

# Clean build artifacts
make clean

# Install to user's ~/bin
make install-user

# Install system-wide (requires root)
make install

# Show all available targets
make help
```

**Important Build Notes:**
- Uses `clang` by default (can override with `CC=gcc`)
- Requires LAPACK and BLAS libraries (`-llapack -lblas`)
- Default library path: `~/lib` (override with `LIBDIR=/path/to/libs`)
- Production builds use `-O2`, debug builds use `-g -O0`

## Running the Program

```bash
# Basic usage
./smooth -m <method> [options] [file|-]

# From file
./smooth -m 2 -l auto data.txt

# From stdin (Unix filter)
cat data.txt | ./smooth -m 1 -n 7 -p 2

# Grid analysis only
./smooth -g data.txt
```

**Key Options:**
- `-m {0|1|2|3}` - Method: polyfit|savgol|tikhonov|butterworth
- `-n N` - Window size (polyfit, savgol)
- `-p P` - Polynomial degree (polyfit, savgol, max 12)
- `-l λ` - Regularization parameter (tikhonov), use `-l auto` for GCV
- `-f fc` - Cutoff frequency (butterworth, 0 < fc < 0.5), use `-f auto` for default
- `-d` - Include first derivatives in output
- `-g` - Show detailed grid uniformity analysis

## Code Architecture

### Modular Structure

The codebase follows a clean modular design with each smoothing method in its own compilation unit:

```
smooth.c              # Main program, CLI parsing, I/O
├─ polyfit.c/h        # Polynomial fitting (local least squares)
├─ savgol.c/h         # Savitzky-Golay filter (pre-computed convolution)
├─ tikhonov.c/h       # Tikhonov regularization (global variational)
├─ butterworth.c/h    # Butterworth filter (frequency-domain)
├─ grid_analysis.c/h  # Grid uniformity analysis (shared utility)
└─ decomment.c/h      # Input comment removal
```

**Design Principles:**
1. **One analysis, shared results:** Grid analysis happens once at startup, results passed to all methods
2. **Result structures:** Each method returns a struct with smoothed values, derivatives, and diagnostics
3. **Memory ownership:** Caller must free results using `free_*_result()` functions
4. **Grid-aware methods:** Methods receive `GridAnalysis*` to make informed decisions

### Core Data Structures

```c
// Grid analysis (shared across all methods)
typedef struct {
    double h_min, h_max, h_avg, h_std;
    double ratio_max_min, cv;
    int is_uniform, n_clusters;
    double uniformity_score;
    double *spacings;  // Optional, must be freed
    // ... warning fields
} GridAnalysis;

// Method result structures (similar pattern)
typedef struct {
    double *y_smooth;     // Smoothed values
    double *y_deriv;      // First derivatives
    int n;                // Number of points
    // ... method-specific fields
} MethodResult;
```

### LAPACK Usage

All methods use LAPACK for linear algebra:

- **polyfit:** `dposv` - Symmetric positive definite solver (normal equations)
- **savgol:** `dposv` - Symmetric positive definite solver (coefficient computation)
- **tikhonov:** `dpbsv` - Banded symmetric positive definite solver (tridiagonal system)
- **butterworth:** `dgesv` - General linear solver (initial conditions via LU decomposition)

**Key:** Methods rely on the symmetric positive definite structure of their formulations for numerical stability.

### Grid Uniformity Philosophy

Grid analysis is central to method selection:

**Coefficient of Variation (CV):**
```
CV = std_dev(spacing) / avg(spacing)

CV < 0.01:  Perfectly uniform
CV < 0.05:  Nearly uniform (Savgol warns but works)
CV < 0.15:  Moderately non-uniform (Tikhonov uses average coefficient method)
CV ≥ 0.15:  Highly non-uniform (Tikhonov uses local spacing method)
CV > 0.05:  Savgol REJECTS with detailed error message
```

**When modifying grid-dependent code:**
- Savitzky-Golay assumes uniform grids (mathematical requirement, not implementation choice)
- Tikhonov automatically switches discretization based on grid uniformity
- Butterworth assumes nearly-uniform sampling (frequency analysis requirement)
- Polyfit is most tolerant of grid non-uniformity

## Testing Framework

Uses **Unity testing framework** (included in `tests/` directory).

**Current test coverage:**
- `test_grid_analysis.c` - Grid analysis module (7 tests)
- `test_polyfit.c` - Polynomial fitting module
- `test_main.c` - Test runner

**To add new tests:**

1. Create `tests/test_module.c`:
```c
#include "unity.h"
#include "../module.h"

void setUp(void) { /* runs before each test */ }
void tearDown(void) { /* runs after each test */ }

void test_module_feature(void) {
    // ARRANGE
    double x[] = {0.0, 1.0, 2.0};

    // ACT
    Result *result = module_function(x, 3);

    // ASSERT
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, expected, result->value);

    // CLEANUP
    free_module_result(result);
}
```

2. Declare test in `tests/test_main.c`:
```c
void test_module_feature(void);  // Add declaration
RUN_TEST(test_module_feature);   // Add to main()
```

3. Update `Makefile`:
```makefile
TEST_SRC = ... tests/test_module.c
TEST_MODULES = ... module.o
```

**Testing Best Practices:**
- Always use `TEST_ASSERT_DOUBLE_WITHIN()` for floating-point comparisons
- Test edge cases: n=1, n=2, NULL pointers, empty arrays
- Always free allocated memory (check with `make test-valgrind`)
- One assertion per logical concept
- Use descriptive test names

## Mathematical Methods Summary

### Method Selection Guide

**Grid is uniform (CV < 0.05):**
- **Savgol:** Best for polynomial signals, need derivatives, computational efficiency critical
- **Butterworth:** Best for frequency-domain interpretation, zero phase distortion needed
- **Tikhonov + auto:** Universal choice, automatic parameter selection via GCV

**Grid is non-uniform (CV ≥ 0.05):**
- **Tikhonov + auto:** Primary recommendation (handles any grid correctly)
- **Polyfit:** Local adaptability, good for variable curvature

**Not sure which method:**
- Use `./smooth -g data.txt` first to analyze grid
- Default to **Tikhonov with `-l auto`** (safest, most automatic)

### Implementation Details

**Polynomial Fitting (polyfit.c):**
- Local polynomial least squares in sliding window
- Solves normal equations at each point: `(X^T X)a = X^T y`
- Asymmetric windows at boundaries with polynomial extrapolation
- O(n·p³) complexity

**Savitzky-Golay (savgol.c):**
- Pre-computes universal convolution coefficients
- Method of undetermined coefficients (moment conditions)
- **Requires uniform grid** - enforced by CV check
- O(p³) + O(n·w) complexity
- Key difference from polyfit: same coefficients applied everywhere (translation invariance)

**Tikhonov (tikhonov.c):**
- Global variational minimization: `J[u] = ||y-u||² + λ||D²u||²`
- **Hybrid discretization (v5.4+):**
  - CV < 0.15: Average coefficient method (harmonic mean)
  - CV ≥ 0.15: Local spacing method (Taylor expansion)
- Natural boundary conditions (D²u = 0 at ends)
- Tridiagonal banded solver for O(n) efficiency
- GCV for automatic λ selection (eigenvalue trace approximation)
- O(n) complexity
- **Boundary artifacts:** On non-uniform grids (CV > 0.15) with small λ, last 2-3 points may show significant oscillations (up to 60% deviation). This is documented behavior, not a bug. See README.md for mitigation strategies.

**Butterworth (butterworth.c):**
- 4th-order digital low-pass filter
- Bilinear transform: s-domain poles → z-domain
- Biquad cascade (two 2nd-order sections)
- **Filtfilt:** Forward-backward filtering for zero phase
- **lfilter_zi:** Scipy-compatible initial conditions via companion matrix + LAPACK dgesv
- Odd reflection padding for edge handling
- O(n) complexity

## Common Development Tasks

### Adding a New Smoothing Method

1. Create `newmethod.c` and `newmethod.h` following existing pattern
2. Define result structure with `y_smooth`, `y_deriv`, `n` fields
3. Implement `newmethod_smooth()` accepting `GridAnalysis*`
4. Implement `free_newmethod_result()`
5. Add to `smooth.c`:
   - Define `METHOD_NEWMETHOD` constant
   - Add case to method selection switch
   - Add call to smoothing function
   - Add output formatting
6. Update `Makefile`: add to `SRC` and `HEAD` variables
7. Update `revision.h` version number
8. Create unit tests in `tests/test_newmethod.c`

### Modifying Grid Analysis

Grid analysis is centralized and shared. Changes affect all methods:

1. **Location:** `grid_analysis.c`, called once from `smooth.c` after data load
2. **Results passed to:** All method functions via `GridAnalysis*` parameter
3. **Update checklist when modifying:**
   - Update `GridAnalysis` struct in `grid_analysis.h`
   - Update `analyze_grid()` computation
   - Update all methods that use the new field
   - Update tests in `test_grid_analysis.c`

### Memory Management Rules

**Allocation patterns:**
```c
// Methods allocate and return structures
MethodResult* method_smooth(...) {
    MethodResult *result = malloc(sizeof(MethodResult));
    result->y_smooth = malloc(n * sizeof(double));
    result->y_deriv = malloc(n * sizeof(double));
    return result;
}

// Caller must free
void free_method_result(MethodResult *result) {
    if (result) {
        free(result->y_smooth);
        free(result->y_deriv);
        free(result);
    }
}
```

**Always verify with:** `make test-valgrind` before committing

## Version History Context

**v5.7.1 (current):** Added polyfit unit tests, small bug fixes
**v5.6:** First unity tests added
**v5.5:** Butterworth filter added, Unix filter support, centralized grid analysis
**v5.4:** Tikhonov hybrid discretization (auto-switch at CV=0.15), GCV improvements
**v5.3:** Savitzky-Golay grid uniformity enforcement
**v5.2:** Grid analysis module with `-g` flag
**v5.1:** Optional derivative output with `-d` flag
**v5.0:** Complete modularization

## File I/O Expectations

**Input format:**
- Two-column ASCII: `x y`
- Comments (lines starting with `#`) are stripped by `decomment.c`
- Supports stdin (`-` or pipe) for Unix filter usage
- Data must have strictly monotonic increasing x-values

**Output format:**
- Header comments with method info and diagnostics
- Without `-d`: `x y_smooth`
- With `-d`: `x y_smooth y_deriv`
- Scientific notation: `%.5E` format

## Dependencies

**Required:**
- LAPACK (linear algebra)
- BLAS (basic linear algebra)
- Standard math library (-lm)
- C99 or later (uses `complex.h` in butterworth.c)

**Development:**
- Unity testing framework (included in tests/)
- Valgrind (optional, for memory checking)

## Notes for Future Development

**Suggested improvements:**
- Automatic cutoff frequency selection for Butterworth (currently returns fixed 0.1)
- Integration tests (full pipeline)
- Performance benchmarks framework
- Coverage reporting (gcov/lcov)

**Do not:**
- Break the modular architecture (keep methods independent)
- Introduce cross-dependencies between method modules
- Change grid analysis behavior without updating all methods
- Commit without running `make test` and checking for memory leaks
