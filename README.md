# Methods for Smoothing Experimental Data in the smooth Program

**Technical Documentation**  
Version 5.3 | October 2025

---

## Contents

1. [Introduction](#introduction)
2. [Polynomial Fitting (POLYFIT)](#polynomial-fitting-polyfit)
3. [Savitzky-Golay Filter (SAVGOL)](#savitzky-golay-filter-savgol)
4. [Tikhonov Regularization (TIKHONOV)](#tikhonov-regularization-tikhonov)
5. [Method Comparison](#method-comparison)
6. [Practical Recommendations](#practical-recommendations)
7. [Usage Examples](#usage-examples)
8. [Grid Analysis Module](#grid-analysis-module)
9. [Compilation and Installation](#compilation-and-installation)

---

## Introduction

The `smooth` program implements three sophisticated methods for smoothing experimental data with the capability of simultaneous derivative computation. Each method has specific properties, advantages, and areas of application.

### General Smoothing Problem

Experimental data often contains random noise:
```
y_obs(x_i) = y_true(x_i) + ε_i
```

The goal of smoothing is to estimate `y_true` while suppressing `ε_i` and preserving physically relevant signal properties.

### Program Structure

```bash
./smooth [options] data_file
```

**Basic Parameters:**
- `-m {0|1|2}` - method selection (polyfit|savgol|tikhonov)
- `-n N` - smoothing window size (polyfit, savgol)
- `-p P` - polynomial degree (polyfit, savgol, max 12)
- `-l λ` - regularization parameter (tikhonov)
- `-l auto` - automatic λ selection using GCV (tikhonov)
- `-d` - display first derivative in output (optional)

**Note on polynomial degree:** Degrees > 6 may generate numerical stability warnings.

**Note on derivatives:** From version 5.1, first derivative output is optional. Without the `-d` switch, the program outputs only smoothed values. With the `-d` switch, it outputs both smoothed values and first derivatives.

---

## What's New in Version 5.3

### Major Improvements

**Tikhonov Regularization:**
- **Simplified API:** Removed `adaptive_weights` option - method now automatically uses correct discretization
- **Mathematically correct D² operator:** Proper implementation for non-uniform grids without user intervention
- **Better accuracy:** Improved results on non-uniform grids by default
- **Cleaner code:** ~150 lines of code removed while improving correctness

**Savitzky-Golay Filter:**
- **Grid uniformity checking:** Method now verifies grid uniformity before processing
- **Automatic rejection:** Non-uniform grids (CV > 0.05) are rejected with helpful error message
- **Recommendations:** When rejected, program suggests appropriate alternative methods
- **Mathematical correctness:** Prevents incorrect application of SG filter to inappropriate data

### API Structure

All smoothing methods have consistent APIs:

```c
// Polyfit
PolyfitResult* polyfit_smooth(double *x, double *y, int n, 
                              int window_size, int poly_degree);
void free_polyfit_result(PolyfitResult *result);

// Savitzky-Golay
SavgolResult* savgol_smooth(double *x, double *y, int n, 
                            int window_size, int poly_degree);
void free_savgol_result(SavgolResult *result);

// Tikhonov - SIMPLIFIED in v5.3
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda);
void free_tikhonov_result(TikhonovResult *result);
```

**Note:** The `tikhonov_smooth_adaptive()` function has been removed in v5.3. The standard `tikhonov_smooth()` function now automatically uses correct discretization for both uniform and non-uniform grids.

---

## Polynomial Fitting (POLYFIT)

### Mathematical Foundations

The POLYFIT method uses local polynomial fitting with least squares method in a sliding window.

**Problem:** For each point `x_i`, we fit a polynomial of degree `p` to the surrounding `n` points:

```
P(x) = a_0 + a_1(x-x_i) + a_2(x-x_i)² + ... + a_p(x-x_i)^p
```

**Optimization criterion:**
```
min Σ[y_j - P(x_j)]²    for j ∈ [i-n/2, i+n/2]
```

### Construction of Normal Equations

For polynomial coefficients, we solve a system of linear equations:

```
[Σ(x-x_i)⁰    Σ(x-x_i)¹    ...  Σ(x-x_i)^p  ] [a_0]   [Σy(x-x_i)⁰]
[Σ(x-x_i)¹    Σ(x-x_i)²    ...  Σ(x-x_i)^{p+1}] [a_1] = [Σy(x-x_i)¹]
[    ⋮            ⋮        ⋱        ⋮        ] [ ⋮ ]   [    ⋮     ]
[Σ(x-x_i)^p Σ(x-x_i)^{p+1} ... Σ(x-x_i)^{2p}] [a_p]   [Σy(x-x_i)^p]
```

where summation is over points in the window around `x_i`.

### Derivative Computation

Derivatives are computed analytically from polynomial coefficients:

```
f(x_i) = a_0
f'(x_i) = a_1
f''(x_i) = 2a_2
```

### Edge Handling

At edges, asymmetric windows are used with extrapolation of the fitted polynomial:

```c
// For point x_k < x_{n/2}
f(x_k) = Σ_{m=0}^p a_m * (x_k - x_{n/2})^m
```

### Efficient Implementation

The program uses LAPACK routine `dposv` for solving symmetric positive definite systems at each point:

```c
// System solution for polynomial coefficients
dposv_(&uplo, &matrix_size, &nrhs, C, &matrix_size, B, &matrix_size, &info);
```

The normal equations matrix is symmetric and positive definite, making `dposv` optimal for this application.

### Modularized Implementation

```c
// polyfit.h
typedef struct {
    double *y_smooth;     // Smoothed values
    double *y_deriv;      // First derivatives
    int n;                // Number of points
    int poly_degree;      // Polynomial degree
    int window_size;      // Window size
} PolyfitResult;
```

### Characteristics

**Advantages:**
- Excellent local approximation
- Analytical computation of derivatives of any order
- Adaptable to changes in curvature
- Good preservation of local extrema
- Works with moderately non-uniform grids

**Disadvantages:**
- Sensitive to outliers
- Boundary effects at edges
- Possible Runge oscillations for high polynomial degrees (p > 6)
- Numerical instability warnings for degrees > 6

---

## Savitzky-Golay Filter (SAVGOL)

### Theoretical Foundations

The Savitzky-Golay filter is an optimal linear filter for smoothing and derivatives based on local polynomial regression. The key innovation is pre-computation of convolution coefficients.

**Fundamental principle:** For given parameters (window size, polynomial degree, derivative order), there exist universal coefficients `c_k` such that:

```
f^(d)(x_i) = Σ_{k=-n_L}^{n_R} c_k · y_{i+k}
```

### Key Difference from POLYFIT Method

While both SAVGOL and POLYFIT use polynomial approximation, they differ fundamentally in their computational approach:

**POLYFIT approach:**
- For each data point, fits a new polynomial to the surrounding window
- Solves the least squares problem individually for each point
- Coefficients of the polynomial change with each window position
- Computationally intensive: O(n·p³)

**SAVGOL approach (Method of Undetermined Coefficients):**
- Recognizes that for equidistant grids, the filter coefficients are translation-invariant
- Uses the **method of undetermined coefficients** to pre-compute universal weights
- These weights depend only on the window geometry, not on the actual data values
- Applies the same weights as a linear convolution across all data points
- Computationally efficient: O(p³) once, then O(n·w) for application

### **CRITICAL: Grid Uniformity Requirement**

**Version 5.3 Important Change:** The Savitzky-Golay method now enforces grid uniformity checking.

The mathematical foundation of SG filter assumes **uniformly spaced data points**. The method is based on fitting polynomials in normalized coordinate space where points are at integer positions: {..., -2, -1, 0, 1, 2, ...}.

**Uniformity Check:**
```
CV = std_dev(spacing) / avg(spacing)

If CV > 0.05:  REJECT - Grid too non-uniform for SG
If CV > 0.01:  WARNING - Nearly uniform, proceed with caution
If CV ≤ 0.01:  OK - Grid sufficiently uniform
```

**What happens when grid is rejected:**
```
========================================
ERROR: Savitzky-Golay method not suitable for non-uniform grid!
========================================
Grid analysis:
  Coefficient of variation (CV) = 0.2341
  Threshold for uniformity = 0.0500
  
RECOMMENDED ALTERNATIVES:
  1. Use Tikhonov method: -m 2 -l auto
     (Works correctly with non-uniform grids)
  2. Use Polyfit method: -m 0 -n 5 -p 2
     (Local fitting, less sensitive to spacing)
  3. Resample your data to uniform grid before smoothing
```

### The Method of Undetermined Coefficients

The Savitzky-Golay method seeks a linear combination of data points:

```
ŷ₀ = c₋ₙₗ·y₋ₙₗ + ... + c₀·y₀ + ... + cₙᵣ·yₙᵣ
```

where the coefficients `c_k` are "undetermined" and must satisfy the condition that the filter exactly reproduces polynomials up to degree `p`.

**The key insight:** For a given window configuration and polynomial degree, these coefficients can be determined once and applied universally - but only on uniform grids!

### Coefficient Derivation

Coefficients are derived from the condition that the filter must exactly reproduce polynomials up to degree `p`.

**Moment conditions:**
```
Σ_{j=-n_L}^{n_R} c_j · j^m = δ_{m,d} · d!    for m = 0,1,...,p
```

where:
- `δ_{m,d}` is the Kronecker delta
- `d` is the derivative order
- `d!` is factorial

This leads to a system of linear equations where the unknowns are the filter coefficients `c_j`.

### Matrix Formulation

We solve a system of linear equations:

```
[1   -n_L    (-n_L)²   ...  (-n_L)^p  ] [c_{-n_L}]   [δ_{0,d}·0!]
[1  -n_L+1  (-n_L+1)²  ... (-n_L+1)^p ] [c_{-n_L+1}] = [δ_{1,d}·1!]
[⋮     ⋮        ⋮      ⋱       ⋮      ] [    ⋮    ]   [    ⋮     ]
[1    n_R     n_R²     ...    n_R^p   ] [  c_{n_R} ]   [δ_{p,d}·p!]
```

### Computational Efficiency

The brilliance of the Savitzky-Golay approach becomes apparent when processing large datasets:

**Example for 10,000 data points, window size 21, polynomial degree 4:**
- **POLYFIT:** Must solve 10,000 separate 5×5 linear systems
- **SAVGOL:** Solves only ONE 5×5 system, then performs 10,000 simple weighted sums

This difference explains why SAVGOL is preferred for real-time signal processing and large datasets, while maintaining the same mathematical accuracy as POLYFIT **for uniform grids**.

### Efficient Implementation

The program uses LAPACK routine `dposv` for solving the symmetric positive definite system when computing filter coefficients:

```c
// Solve linear system for Savitzky-Golay coefficients
dposv_(&uplo, &matrix_size, &nrhs, A, &matrix_size, B, &matrix_size, &info);
```

### Modularized Implementation

```c
// savgol.h
typedef struct {
    double *y_smooth;     // Smoothed values
    double *y_deriv;      // First derivatives
    int n;                // Number of points
    int poly_degree;      // Polynomial degree
    int window_size;      // Window size
} SavgolResult;

// Coefficient computation
void savgol_coefficients(int nl, int nr, int poly_degree, 
                        int deriv_order, double *c);
```

### Optimal Properties

The Savitzky-Golay filter minimizes approximation error in the least squares sense and maximizes signal-to-noise ratio for polynomial signals **on uniform grids**.

### Characteristics

**Advantages:**
- Optimal for polynomial signals on uniform grids
- Excellent preservation of moments and peak areas
- Efficient implementation (convolution)
- Minimal phase distortion
- Simultaneous computation of functions and derivatives

**Disadvantages:**
- **Requires uniform grid** - automatically rejected if CV > 0.05
- Fixed coefficients for entire window
- May introduce oscillations at sharp edges
- Limited adaptability
- Numerical warnings for degrees > 6

---

## Tikhonov Regularization (TIKHONOV)

### Theoretical Foundation

Tikhonov regularization solves the ill-posed inverse smoothing problem using a variational approach. We seek a function minimizing the functional:

**Continuous formulation:**
```
J[u] = ∫(y(x) - u(x))² dx + λ∫(u''(x))² dx
```

**Discrete formulation:**
```
J[u] = ||y - u||² + λ||D²u||²
```

where:
- `||y - u||²` = Σ(y_i - u_i)² is the data fidelity term
- `||D²u||²` = Σ(D²u_i)² is the regularization term
- `λ` is the regularization parameter controlling smoothness
- `D²` is the discrete second derivative operator

### Second Derivative Discretization

**Version 5.3 Improvement:** Correct implementation for both uniform and non-uniform grids.

**For uniform grid with step `h`:**
```
D²u_i ≈ (u_{i-1} - 2u_i + u_{i+1})/h²
```

**For non-uniform grid (NEW in v5.3):**
For point `i` with left spacing `h₁ = x[i] - x[i-1]` and right spacing `h₂ = x[i+1] - x[i]`:

```
D²u_i ≈ 2/(h₁+h₂) · [u_{i-1}/h₁ - u_i·(1/h₁+1/h₂) + u_{i+1}/h₂]
```

This leads to a **symmetric tridiagonal matrix** with variable coefficients:

```
A[i,i-1] = 2λ/(h₁·(h₁+h₂))
A[i,i]   = 1 + 2λ/(h₁+h₂)·(1/h₁+1/h₂)
A[i,i+1] = 2λ/(h₂·(h₁+h₂))
```

**Key improvement:** The method automatically uses the correct discretization. No user intervention or parameter selection needed - it "just works" for any grid spacing.

### Variational Approach

The minimum of functional J[u] satisfies the Euler-Lagrange equation:

```
∂J/∂u_i = 0  ⟹  -2(y_i - u_i) + 2λ(D^T D u)_i = 0
```

which leads to the linear system:
```
(I + λD^T D)u = y
```

### Matrix Representation

Matrix `A = I + λD^T D` is:
- Symmetric
- Positive definite
- Tridiagonal (banded with bandwidth 1)

This structure allows efficient solution using LAPACK's banded solver `dpbsv`.

### Boundary Conditions

The program implements natural boundary conditions:
```
u''(0) = u''(L) = 0
```

This is realized by matrix modification at edges using first-order differences.

### Automatic Grid Handling (NEW in v5.3)

**Previous versions (≤5.2):** Required user to choose between:
- Average coefficient method (default)
- Adaptive weights method (`-a` flag)

**Version 5.3:** Automatic correct discretization!
- No user choice needed
- Mathematically correct for uniform grids
- Mathematically correct for non-uniform grids
- Same code path for both cases
- Simplified API

### Generalized Cross Validation (GCV)

For automatic λ selection, we minimize the GCV criterion:

```
GCV(λ) = n·RSS(λ) / (n - tr(H_λ))²
```

where:
- `RSS(λ) = ||y - u_λ||²` is the residual sum of squares
- `H_λ = (I + λD^T D)^{-1}` is the influence matrix
- `tr(H_λ)` is the matrix trace

**Trace estimation using eigenvalues (corrected in v5.3):**
```
tr(H_λ) ≈ Σ_{k=1}^n 1/(1 + λμ_k)
```

For natural boundary conditions:
```
θ_k = πk/n  (corrected from πk/(n+1))
μ_k = 4·sin²(θ_k/2) / h²
```

### Efficient Implementation

The program uses LAPACK routine `dpbsv` for solving symmetric positive definite banded systems:

```c
// Banded matrix storage (LAPACK format)
AB[0,j] = superdiagonal elements
AB[1,j] = diagonal elements

// System solution
dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
```

### Simplified Implementation (v5.3)

```c
typedef struct {
    double *y_smooth;            // Smoothed values
    double *y_deriv;             // First derivatives
    double lambda;               // Used parameter
    int n;                       // Number of points
    double data_term;            // ||y - u||²
    double regularization_term;  // λ||D²u||²
    double total_functional;     // J[u]
} TikhonovResult;

// Single simple function - no more adaptive_weights parameter!
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda);
```

### Characteristics

**Advantages:**
- Global optimization with theoretical foundation
- Flexible balance between data fidelity and smoothness
- Robust to outliers
- Efficient for large datasets (O(n) memory and time)
- Automatic parameter selection via GCV
- **NEW v5.3:** Excellent for non-uniform grids automatically
- **NEW v5.3:** Simplified usage - no complex parameter choices

**Disadvantages:**
- Single global parameter λ
- May suppress local details
- GCV may fail for some data types
- Requires LAPACK library

---

## Method Comparison

### Computational Complexity

| Method | Time | Memory | Scalability |
|--------|------|--------|-------------|
| POLYFIT | O(n·p³) | O(p²) | Good for small p |
| SAVGOL | O(p³) + O(n·w) | O(w) | Excellent for large n |
| TIKHONOV | O(n) | O(n) | Excellent |

*Note: w = window size, p = polynomial degree (≤12), n = number of data points.*

### Smoothing Quality

| Property | POLYFIT | SAVGOL | TIKHONOV |
|----------|---------|--------|----------|
| Local adaptability | ***** | **** | ** |
| Extreme preservation | **** | ***** | *** |
| Noise robustness | *** | **** | ***** |
| Derivative quality | ***** | ***** | *** |
| Boundary behavior | ** | *** | **** |
| Non-uniform grids | *** | ❌ | ***** |
| Ease of use (v5.3) | **** | **** | ***** |

**Key:** ❌ = Not suitable (automatically rejected)

### Grid Type Compatibility

| Grid Type | POLYFIT | SAVGOL | TIKHONOV |
|-----------|---------|--------|----------|
| Perfectly uniform (CV < 0.01) | ✓ | ✓ | ✓ |
| Nearly uniform (CV < 0.05) | ✓ | ⚠️ | ✓ |
| Moderately non-uniform (CV < 0.2) | ✓ | ❌ | ✓ |
| Highly non-uniform (CV > 0.2) | ⚠️ | ❌ | ✓ |

**Legend:**
- ✓ = Recommended
- ⚠️ = Usable with caution
- ❌ = Rejected or not recommended

---

## Practical Recommendations

### Method Selection by Data Type

#### POLYFIT - when:
- Data has variable curvature
- You need to preserve local details
- You have moderately noisy data
- You want highest quality derivatives
- Grid has moderate spacing variations

#### SAVGOL - when:
- **Grid is uniform (CV < 0.05)** - automatically checked!
- You want mathematically optimal linear smoothing for polynomial signals
- Data contains periodic or oscillatory components that need preservation
- You need excellent peak shape preservation (areas, moments)
- You want minimal phase distortion in the smoothed signal
- You're processing time series or spectroscopic data on uniform grids
- You need simultaneous high-quality function and derivative estimation

#### TIKHONOV - when:
- **Grid is non-uniform** - now works perfectly automatically!
- Data is very noisy
- You need global consistency
- You want automatic parameter selection
- You prefer global optimization approaches over local fitting
- You want the simplest interface (just λ parameter)

### Parameter Selection

#### Window size (n) for POLYFIT/SAVGOL:
```
n = 2*k + 1    (odd number)

Recommendations:
- Low noise: n = 5-9
- Medium noise: n = 9-15  
- High noise: n = 15-25
```

#### Polynomial degree (p):
```
- Linear trends: p = 1-2
- Smooth curves: p = 2-3
- Complex signals: p = 3-4
- Advanced applications: p = 5-8
- Maximum: p ≤ 12
- Recommended maximum: p < n/2

Note: Degrees > 6 may cause numerical instability warnings.
```

#### Lambda (λ) for TIKHONOV:
```
- Auto selection: -l auto  (recommended!)
- Manual start: λ = 0.1
- More smoothing: λ > 0.1  
- Less smoothing: λ < 0.1
- Range: 10⁻⁶ to 10³
```

---

## Usage Examples

### Basic Syntax

```bash
# Polynomial fitting (smoothed values only)
./smooth -m 0 -n 7 -p 2 data.txt

# Polynomial fitting with derivatives
./smooth -m 0 -n 7 -p 2 -d data.txt

# Savitzky-Golay (smoothed values only)
# NOTE: Will be rejected if grid is non-uniform!
./smooth -m 1 -n 9 -p 3 data.txt

# Savitzky-Golay with derivatives
./smooth -m 1 -n 9 -p 3 -d data.txt

# Tikhonov with automatic λ (without derivatives)
./smooth -m 2 -l auto data.txt

# Tikhonov with automatic λ and derivatives
./smooth -m 2 -l auto -d data.txt

# Tikhonov with manual λ
./smooth -m 2 -l 0.01 data.txt

# Tikhonov with manual λ and derivatives
./smooth -m 2 -l 0.01 -d data.txt
```

### Output Format

**Without `-d` flag:**
```
# Data smooth - aprox. pol. 2dg from 7 points of moving window (least square)
#    x          y
  0.00000E+00  1.00000E+00
  1.00000E+00  2.71828E+00
  ...
```

**With `-d` flag:**
```
# Data smooth - aprox. pol. 2dg from 7 points of moving window (least square)
#    x          y          y'
  0.00000E+00  1.00000E+00  1.00000E+00
  1.00000E+00  2.71828E+00  2.71828E+00
  ...
```

### Working with Non-uniform Grids (v5.3)

```bash
# Try Savitzky-Golay on non-uniform grid
./smooth -m 1 -n 7 -p 2 nonuniform_data.txt

# If grid is non-uniform (CV > 0.05), you'll see:
# ========================================
# ERROR: Savitzky-Golay method not suitable for non-uniform grid!
# ========================================
# RECOMMENDED ALTERNATIVES:
#   1. Use Tikhonov method: -m 2 -l auto

# Use recommended Tikhonov method instead
./smooth -m 2 -l auto nonuniform_data.txt

# Works perfectly! No additional flags needed in v5.3
```

### Typical Workflow

1. **Quick data exploration:**
   ```bash
   # Try your favorite method
   ./smooth -m 1 data.txt > smooth_data.txt
   
   # If rejected due to non-uniform grid:
   ./smooth -m 2 -l auto data.txt > smooth_data.txt
   ```

2. **Analysis with derivatives:**
   ```bash
   # Add -d flag
   ./smooth -m 2 -l auto -d data.txt > smooth_with_deriv.txt
   ```

3. **Optimal smoothing for very noisy data:**
   ```bash
   # Tikhonov with automatic parameter
   ./smooth -m 2 -l auto data.txt > tikhonov_smooth.txt
   ```

4. **For publication graphics:**
   ```bash
   # Savitzky-Golay if grid is uniform
   ./smooth -m 1 -n 9 -p 3 -d data.txt > publication_data.txt
   
   # Or Tikhonov for any grid type
   ./smooth -m 2 -l auto -d data.txt > publication_data.txt
   ```

---

## Grid Analysis Module

The `grid_analysis` module provides comprehensive analysis of input data and helps optimize smoothing parameters.

### Main Functions

```c
// Complete grid analysis
GridAnalysis* analyze_grid(double *x, int n, int store_spacings);

// Quick uniformity check
int is_uniform_grid(double *x, int n, double *h_avg, double tolerance);

// Method recommendation
const char* get_grid_recommendation(GridAnalysis *analysis);

// Optimal window size
int optimal_window_size(GridAnalysis *analysis, int min_window, int max_window);
```

### GridAnalysis Structure

```c
typedef struct {
    double h_min;           // Minimum spacing
    double h_max;           // Maximum spacing
    double h_avg;           // Average spacing
    double h_std;           // Standard deviation
    double ratio_max_min;   // h_max/h_min ratio
    double cv;              // Coefficient of variation
    double uniformity_score;// Uniformity score (0-1)
    int is_uniform;         // 1 = uniform, 0 = non-uniform
    int n_clusters;         // Number of detected clusters
    int reliability_warning;// Reliability warning
    char warning_msg[512];  // Warning text
} GridAnalysis;
```

### Example Analysis Output

```
# Grid uniformity analysis:
#   n = 1000 points
#   h_min = 1.000000e-02, h_max = 1.000000e-01, h_avg = 5.500000e-02
#   h_max/h_min = 10.00, CV = 0.450
#   Grid type: NON-UNIFORM
#   Uniformity score: 0.35
#   Recommendation: High non-uniformity - use Tikhonov method
```

### Grid Uniformity Thresholds

```
CV < 0.01  : Perfectly uniform - all methods work optimally
CV < 0.05  : Nearly uniform - SAVGOL works with warning, others work fine
CV < 0.20  : Moderately non-uniform - SAVGOL rejected, use POLYFIT or TIKHONOV
CV ≥ 0.20  : Highly non-uniform - TIKHONOV recommended
```

---

## Compilation and Installation

### Requirements

- C compiler (gcc, clang)
- LAPACK and BLAS libraries
- Make (optional)

### Compilation using Make

```bash
# Standard compilation
make

# Debug build
make debug

# Clean
make clean

# Install to user's home directory
make install-user

# Install to system (requires root)
make install
```

### Manual Compilation

```bash
# Standard compilation
gcc -o smooth smooth.c polyfit.c savgol.c tikhonov.c \
    grid_analysis.c decomment.c -llapack -lblas -lm -O2
```

### File Structure

```
smooth/
├── smooth.c           # Main program
├── polyfit.c/h        # Polynomial fitting module
├── savgol.c/h         # Savitzky-Golay module (v5.3: with uniformity check)
├── tikhonov.c/h       # Tikhonov module (v5.3: simplified, correct D²)
├── grid_analysis.c/h  # Grid analysis
├── decomment.c/h      # Comment removal
├── revision.h         # Program version
├── Makefile           # Build system
└── README.md          # This documentation
```

---

## Conclusion

The `smooth` program v5.3 provides three complementary smoothing methods in a modular architecture with advanced input data analysis:

- **POLYFIT** - local polynomial approximation using least squares method
- **SAVGOL** - optimal linear filter with pre-computed coefficients (uniform grids only)
- **TIKHONOV** - global variational method with correct second derivative regularization
- **GRID_ANALYSIS** - automatic analysis and method recommendation

### Version 5.3 Highlights

**Major Improvements:**
1. **Tikhonov method simplified** - removed `adaptive_weights` option, now automatically correct for all grid types
2. **Savitzky-Golay protected** - automatic rejection of non-uniform grids with helpful recommendations
3. **Better mathematical correctness** - proper D² discretization for non-uniform grids
4. **Cleaner API** - removed `tikhonov_smooth_adaptive()`, simplified parameter set
5. **Improved user experience** - automatic method selection guidance

Each method has a strong mathematical foundation and is optimized for specific data types and applications. The program now provides better guidance on which method to use and automatically prevents mathematically incorrect usage patterns.

---

**Document revision:** 2025-10-04  
**Program version:** smooth v5.3  
**Dependencies:** LAPACK, BLAS  
**License:** See source files
