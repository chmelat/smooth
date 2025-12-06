# Methods for Smoothing Experimental Data in the smooth Program

**Technical Documentation**
Version 5.9.2 | December 3, 2025

---

## Contents

1. [Introduction](#introduction)
2. [Polynomial Fitting (POLYFIT)](#polynomial-fitting-polyfit)
3. [Savitzky-Golay Filter (SAVGOL)](#savitzky-golay-filter-savgol)
4. [Tikhonov Regularization (TIKHONOV)](#tikhonov-regularization-tikhonov)
5. [Butterworth Filter (BUTTERWORTH)](#butterworth-filter-butterworth)
6. [Method Comparison](#method-comparison)
7. [Practical Recommendations](#practical-recommendations)
8. [Usage Examples](#usage-examples)
9. [Grid Analysis Module](#grid-analysis-module)
10. [Compilation and Installation](#compilation-and-installation)

---

## Introduction

The `smooth` program implements four sophisticated methods for smoothing experimental data with the capability of simultaneous derivative computation. Each method has specific properties, advantages, and areas of application.

### General Smoothing Problem

Experimental data often contains random noise:
```
y_obs(x_i) = y_true(x_i) + ε_i
```

The goal of smoothing is to estimate `y_true` while suppressing `ε_i` and preserving physically relevant signal properties.

### Program Structure

```bash
./smooth [options] [data_file|-]
```

**Input Options:**
- `data_file` - read data from file
- `-` - read data from stdin (standard input)
- omit argument - read data from stdin (default when no file specified)

**Unix Filter Usage:**
The program can be used as a standard Unix filter in pipe chains:
```bash
cat data.txt | smooth -m 1 -n 5 -p 2       # Pipe input
smooth -m 2 -l 0.01 < input.txt > out.txt  # Redirection
command | smooth -m 3 -f 0.15 | gnuplot    # Pipeline
```

**Basic Parameters:**
- `-m {0|1|2|3}` - method selection (polyfit|savgol|tikhonov|butterworth)
- `-n N` - smoothing window size (polyfit, savgol)
- `-p P` - polynomial degree (polyfit, savgol, max 12)
- `-l λ` - regularization parameter (tikhonov)
- `-l auto` - automatic λ selection using GCV (tikhonov)
- `-f fc` - normalized cutoff frequency (butterworth, 0 < fc < 1.0)
- `-f auto` - automatic cutoff selection (butterworth, currently returns 0.1)
- `-T` - timestamp mode: first column is RFC3339-style timestamp, second is y-value
- `-d` - display first derivative in output (optional, not available for butterworth)
- `-g` - show detailed grid uniformity analysis (optional)

**Note on polynomial degree:** Degrees > 6 may generate numerical stability warnings.

**Note on derivatives:** First derivative output is optional. Without the `-d` switch, the program outputs only smoothed values. With the `-d` switch, it outputs both smoothed values and first derivatives.

**Note on grid analysis:** The `-g` flag provides detailed grid uniformity statistics helpful for understanding your data and choosing appropriate smoothing parameters.

**Note on timestamp mode:** The `-T` flag enables smoothing of time-series data with RFC3339-style timestamps. Timestamps are converted to relative time in seconds for smoothing computations, but the original timestamp format is preserved in output. When combined with `-d`, derivatives are output as dy/dt where t is in seconds.

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

### Least Squares Solution via SVD

The polynomial coefficients are found by solving an **overdetermined linear system** using **Singular Value Decomposition (SVD)**:

```
V · a = y_window

where V is the Vandermonde matrix:

    V[j,k] = (x_j - x_i)^k    for j = 0, …, window_size-1
                                   k = 0, …, poly_degree

         [ 1    (x_0-x_i)    (x_0-x_i)²  ...  (x_0-x_i)^p ]
         [ 1    (x_1-x_i)    (x_1-x_i)²  ...  (x_1-x_i)^p ]
    V =  [ 1    (x_2-x_i)    (x_2-x_i)²  ...  (x_2-x_i)^p ]
         [ ...    ...          ...        ...      ...     ]
         [ 1    (x_n-x_i)    (x_n-x_i)²  ...  (x_n-x_i)^p ]

  a = [a_0, a_1, …, a_p]^T           (polynomial coefficients)
  y_window = [y_{i-n/2}, …, y_{i+n/2}]^T  (data in window)
```

**Why SVD instead of Normal Equations?**

The implementation uses LAPACK's `dgelss` (SVD decomposition) rather than forming normal equations (V^T V)a = V^T y:

1. **Numerical stability:** SVD avoids squaring the condition number (κ(V^T V) = κ(V)²)
2. **Automatic regularization:** Singular values below rcond·σ_max are truncated
3. **Rank detection:** Provides effective rank for diagnosing ill-conditioning
4. **Robustness:** Handles high polynomial degrees (p > 6) more reliably

**SVD truncation parameter:** rcond = 10^{-10}
- Singular values σ_i < 10^{-10}·σ_max are treated as zero
- Provides implicit Tikhonov-style regularization
- Conservative threshold ensures stability without over-regularization

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

The program uses LAPACK routine `dgelss` for solving least squares via SVD at each point:

```c
// Build Vandermonde matrix
build_vandermonde(x, i - offset, i + offset, x[i], poly_degree, V, window_size);

// Solve least squares using SVD decomposition
dgelss_(&window_size, &matrix_cols, &nrhs, V, &window_size,
        rhs, &rhs_size, sing_vals, &rcond, &effective_rank,
        work, &lwork, &info);

// Extract solution: rhs[0] = a_0 (value), rhs[1] = a_1 (derivative)
result->y_smooth[i] = rhs[0];
result->y_deriv[i] = (poly_degree > 0) ? rhs[1] : 0.0;
```

**Numerical diagnostics:**
- On first window, reports condition number: κ = σ_max / σ_min
- If κ > 10^8, issues warning about potential numerical issues
- Reports effective rank if matrix is rank-deficient
- Fallback to original value if SVD fails

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
- **Numerically stable:** SVD decomposition handles ill-conditioned systems
- **Automatic regularization:** Implicit truncation of small singular values
- **Diagnostic feedback:** Reports condition number and effective rank

**Disadvantages:**
- Sensitive to outliers
- Boundary effects at edges
- Possible Runge oscillations for high polynomial degrees (p > 6)
- Numerical instability warnings for degrees > 6 (but handled gracefully by SVD)
- Computationally expensive: O(n·p³) due to per-point SVD

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
Σ_{j=-n_L}^{n_R} c_j · j^m = δ_{m,d} · d!    for m = 0, 1, …, p
```

where:
- `δ_{m,d}` is the Kronecker delta
- `d` is the derivative order
- `d!` is factorial

This leads to a system of linear equations where the unknowns are the filter coefficients `c_j`.

### Matrix Formulation

The coefficients are found by solving a **normal equations system** (not a Vandermonde system):

```
A · β = b

where A is a symmetric (p+1)×(p+1) moment matrix:

    A[i,j] = Σ_{k=-n_L}^{n_R} k^{i+j}    for i,j = 0, 1, …, p

         [ Σk⁰  Σk¹  Σk²  ...  Σk^p   ]
         [ Σk¹  Σk²  Σk³  ...  Σk^(p+1) ]
    A =  [ Σk²  Σk³  Σk⁴  ...  Σk^(p+2) ]
         [ ...  ...  ...  ...  ...     ]
         [ Σk^p Σk^(p+1) ... Σk^(2p) ]

and the right-hand side vector:

    b[j] = δ_{j,d}·d!    (Kronecker delta: 1 if j=d, else 0)
```

This results in a symmetric positive definite (p+1)×(p+1) matrix. The filter coefficients are then:

```
c_k = Σ_{j=0}^p β_j · k^j    for k = -n_L, …, n_R
```

**Note:** This formulation through normal equations is mathematically equivalent to least-squares polynomial fitting but more efficient computationally.

### Computational Efficiency

The brilliance of the Savitzky-Golay approach becomes apparent when processing large datasets:

**Example for 10,000 data points, window size 21, polynomial degree 4:**
- **POLYFIT:** Must solve 10,000 separate 5×5 linear systems
- **SAVGOL:**
  - Solves ONE 5×5 system for central points (pre-computed coefficients)
  - Performs 9,980 simple weighted sums (fast convolution)
  - Solves 20 boundary systems (asymmetric windows at edges)
  - **Net result:** ~500× faster for large datasets

This difference explains why SAVGOL is preferred for real-time signal processing and large datasets, while maintaining the same mathematical accuracy as POLYFIT **for uniform grids**.

**Implementation optimization:** The code pre-computes coefficients for the symmetric window once, then applies them via fast convolution to all central points. Only boundary points require per-point coefficient computation.

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

### Derivative Scaling

**IMPORTANT:** The derivative coefficients computed by `savgol_coefficients()` assume **unit spacing** (normalized integer coordinates). For physical derivatives on real grids, the results must be scaled:

```
dy/dx_physical = (dy/dx_normalized) / h_avg

where h_avg = average grid spacing
```

The implementation automatically performs this scaling:
```c
result->y_deriv[i] = deriv / h_avg;  // Scale to physical units
```

This is **essential** for correct derivative values on uniform grids with spacing ≠ 1.

### Boundary Handling

At data boundaries where a full symmetric window cannot be used, the method employs **asymmetric windows**:

**Central points (i = offset to n-offset):**
- Use symmetric window: nl = nr = offset
- Pre-computed coefficients applied via fast convolution

**Boundary points (left and right edges):**
- Use asymmetric windows: nl ≠ nr
- Coefficients computed per-point for each boundary configuration
- Ensures enough points available for polynomial degree
- Example: leftmost point uses nl=0, nr=window_size-1

**Edge cases:**
- If insufficient points for polynomial degree, falls back to original value
- Maintains polynomial exactness property at boundaries
- More computationally expensive than central points (acceptable for small boundary regions)

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
       _==================     _=============
       Data fidelity term      Smoothness penalty
```

**Discrete formulation:**
```
J[u] = ||y - u||² + λ||D²u||²
```

where:
- `||y - u||²` = Σ(y_i - u_i)² is the **data fidelity term**
- `||D²u||²` = Σ(D²u_i)² is the **regularization term** (smoothness penalty)
- `λ` is the **regularization parameter** controlling the balance
- `D²` is the discrete second derivative operator

### The Regularization Parameter λ

The parameter λ is the **heart of Tikhonov regularization** - it controls the balance between fitting the data and smoothing the result.

#### Physical Interpretation

```
λ = 0:     No smoothing, u = y (exact data fit)
           J[u] = ||y - u||² only

λ → ∞:     Maximum smoothing, u → straight line
           J[u] ≈ λ||D²u||² dominates

λ optimal: Balanced between data fit and smoothness
           Both terms contribute meaningfully
```

#### Mathematical Role

The minimization of J[u] leads to:
```
(I + λD^TD)u = y
```

**Effect of λ on the solution:**
- **Small λ (< 0.01):** Matrix ≈ I → solution u ≈ y (minimal smoothing)
- **Large λ (> 1.0):** Matrix ≈ λD^TD → strong curvature penalty (heavy smoothing)
- **Optimal λ:** Matrix components balanced → noise removed, signal preserved

#### Frequency Domain Interpretation

In Fourier space, Tikhonov acts as a low-pass filter:
```
Ĥ(ω) = 1 / (1 + λω⁴)
```

where ω is spatial frequency. 

**Effect:**
- **Low frequencies (slow variations):** Ĥ ≈ 1 → preserved
- **High frequencies (noise, rapid variations):** Ĥ ≈ 1/(λω⁴) → attenuated
- **Cutoff frequency:** ω_c ~ λ^(-1/4)

**This means:**
```
Larger λ  → Lower cutoff → More aggressive low-pass filtering → Smoother result
Smaller λ → Higher cutoff → Less filtering → Result closer to data
```

#### Practical Guidelines for λ Selection

**1. Automatic Selection (RECOMMENDED):**
```bash
./smooth -m 2 -l auto data.txt
```
Uses Generalized Cross Validation (GCV) to find optimal λ.

**2. Manual Selection:**

| Data Characteristics | Recommended λ | Reasoning |
|---------------------|---------------|-----------|
| Low noise, important details | 0.001 - 0.01 | Preserve features |
| Moderate noise | 0.01 - 0.1 | Balanced (default: 0.1) |
| High noise | 0.1 - 1.0 | Strong smoothing |
| Very noisy, global trends | 1.0 - 10.0 | Maximum smoothing |

**3. Iterative Refinement:**
```bash
# Start with automatic
./smooth -m 2 -l auto data.txt

# If result is over-smoothed (details lost):
./smooth -m 2 -l 0.01 data.txt

# If result is under-smoothed (still noisy):
./smooth -m 2 -l 1.0 data.txt
```

**4. Diagnostic Criteria:**

The program outputs functional components:
```
# Functional J = 1.234e+02 (Data: 5.67e+01 + Regularization: 6.67e+01)
# Data/Total ratio = 0.460, Regularization/Total ratio = 0.540
```

**Good balance indicators:**
- Data term: 30-70% of total functional
- Regularization term: 30-70% of total functional

**Warning signs:**
- Data term > 95%: Under-smoothed (λ too small)
- Regularization term > 95%: Over-smoothed (λ too large)

**5. Grid-Dependent Considerations:**

For highly non-uniform grids (CV > 0.2):
```bash
# Start with more conservative (larger) λ
./smooth -m 2 -l 0.5 nonuniform_data.txt

# GCV may be less accurate - check results visually
./smooth -m 2 -l auto nonuniform_data.txt
```

#### λ and Grid Spacing

The effective regularization strength depends on grid spacing:
```
Effective strength ~ λ / h²_avg

Same λ on finer grid   → weaker smoothing
Same λ on coarser grid → stronger smoothing
```

For dimensional consistency, λ has units [Length²].

### Second Derivative Discretization (Hybrid Method)

Automatic selection between two discretization schemes based on grid uniformity.

#### Grid Uniformity Detection
```
CV = coefficient of variation (h_std / h_avg)

CV < 0.15:  Nearly uniform    → Average Coefficient Method
CV ≥ 0.15:  Highly non-uniform → Local Spacing Method
```

#### Method 1: Average Coefficient (for CV < 0.15)

Used for uniform and mildly non-uniform grids. More robust numerically.

**Discretization:**
For interior point `i` with neighbors at spacing h_left and h_right, use **harmonic mean**:
```
h_harm = 2·h_left·h_right / (h_left + h_right)

D²u_i ≈ (u_{i-1} - 2u_i + u_{i+1}) / h_harm²
```

**Why harmonic mean?** 
- More accurate than arithmetic mean for averaging intervals
- Gives greater weight to smaller spacing (physically correct)
- For h_left = h_right, reduces to standard formula

**Matrix construction:**
```
c = λ · Σ(1/h_i²) / (n-1)    (average coefficient)

Tridiagonal matrix A = I + λD^TD:

    A[i,i]   = 1 + 2c     (diagonal, interior points)
    A[i,i±1] = -c         (off-diagonals)

Example for n=5 points:

         [ 1+2c   -c     0     0     0  ]
         [  -c   1+2c   -c     0     0  ]
    A =  [  0     -c   1+2c   -c     0  ]
         [  0     0     -c   1+2c   -c  ]
         [  0     0     0     -c   1+2c ]
```

#### Method 2: Local Spacing (for CV ≥ 0.15)

Used for highly non-uniform grids. More accurate for variable spacing.

**Discretization:**
For point `i` with left spacing `h₁ = x[i] - x[i-1]` and right spacing `h₂ = x[i+1] - x[i]`:

```
D²u_i ≈ (2/(h₁+h₂)) · [u_{i-1}/h₁ - u_i·(1/h₁+1/h₂) + u_{i+1}/h₂]
```

This is the **correct second derivative formula** for non-uniform grids derived from Taylor expansion.

**Matrix construction:**
```
For each interior point i, compute weight:
    w_i = 2λ / (h_{i-1} + h_i)

where h_{i-1} = x[i] - x[i-1], h_i = x[i+1] - x[i]

Tridiagonal matrix A = I + λD^TD:

    A[i,i]   = 1 + w_i·(1/h_{i-1} + 1/h_i)    (diagonal)
    A[i,i-1] = -w_i/h_{i-1}                   (lower diagonal)
    A[i,i+1] = -w_i/h_i                       (upper diagonal)

Example structure for n=5 points (non-uniform spacing):

         [  d0    u0     0      0      0   ]
         [  l1    d1    u1      0      0   ]
    A =  [  0     l2    d2     u2      0   ]
         [  0     0     l3     d3     u3   ]
         [  0     0      0     l4     d4   ]

where d_i, u_i, l_i vary with local spacing
```

**The resulting matrix is:**
- Symmetric (u_i = l_{i+1})
- Positive definite
- Tridiagonal (bandwidth = 1)

#### Boundary Conditions

Natural boundary conditions (second derivative = 0 at ends):

**Matrix construction uses simplified two-point formula:**

**Left boundary (i=0):**
```
Penalty: λ·[(u_1 - u_0) / h_0²]²

A[0,0] += λ/h_0²
A[0,1] += -λ/h_0²
```

**Right boundary (i=n-1):**
```
Penalty: λ·[(u_{n-1} - u_{n-2}) / h_{n-1}²]²

A[n-1,n-1] += λ/h_{n-1}²
```

**Implementation note:** The boundary superdiagonal element A[0,1] is properly included to prevent isolation of the first point.

**Functional computation uses more accurate three-point formula** (as shown in the Functional Computation section above) for better accuracy when evaluating the objective function value.

#### Boundary Effects and Edge Artifacts

**Important:** Natural boundary conditions can cause **oscillations near data endpoints**, especially on non-uniform grids with small λ values.

**Observed behavior:**
- **Uniform grids (CV < 0.15):** Minimal boundary effects, typically < 5% deviation
- **Non-uniform grids (CV > 0.15):** Significant boundary artifacts possible:
  - Last 2-3 points may show deviations up to 30-60% from expected values
  - Effect increases with grid non-uniformity and decreases with larger λ
  - More pronounced with Local Spacing Method discretization

**Example (documented in unit tests):**
```
Grid: CV = 0.176, N = 100, λ = 0.01
Expected y_max ≈ 32 (parabolic function)
Observed y_max ≈ 52 at last point (+61% overshoot)
```

**Practical recommendations:**

1. **Discard edge points in analysis:**
   ```bash
   # Process data and remove first/last 3 points
   ./smooth -m 2 -l 0.01 data.txt | tail -n +4 | head -n -3
   ```

2. **Use larger λ for stability:**
   ```bash
   # Increase λ to reduce boundary oscillations
   ./smooth -m 2 -l 0.1 data.txt    # Instead of 0.01
   ```

3. **Add padding data:**
   - Extend dataset by extrapolating 5-10 points at each end
   - Apply smoothing to extended data
   - Use only interior region of smoothed result

4. **Use GCV with caution on non-uniform grids:**
   ```bash
   # GCV may select λ too small for stable boundaries
   ./smooth -m 2 -l auto data.txt

   # Verify boundary behavior visually
   # Consider manually increasing λ if needed
   ```

5. **Grid-specific guidelines:**
   - **CV < 0.05:** Edge artifacts negligible (< 2%)
   - **CV 0.05-0.15:** Monitor last 1-2 points (< 10% typical)
   - **CV > 0.15:** **Discard last 2-3 points** (artifacts up to 60%)

**Why this happens:**
- Natural boundary conditions impose `u'' = 0` at endpoints
- On non-uniform grids, this constraint conflicts with data fidelity
- Small λ amplifies this conflict (weak regularization)
- Result: The solver "overshoots" to satisfy both constraints

**This is not a bug** - it is inherent to the variational formulation with natural boundary conditions on non-uniform domains.

### Functional Computation

The actual value of the minimized functional is computed for diagnostic purposes:

**Data term:**
```
||y - u||² = Σ(y_i - u_i)²
```

**Regularization term:**

For **average coefficient method:**
```
||D²u||² = Σ_{interior} [(u_{i-1} - 2u_i + u_{i+1})/h_harm²]²
         + 0.5·[D²u_left]²
         + 0.5·[D²u_right]²
```

where boundary second derivatives use forward/backward three-point formulas:
```
D²u_left  = 2·(h₁·u₀ - (h₀+h₁)·u₁ + h₀·u₂) / (h₀·h₁·(h₀+h₁))
D²u_right = 2·(h_{n-2}·u_{n-1} - (h_{n-3}+h_{n-2})·u_{n-2} + h_{n-3}·u_{n-3}) / (h_{n-3}·h_{n-2}·(h_{n-3}+h_{n-2}))
```

For **local spacing method:**
```
||D²u||² = Σ_{interior} [D²u_i]² · (h₁+h₂)/2
         + 0.5·[D²u_left]² · h₀
         + 0.5·[D²u_right]² · h_{n-1}
```

where the boundary terms use the same three-point formulas as above, and weighting factors ensure proper integration over non-uniform grid.

**Total functional:**
```
J[u] = ||y - u||² + λ||D²u||²
```

### Variational Approach

The minimum of functional J[u] satisfies the Euler-Lagrange equation:

```
∂J/∂u_i = 0  ==>  -2(y_i - u_i) + 2λ(D^T D u)_i = 0
```

which leads to the linear system:
```
(I + λD^T D)u = y
```

### Matrix Representation

The linear system `(I + λD^T D)u = y` has matrix `A = I + λD^T D` with structure:

**Properties:**
- Symmetric
- Positive definite
- Tridiagonal (banded with bandwidth 1)

**General tridiagonal form:**
```
         [  d0    c0     0      0     ...   0   ]
         [  c0    d1    c1      0     ...   0   ]
         [   0    c1    d2     c2     ...   0   ]
    A =  [   0     0    c2     d3     ...   0   ]
         [ ...   ...   ...    ...     ...  ...  ]
         [   0     0     0      0     c_{n-2} d_{n-1} ]

where d_i = diagonal, c_i = off-diagonal elements
```

This structure allows efficient solution using LAPACK's banded solver `dpbsv`.

### Generalized Cross Validation (GCV)

For automatic λ selection (`-l auto`), we minimize the GCV criterion:

```
GCV(λ) = n·RSS(λ) / (n - tr(H_λ))²
```

where:
- `RSS(λ) = ||y - u_λ||²` is the residual sum of squares
- `H_λ = (I + λD^T D)^{-1}` is the influence matrix (smoother matrix)
- `tr(H_λ)` is the trace (effective number of parameters)

**Interpretation:**
- `tr(H_λ)` measures model complexity (degrees of freedom)
- Small λ: tr(H) ≈ n (interpolation, overfitting)
- Large λ: tr(H) ≈ 2 (straight line, underfitting)
- Optimal λ: minimizes prediction error

**Trace estimation using eigenvalues:**

For uniform grids with natural boundary conditions:
```
tr(H_λ) ≈ Σ_{k=1}^n 1/(1 + λμ_k)

where eigenvalues:
θ_k = πk/n
μ_k = 4·sin²(θ_k/2) / h²
```

**Note:** This approximation is exact for uniform grids but approximate for non-uniform grids. For highly non-uniform grids (CV > 0.2), the program issues a warning.

#### Enhanced GCV

**Over-fitting penalty:**
```
If tr(H)/n > 0.7:
    GCV_modified = GCV · exp(10·(tr(H)/n - 0.7))
```

This exponential penalty prevents selection of too-small λ that would lead to overfitting.

**L-curve backup (for n > 20000):**

For very large datasets, GCV trace approximation may be inaccurate. The program also computes the L-curve (plot of ||D²u||² vs ||y-u||²) and finds the corner point with maximum curvature:

```
κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)

where:
x = log(||y - u||²)
y = log(||D²u||²)
```

If GCV and L-curve disagree significantly, the program uses the more conservative (larger) λ.

### Efficient Implementation

The program uses LAPACK routine `dpbsv` for solving symmetric positive definite banded systems:

```c
// Banded matrix storage (LAPACK column-major format)
// For tridiagonal symmetric matrix A with bandwidth kd=1:
//
//     [ AB[0,j] ]  = superdiagonal elements  (a[i,i+1])
//     [ AB[1,j] ]  = diagonal elements       (a[i,i])
//
// Storage layout for tridiagonal matrix:
//
//         [ *    a01  a12  a23  a34  ... ]   <- row 0 (superdiagonal)
//    AB = [ a00  a11  a22  a33  a44  ... ]   <- row 1 (diagonal)
//
// System solution
dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
```

**Complexity:**
- Memory: O(n) for banded storage
- Time: O(n) for factorization and back-substitution

This is **optimal** for tridiagonal systems.

### Implementation Details

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

// Main function
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda);

// Automatic λ selection
double find_optimal_lambda_gcv(double *x, double *y, int n);

// Memory cleanup
void free_tikhonov_result(TikhonovResult *result);
```

### Characteristics

**Advantages:**
- Global optimization with theoretical foundation
- Flexible balance between data fidelity and smoothness (controlled by λ)
- Robust to outliers (quadratic penalty less sensitive than least squares)
- Efficient for large datasets (O(n) memory and time)
- **Automatic λ selection via GCV** - no guessing needed
- **Excellent for non-uniform grids** - correct discretization automatic
- **Unified approach** - same algorithm for uniform and non-uniform grids
- Works well for noisy data with global trends

**Disadvantages:**
- Single global parameter λ (cannot vary locally)
- May suppress local details if λ too large
- GCV may fail for some data types (especially highly non-uniform grids)
- Requires LAPACK library
- Boundary effects if data has discontinuities at edges

---

## Butterworth Filter (BUTTERWORTH)

### Theoretical Foundation

The Butterworth filter is a classical **low-pass frequency filter** in digital signal processing (DSP). It removes high-frequency noise while preserving low-frequency signal trends.

**What does "low-pass" mean?**
- **Passes low frequencies:** Slow variations in your data pass through unchanged
- **Blocks high frequencies:** Rapid fluctuations (noise) are removed
- **The cutoff frequency (fc)** determines the boundary between "low" and "high"
  - Lower fc → more aggressive smoothing (removes more detail)
  - Higher fc → gentler smoothing (preserves more detail)

The filter is characterized by a **maximally flat magnitude response** in the passband and provides zero phase distortion when implemented as filtfilt.

**Filter Transfer Function:**

In the analog domain (s-domain), the Butterworth filter has magnitude response:
```
|H(jω)|² = 1 / (1 + (ω/ωc)^(2N))
```

where:
- `N` = filter order (4 in our implementation)
- `ωc` = cutoff frequency (3dB point)
- `ω` = frequency

**Key Properties:**
- **Maximally flat passband:** No ripples for ω < ωc
- **Monotonic rolloff:** Smooth transition from passband to stopband
- **-3dB at cutoff:** |H(jωc)| = 1/√2 ≈ 0.707
- **Rolloff rate:** -20N dB/decade (for N=4: -80 dB/decade)

### Digital Implementation

The smooth program implements a **4th-order digital Butterworth low-pass filter** using the following algorithm:

**Step 1: Pole Calculation**

Butterworth poles lie on unit circle in s-domain at angles:
```
θk = π/2 + π(2k+1)/(2N),  k = 0, 1, …, N−1
```

For N=4:
```
s_poles[k] = exp(j·θk)  where θ = {5π/8, 7π/8, 9π/8, 11π/8}
```

**Step 2: Frequency Scaling**

Scale poles by prewarped cutoff frequency:
```
wc = tan(π·fc/2)    (prewarp for bilinear transform)
s_poles_scaled = wc · s_poles
```

**Prewarping correction:** The bilinear transform introduces frequency warping. The factor `tan(π·fc/2)` compensates for this, ensuring the digital filter's cutoff matches the desired normalized frequency `fc`.

**Step 3: Bilinear Transform**

Convert analog poles to digital domain:
```
z_poles = (2 + s_poles_scaled) / (2 - s_poles_scaled)
```

The bilinear transformation maps:
- Left half of s-plane → inside unit circle in z-plane
- jω axis → unit circle in z-plane
- Preserves stability

**Step 4: Biquad Cascade**

Form two 2nd-order sections (biquads) from conjugate pole pairs:
```
H(z) = H1(z) · H2(z)

Each biquad: H_i(z) = (b0 + b1·z⁻¹ + b2·z⁻²) / (1 + a1·z⁻¹ + a2·z⁻²)
```

This approach provides better numerical stability than direct 4th-order implementation.

### Filtfilt Algorithm

The **filtfilt** (forward-backward filtering) eliminates phase distortion:

**Algorithm:**
1. **Pad signal:** Reflect signal at boundaries (3×order length)
2. **Forward filter:** Apply H(z) from left to right → y_fwd
3. **Reverse:** y_rev = reverse(y_fwd)
4. **Backward filter:** Apply H(z) to y_rev → y_bwd
5. **Reverse back:** y_final = reverse(y_bwd)
6. **Extract:** Remove padding to get final result

**Effect:**
- **Zero phase lag:** No signal delay
- **Effective order:** 2N = 8 (squared magnitude response)
- **Steeper rolloff:** |H_eff(jω)|² = |H(jω)|⁴

### Initial Conditions (lfilter_zi)

To minimize edge transients, we compute initial filter state using **scipy's lfilter_zi algorithm**:

**Problem:** Find initial state `zi` such that for constant input `x = c`:
```
zi = A·zi + B·c
```

This ensures the filter starts in steady-state, eliminating startup transients.

**Solution:** Solve linear system using companion matrix:
```
(I - A^T)·zi = B

where:
  A = companion(a).T    (companion matrix transposed)
  B = b[1:] - a[1:]·b[0]
```

The companion matrix for filter coefficients `[1, a1, a2, a3, a4]` is:

```
         [ -a1  -a2  -a3  -a4 ]
         [  1    0    0    0  ]
    C =  [  0    1    0    0  ]
         [  0    0    1    0  ]
```

And its transpose `A^T = C^T`:

```
          [ -a1   1    0    0 ]
          [ -a2   0    1    0 ]
    A^T = [ -a3   0    0    1 ]
          [ -a4   0    0    0 ]
```

Therefore the system matrix `I - A^T`:

```
              [ 1+a1  -1    0    0 ]
              [  a2    1   -1    0 ]
    I - A^T = [  a3    0    1   -1 ]
              [  a4    0    0    1 ]
```

**Implementation:** The linear system is solved using **LAPACK's dgesv** routine for robustness:
```c
// Uses LU decomposition with partial pivoting
dgesv_(&n, &nrhs, A_colmajor, &lda, ipiv, zi, &ldb, &info);
```

This approach ensures numerical stability even for challenging filter coefficients (e.g., very low cutoff frequencies where coefficients can be extremely small ~10⁻⁸).

### Normalized Cutoff Frequency

The cutoff frequency `fc` is **normalized** to the sampling rate and is the **most important parameter** for Butterworth filtering.

**Simple explanation:**
- `fc` controls how much smoothing you get
- **Smaller fc (e.g., 0.05)** → heavy smoothing, only very slow trends preserved
- **Larger fc (e.g., 0.30)** → light smoothing, more detail preserved
- Valid range: `0 < fc < 1.0` (Nyquist limit)

**Technical details:**
```
fc = f_cutoff / f_Nyquist = f_cutoff / (f_sample / 2)

where:
  f_cutoff = desired cutoff frequency in physical units
  f_sample = 1 / h_avg  (h_avg = average data spacing)
  f_Nyquist = f_sample / 2  (Nyquist frequency)
```

**Nyquist Constraint:** `0 < fc < 1.0`
- `fc = 1.0` corresponds to Nyquist frequency (f_sample/2) - maximum possible
- Higher fc -> less filtering (more high frequencies pass)
- Lower fc -> more filtering (smoother result)

**Physical Interpretation:**

Example: Data with spacing h_avg = 0.1 seconds
- Sample rate: f_sample = 1/0.1 = 10 Hz
- Nyquist frequency: f_Nyquist = 5 Hz
- If fc = 0.2, then f_cutoff = fc × f_Nyquist = 0.2 × 5 = 1 Hz
- Filter removes frequencies above ~1 Hz

**Practical Guidelines for Choosing fc:**

| fc Value | Smoothing Strength | When to Use |
|----------|-------------------|-------------|
| 0.01 - 0.05 | **Very strong** | Extremely noisy data, only global trends matter |
| 0.05 - 0.15 | **Moderate** | Typical experimental data with noise |
| 0.15 - 0.30 | **Light** | Good quality data, preserve features |
| > 0.30 | **Minimal** | Low noise, want to keep almost everything |

**Quick Start Recommendations:**
- **Not sure? Start with fc = 0.15** - good balance for most data
- **Too noisy after smoothing?** Decrease fc (e.g., try 0.10)
- **Lost important details?** Increase fc (e.g., try 0.25)
- **Extreme noise?** Try fc = 0.05
- **High quality data?** Try fc = 0.25 - 0.30

### IIR Filter Implementation

Uses **Transposed Direct Form II** for numerical stability:

```
For each sample n:
  y[n] = b[0]·x[n] + z[0]
  z[0] = b[1]·x[n] - a[1]·y[n] + z[1]
  z[1] = b[2]·x[n] - a[2]·y[n] + z[2]
  z[2] = b[3]·x[n] - a[3]·y[n] + z[3]
  z[3] = b[4]·x[n] - a[4]·y[n]
```

where `z[]` is the filter state (4 elements for 4th order).

### Modularized Implementation

```c
// butterworth.h
typedef struct {
    double *y_smooth;     // Smoothed values
    int n;                // Number of points
    int order;            // Filter order (4)
    double cutoff_freq;   // Normalized cutoff frequency
    double sample_rate;   // Effective sample rate (1/h_avg)
} ButterworthResult;

// Main function
ButterworthResult* butterworth_filtfilt(double *x, double *y, int n,
                                        double cutoff_freq, int auto_cutoff);

// Automatic cutoff selection (currently returns 0.1)
double estimate_cutoff_frequency(double *x, double *y, int n);

// Memory cleanup
void free_butterworth_result(ButterworthResult *result);
```

### Grid Requirements

**IMPORTANT:** Butterworth filter works best with **uniform or nearly-uniform grids**.

The filter assumes uniform sampling when computing the cutoff frequency. For highly non-uniform grids, the program checks grid uniformity before applying the filter.

**Why uniform grids?**
- Frequency analysis assumes constant sampling rate
- fc is defined relative to sample rate
- Non-uniform sampling distorts frequency response

**For non-uniform grids:** Use Tikhonov method (`-m 2 -l auto`) which handles arbitrary spacing correctly.

### Characteristics

**Advantages:**
- **Zero phase distortion** (filtfilt eliminates all phase lag)
- **Maximally flat frequency response** in passband
- **Classical DSP approach** with extensive literature and understanding
- **Predictable frequency-domain behavior** - easy to interpret cutoff frequency
- **No ringing** (unlike Chebyshev or elliptic filters)
- **Efficient implementation** - O(n) time complexity
- **Smooth monotonic rolloff** - natural attenuation curve

**Disadvantages:**
- **Requires uniform/nearly-uniform grid** (CV < 0.15 recommended)
- **No derivative output** (Butterworth is smoothing-only)
- **Less local adaptability** than polynomial methods
- **Cutoff selection not automatic** (currently manual tuning needed)
- **Edge effects** despite padding
- **Frequency interpretation** may be less intuitive than λ for some users

### Comparison with Other Methods

**BUTTERWORTH vs SAVITZKY-GOLAY:**
- Both assume uniform grids
- **Butterworth:** True frequency-domain filtering, maximally flat passband
- **Savitzky-Golay:** Polynomial approximation in time domain
- **Choose Butterworth for:** Periodic signals, spectral data, frequency-domain interpretation
- **Choose Savitzky-Golay for:** Polynomial trends, peak detection, derivative estimation

**BUTTERWORTH vs TIKHONOV:**
- **Butterworth:** Classical signal processing, frequency-domain control
- **Tikhonov:** Variational optimization, works with non-uniform grids
- **Choose Butterworth for:** Uniform data, need frequency-domain understanding
- **Choose Tikhonov for:** Non-uniform grids, mathematical optimization approach

**BUTTERWORTH vs POLYFIT:**
- **Butterworth:** Global frequency filtering, uniform smoothing
- **Polyfit:** Local polynomial fitting, adapts to curvature changes
- **Choose Butterworth for:** Stationary signals, spectroscopic data
- **Choose Polyfit for:** Variable curvature, local feature preservation

---

## Method Comparison

### Computational Complexity

| Method | Time | Memory | Scalability |
|--------|------|--------|-------------|
| POLYFIT | O(n·p³) | O(p²) | Good for small p |
| SAVGOL | O(p³) + O(n·w) | O(w) | Excellent for large n |
| TIKHONOV | O(n) | O(n) | Excellent |
| BUTTERWORTH | O(n) | O(n) | Excellent |

*Note: w = window size, p = polynomial degree (≤12), n = number of data points.*

### Smoothing Quality

| Property | POLYFIT | SAVGOL | TIKHONOV | BUTTERWORTH |
|----------|---------|--------|----------|-------------|
| Local adaptability | ***** | **** | ** | ** |
| Extreme preservation | **** | ***** | *** | *** |
| Noise robustness | *** | **** | ***** | ***** |
| Derivative quality | ***** | ***** | *** | N/A |
| Boundary behavior | ** | *** | **** | *** |
| Non-uniform grids | *** | [X] | ***** | ** |
| Ease of use | **** | **** | ***** | **** |
| Parameter selection | Manual | Manual | Auto (GCV) | Manual |
| Frequency control | No | No | No | Yes |
| Phase distortion | N/A | N/A | N/A | Zero |

**Key:** [X] = Not suitable (automatically rejected)

### Grid Type Compatibility

| Grid Type | POLYFIT | SAVGOL | TIKHONOV | BUTTERWORTH |
|-----------|---------|--------|----------|-------------|
| Uniform (CV <= 0.01) | [OK] | [OK] | [OK] | [OK] |
| Nearly uniform (CV < 0.05) | [OK] | [WARNING] | [OK] | [OK] |
| Moderately non-uniform (0.05 < CV < 0.2) | [OK] | [X] | [OK] | [WARNING] |
| Highly non-uniform (CV > 0.2) | [WARNING] | [X] | [OK] | [WARNING] |

**Legend:**
- [OK] = Recommended
- [WARNING] = Usable with caution
- [X] = Rejected or not recommended
- \* = Uses local spacing method automatically

---

## Practical Recommendations

### Method Selection by Data Type

#### POLYFIT - when:
- Data has variable curvature
- You need to preserve local details
- You have moderately noisy data
- You want highest quality derivatives
- Grid has moderate spacing variations
- You need local adaptability

#### SAVGOL - when:
- **Grid is uniform (CV < 0.05)** - automatically checked!
- You want mathematically optimal linear smoothing for polynomial signals
- Data contains periodic or oscillatory components that need preservation
- You need excellent peak shape preservation (areas, moments)
- You want minimal phase distortion in the smoothed signal
- You're processing time series or spectroscopic data on uniform grids
- You need simultaneous high-quality function and derivative estimation
- Computational efficiency is critical (large datasets)

#### TIKHONOV - when:
- **Grid is non-uniform** - works perfectly automatically!
- Data is very noisy
- You need global consistency
- **You want automatic parameter selection (λ auto)** - highly recommended!
- You prefer global optimization approaches over local fitting
- You want robust handling of outliers
- You want the simplest workflow (one parameter, automatic selection)
- You need to process very large datasets efficiently
- **You're not sure which method to use** - Tikhonov with `-l auto` is safest!

#### BUTTERWORTH - when:
- **Grid is uniform or nearly-uniform (CV < 0.15)** - essential requirement!
- You want to **remove high-frequency noise** while keeping slow trends
- You need **simple frequency-based smoothing** - just set cutoff frequency fc
- You want **zero phase distortion** (no signal delay)
- Data is periodic, oscillatory, or spectroscopic
- You need **frequency-domain interpretation** of filtering
- Data is from instrumentation with known sampling rate
- You need **predictable frequency response** (maximally flat passband)
- Working with time-series data at constant sampling
- You understand or want to learn about cutoff frequency concept

### Parameter Selection

#### Window size (n) for POLYFIT/SAVGOL:
```
n = 2*k + 1    (odd number)

Recommendations:
- Low noise: n = 5-9
- Medium noise: n = 9-15  
- High noise: n = 15-25

Rule of thumb: n ≈ 2p + 3
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
**RECOMMENDED:**
- Auto selection: -l auto  (uses GCV optimization)

**MANUAL SELECTION:**
Starting points by noise level:
- Low noise:        λ = 0.001 - 0.01
- Medium noise:     λ = 0.01 - 0.1   (default: 0.1)
- High noise:       λ = 0.1 - 1.0
- Very noisy data:  λ = 1.0 - 10.0

Full range: 10⁻⁶ to 10³

**ITERATIVE REFINEMENT:**
1. Start with -l auto
2. Check functional balance (should be 30-70% each)
3. If over-smoothed: decrease λ by factor of 10
4. If under-smoothed: increase λ by factor of 10
5. Repeat until satisfied

**GRID-DEPENDENT:**
For non-uniform grids (ratio > 5):
- Start more conservative (larger λ)
- GCV may be less accurate - check visually
```

#### Cutoff frequency (fc) for BUTTERWORTH:
```
**RECOMMENDED:**
Start with manual selection: fc = 0.15 - 0.20

**MANUAL SELECTION by noise level:**
- Low noise:        fc = 0.20 - 0.30  (preserve details)
- Medium noise:     fc = 0.15 - 0.20  (typical, recommended)
- High noise:       fc = 0.05 - 0.15  (aggressive smoothing)
- Very noisy data:  fc = 0.01 - 0.05  (heavy smoothing)

Full range: 0 < fc < 1.0 (Nyquist limit)

**AUTOMATIC SELECTION:**
- Use -f auto (currently returns default fc = 0.1)
- Note: Automatic selection not yet fully implemented
- Manual tuning recommended for best results

**PHYSICAL INTERPRETATION:**
fc = f_cutoff / f_Nyquist = f_cutoff / (f_sample / 2)
where f_sample = 1 / h_avg

Example: h_avg = 0.1 sec -> f_sample = 10 Hz, f_Nyquist = 5 Hz
         fc = 0.2 -> f_cutoff = fc × f_Nyquist = 1 Hz (removes freq > 1 Hz)

**ITERATIVE REFINEMENT:**
1. Start with fc = 0.15 or fc = 0.20
2. If result too smooth (details lost): increase fc
3. If result too noisy (not smooth enough): decrease fc
4. Typical adjustment: ±0.05
5. Repeat until satisfied

**GRID-DEPENDENT:**
For non-uniform grids (CV > 0.05):
- Results may be suboptimal
- Consider using Tikhonov instead
- If CV > 0.15: use caution, Tikhonov recommended
```

---

## Usage Examples

### Basic Syntax

```bash
# Read from file
./smooth -m 0 -n 7 -p 2 data.txt

# Read from stdin (pipe)
cat data.txt | ./smooth -m 0 -n 7 -p 2

# Read from stdin (explicit)
./smooth -m 0 -n 7 -p 2 -

# Read from stdin (redirection)
./smooth -m 0 -n 7 -p 2 < data.txt
```

### Method Examples

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

# Tikhonov with automatic λ (RECOMMENDED)
./smooth -m 2 -l auto data.txt

# Tikhonov with automatic λ and derivatives
./smooth -m 2 -l auto -d data.txt

# Tikhonov with manual λ
./smooth -m 2 -l 0.01 data.txt

# Tikhonov with manual λ and derivatives
./smooth -m 2 -l 0.01 -d data.txt

# Butterworth with manual cutoff frequency
./smooth -m 3 -f 0.15 data.txt

# Butterworth with automatic cutoff (currently returns 0.1)
./smooth -m 3 -f auto data.txt

# Grid analysis only (exits after analysis)
./smooth -g data.txt
```

### Timestamp Mode Examples

```bash
# Tikhonov smoothing with automatic lambda selection
./smooth -T -m 2 -l auto timeseries.dat

# Tikhonov with derivatives (dy/dt in seconds)
./smooth -T -m 2 -l 0.01 -d timeseries.dat

# Savitzky-Golay filter (requires nearly uniform time spacing)
./smooth -T -m 1 -n 5 -p 2 sensor_data.txt

# Polyfit with derivatives
./smooth -T -m 0 -n 7 -p 3 -d measurements.csv

# Butterworth filter (no derivatives available)
./smooth -T -m 3 -f 0.15 signal.dat
```

**Input format for timestamp mode:**
```
# Space separator format
2025-09-25 14:06:06.390  0.02128
2025-09-25 14:06:06.391  0.02110
2025-09-25 14:06:06.763  0.02230

# T separator format (RFC3339)
2025-09-25T14:06:06.390  0.02128
2025-09-25T14:06:06.391  0.02110
2025-09-25T14:06:06.763  0.02230
```

### Output Format

**Without `-d` flag:**
```
# Data smooth - Tikhonov regularization with lambda = 1e-01
# Functional J = 1.234e+02 (Data: 5.67e+01 + Regularization: 6.67e+01)
# Data/Total ratio = 0.460, Regularization/Total ratio = 0.540
#    x          y
  0.00000E+00  1.00000E+00
  1.00000E+00  2.71828E+00
  ...
```

**With `-d` flag:**
```
# Data smooth - Tikhonov regularization with lambda = 1e-01
# Functional J = 1.234e+02 (Data: 5.67e+01 + Regularization: 6.67e+01)
# Data/Total ratio = 0.460, Regularization/Total ratio = 0.540
#    x          y          y'
  0.00000E+00  1.00000E+00  1.00000E+00
  1.00000E+00  2.71828E+00  2.71828E+00
  ...
```

**Timestamp mode output (with `-T` flag):**
```
# Data smooth - Tikhonov regularization with lambda = 1e-02
# Functional J = 1.07e-03 (Data: 8.33e-04 + Regularization: 2.39e-04)
# Data/Total ratio = 0.777, Regularization/Total ratio = 0.223
# Derivative units: dy/dt (t in seconds)
#    timestamp          y          y'
2025-09-25 14:06:06.390 0.000816394 -0.00204621
2025-09-25 14:06:06.391 0.000814348  0.0542818
2025-09-25 14:06:06.763  0.0210635  0.0284901
  ...
```

**Note:** In timestamp mode, the original timestamp format from input is preserved exactly in output.

**With `-g` flag (grid analysis):**
```
# ========================================
# GRID UNIFORMITY ANALYSIS
# ========================================
# Grid uniformity analysis:
#   n = 1000 points
#   h_min = 9.500000e-03, h_max = 1.200000e-02, h_avg = 1.000000e-02
#   CV = 0.052
#   Grid type: NON-UNIFORM
#   Uniformity score: 0.94
#   Standard deviation: 5.200000e-04
#   Detected clusters: 0
#   Recommendation: Grid is nearly uniform - standard methods work well
# ========================================
...
```

### Working with Non-uniform Grids

```bash
# First, analyze your grid
./smooth -g nonuniform_data.txt

# Based on the grid uniformity (CV value), the program will automatically select:
# - Average coefficient method for nearly uniform grids (CV < 0.15)
# - Local spacing method for non-uniform grids (CV ≥ 0.15)

# After grid analysis, apply smoothing with automatic parameter selection
./smooth -m 2 -l auto nonuniform_data.txt

# For highly non-uniform grids (CV > 0.2), you may see:
# "WARNING: Highly non-uniform grid detected!"
# "GCV trace approximation may be less accurate."
# In this case, try manual λ or check results visually.
```

### Unix Filter Examples

The program can be seamlessly integrated into Unix pipelines:

```bash
# Simple pipe from cat
cat noisy_data.txt | ./smooth -m 1 -n 5 -p 2 > smoothed.txt

# Extract columns, smooth, and plot
awk '{print $1, $3}' experiment.dat | ./smooth -m 2 -l auto | gnuplot -p plot.gp

# Process multiple files
for f in data_*.txt; do
  cat "$f" | ./smooth -m 3 -f 0.15 > "smooth_$f"
done

# Filter out comments, smooth, extract columns
grep -v '^#' raw.txt | ./smooth -m 2 -l 0.01 | awk '{print $1, $2}' > final.txt

# Combine with other tools
./generate_data | ./smooth -m 1 -n 7 -p 3 | ./analyze_results

# Use in complex pipeline
curl https://example.com/data.txt | \
  grep -v '^#' | \
  ./smooth -m 2 -l auto | \
  awk '{if($2>threshold) print}' | \
  sort -k2 -n > filtered_smooth.txt

# Standard input/output redirection
./smooth -m 0 -n 5 -p 2 < input.dat > output.dat 2> errors.log

# Combine smoothing methods (not recommended, just for demo)
cat data.txt | ./smooth -m 2 -l 0.1 | ./smooth -m 1 -n 5 -p 2
```

### Typical Workflow

1. **Quick data exploration with grid analysis:**
   ```bash
   # Check grid uniformity only (program exits after analysis)
   ./smooth -g data.txt

   # Review output for:
   # - Grid uniformity (CV, ratio)
   # - Method recommendations
   ```

2. **Choose method based on grid:**
   ```bash
   # For uniform grids (CV < 0.05):
   ./smooth -m 1 -n 9 -p 3 -d data.txt

   # For non-uniform grids:
   ./smooth -m 2 -l auto -d data.txt
   ```

3. **Refine λ if needed:**
   ```bash
   # If automatic λ gives over-smoothing:
   ./smooth -m 2 -l 0.01 -d data.txt

   # If under-smoothing:
   ./smooth -m 2 -l 1.0 -d data.txt
   ```

4. **For publication graphics:**
   ```bash
   # Final smoothing with derivatives
   ./smooth -m 2 -l auto -d data.txt > publication_data.txt

   # Check functional balance in output comments
   # Ideal: both terms contribute 30-70%
   ```

---

## Grid Analysis Module

The `grid_analysis` module provides comprehensive analysis of input data and helps optimize smoothing parameters.

**Architecture:** Grid analysis is performed once at program startup (after data loading) and the results are shared across all smoothing methods. This eliminates redundant computation while ensuring all methods have access to consistent grid uniformity information. Methods that require uniform grids (Savitzky-Golay, Butterworth) receive pre-computed analysis results and can immediately reject unsuitable data with detailed recommendations.

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
    double *spacings;       // Array of spacings (optional)
    int n_points;           // Number of points
    int n_intervals;        // Number of intervals (n-1)
} GridAnalysis;
```

### Example Analysis Output

```
# ========================================
# GRID UNIFORMITY ANALYSIS
# ========================================
# Grid uniformity analysis:
#   n = 1000 points
#   h_min = 1.000000e-02, h_max = 1.000000e-01, h_avg = 5.500000e-02
#   CV = 0.450
#   Grid type: NON-UNIFORM
#   Uniformity score: 0.35
#   Standard deviation: 2.475000e-02
#   Detected clusters: 2
#   Recommendation: High non-uniformity - adaptive methods recommended
# WARNING: Significant spacing variation detected: CV = 0.45
# Adaptive methods may improve results.
#
# WARNING: 2 abrupt spacing changes detected (possible data clustering).
# Standard methods may over-smooth clustered regions.
# ========================================
```

### Grid Uniformity Thresholds

```
CV <= 0.01:  Uniform grid (is_uniform flag = 1) - all methods work optimally
             -> POLYFIT, SAVGOL, TIKHONOV all excellent

CV < 0.05:   Nearly uniform - SAVGOL works with warning
             -> SAVGOL may show warning but works
             -> POLYFIT and TIKHONOV work fine

0.05 ≤ CV < 0.15: Moderately non-uniform
             -> SAVGOL rejected automatically
             -> POLYFIT usable
             -> TIKHONOV uses average coefficient method

0.15 ≤ CV < 0.20: Non-uniform
             -> SAVGOL rejected
             -> POLYFIT usable with caution
             -> TIKHONOV uses local spacing method

CV ≥ 0.20:   Highly non-uniform
             -> SAVGOL rejected
             -> POLYFIT with caution
             -> TIKHONOV uses local spacing method (warning issued)

TIKHONOV METHOD SELECTION:
CV < 0.15:  Uses average coefficient method (robust, efficient)
CV ≥ 0.15:  Uses local spacing method (more accurate for non-uniform grids)
```

---

## Compilation and Installation

### Requirements

**Runtime:**
- C compiler (gcc, clang)
- LAPACK and BLAS libraries
- Make (optional, but recommended)

**Development/Testing:**
- Unity testing framework (included in `tests/` directory)
- Valgrind (optional, for memory leak detection)

### Compilation using Make

```bash
# Standard compilation
make

# Debug build
make debug

# Run unit tests
make test

# Run tests with Valgrind (memory leak detection)
make test-valgrind

# Clean build artifacts
make clean

# Clean test artifacts
make test-clean

# Install to user's home directory
make install-user

# Install to system (requires root)
make install

# Show all available targets
make help
```

### Manual Compilation

```bash
# Standard compilation
gcc -o smooth smooth.c polyfit.c savgol.c tikhonov.c butterworth.c \
    grid_analysis.c decomment.c -llapack -lblas -lm -O2

# With warnings
gcc -Wall -Wextra -pedantic -o smooth smooth.c polyfit.c savgol.c \
    tikhonov.c butterworth.c grid_analysis.c decomment.c -llapack -lblas -lm -O2
```

### File Structure

```
smooth/
|--- smooth.c           # Main program
|--- polyfit.c/h        # Polynomial fitting module
|--- savgol.c/h         # Savitzky-Golay module
|--- tikhonov.c/h       # Tikhonov regularization module
|--- butterworth.c/h    # Butterworth filter module
|--- grid_analysis.c/h  # Grid analysis module
|--- decomment.c/h      # Comment removal utility
|--- timestamp.c/h      # Timestamp parsing module
|--- revision.h         # Program version
|--- Makefile           # Build system with test targets
|--- README.md          # This documentation
+--- tests/             # Unit testing framework (Unity)
    |--- unity.c/h                # Unity testing framework
    |--- unity_internals.h        # Unity internals
    |--- test_main.c              # Test runner (50+ tests)
    |--- test_grid_analysis.c     # Grid analysis tests (7 tests)
    |--- test_polyfit.c           # Polyfit module tests (18 tests)
    |--- test_savgol.c            # Savgol module tests (16 tests)
    +--- test_timestamp.c         # Timestamp module tests (15 tests)
```

---

## Conclusion

The `smooth` program provides four complementary smoothing methods in a modular architecture with advanced input data analysis and comprehensive testing:

- **POLYFIT** - local polynomial approximation using least squares method
- **SAVGOL** - optimal linear filter with pre-computed coefficients (uniform grids only)
- **TIKHONOV** - global variational method with hybrid automatic discretization
- **BUTTERWORTH** - digital low-pass filter with zero-phase filtfilt
- **GRID_ANALYSIS** - automatic analysis and method recommendation

### Key Features

**Robust Testing Infrastructure:**
- 50+ unit tests using Unity testing framework
- AAA pattern (Arrange-Act-Assert) for all tests
- Memory leak detection via Valgrind integration
- Edge case coverage for robust production use
- Test modules: grid_analysis (7 tests), polyfit (18 tests), savgol (16 tests), timestamp (15 tests)

**Advanced Capabilities:**
- Centralized grid analysis performed once at startup
- Unix filter support for integration into pipelines
- Timestamp mode for RFC3339 time-series data
- Automatic parameter selection (GCV for Tikhonov)
- Scipy-compatible algorithms (Butterworth filtfilt, lfilter_zi)
- Zero-phase filtering for Butterworth method
- Hybrid discretization for Tikhonov on non-uniform grids

### When to Use Each Method

**Quick Decision Tree:**

```
Is your grid uniform (CV < 0.05)?
|-- YES: Multiple good options:
|   |-- SAVGOL: Best for polynomial signals with derivatives
|   |   smooth -m 1 -n 9 -p 3 -d data.txt
|   |-- BUTTERWORTH: Best for frequency-domain interpretation
|   |   smooth -m 3 -f 0.15 data.txt
|   +-- TIKHONOV: Universal choice with auto parameters
|       smooth -m 2 -l auto -d data.txt
|
+-- NO (non-uniform): Use TIKHONOV for correct handling
    smooth -m 2 -l auto -d data.txt

Need frequency-domain control?
+-- Use BUTTERWORTH (requires uniform grid)
    smooth -m 3 -f 0.15 data.txt

Need local adaptability?
+-- Use POLYFIT regardless of grid
    smooth -m 0 -n 7 -p 2 -d data.txt

Not sure?
+-- Use TIKHONOV with automatic λ - safest choice!
    smooth -m 2 -l auto -d data.txt
```

Each method has a strong mathematical foundation and is optimized for specific data types. The program provides automatic guidance on method selection and parameters, with extensive diagnostics to ensure correct usage.

### Best Practices

**For Users:**
1. **Always check grid first:** Use `-g` flag to understand your data
2. **Start with automatic:** Use `-l auto` for Tikhonov, let GCV find optimal λ
3. **Check functional balance:** Look for 30-70% split between data and regularization terms
4. **Iterate if needed:** Adjust λ manually if automatic selection doesn't satisfy requirements
5. **Use derivatives wisely:** Add `-d` only when needed - cleaner output without it
6. **Understand the trade-off:** More smoothing (larger λ) = more noise reduction but less detail

**For Developers:**
7. **Run tests before committing:** Always run `make test` to verify no regressions
8. **Check for memory leaks:** Use `make test-valgrind` to ensure clean memory management
9. **Write tests for new features:** Follow AAA pattern (Arrange-Act-Assert) used in existing tests
10. **Maintain test coverage:** Add edge cases and ensure numerical tolerances are appropriate

---

**Document revision:** 2025-12-03
**Program version:** smooth v5.9.2
**Dependencies:** LAPACK, BLAS
**Testing framework:** Unity (included in tests/)
**License:** MIT License
