# Methods for Smoothing Experimental Data in the smooth Program

**Technical Documentation**
Version 5.11.0 | February 7, 2026

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

$$y_{\text{obs}}(x_i) = y_{\text{true}}(x_i) + \varepsilon_i$$

The goal of smoothing is to estimate $y_{\text{true}}$ while suppressing $\varepsilon_i$ and preserving physically relevant signal properties.

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
- `-f auto` - automatic cutoff selection (butterworth, currently returns 0.2)
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

**Problem:** For each point $x_i$, we fit a polynomial of degree $p$ to the surrounding $n$ points:

$$P(x) = a_0 + a_1(x-x_i) + a_2(x-x_i)^2 + \cdots + a_p(x-x_i)^p$$

**Optimization criterion:**

$$\min \sum_{j \in [i-n/2,\, i+n/2]} \left[ y_j - P(x_j) \right]^2$$

### Least Squares Solution via SVD

The polynomial coefficients are found by solving an **overdetermined linear system** using **Singular Value Decomposition (SVD)**:

$$V \cdot \mathbf{a} = \mathbf{y}_{\text{window}}$$

where $V$ is the Vandermonde matrix with $V_{j,k} = (x_j - x_i)^k$:

$$V = \begin{pmatrix} 1 & (x_0-x_i) & (x_0-x_i)^2 & \cdots & (x_0-x_i)^p \\ 1 & (x_1-x_i) & (x_1-x_i)^2 & \cdots & (x_1-x_i)^p \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 1 & (x_n-x_i) & (x_n-x_i)^2 & \cdots & (x_n-x_i)^p \end{pmatrix}$$

$$\mathbf{a} = [a_0, a_1, \ldots, a_p]^T \quad\text{(polynomial coefficients)}, \qquad \mathbf{y}_{\text{window}} = [y_{i-n/2}, \ldots, y_{i+n/2}]^T$$

**Why SVD instead of Normal Equations?**

The implementation uses LAPACK's `dgelss` (SVD decomposition) rather than forming normal equations $(V^T V)\mathbf{a} = V^T \mathbf{y}$:

1. **Numerical stability:** SVD avoids squaring the condition number ($\kappa(V^T V) = \kappa(V)^2$)
2. **Automatic regularization:** Singular values below $\text{rcond} \cdot \sigma_{\max}$ are truncated
3. **Rank detection:** Provides effective rank for diagnosing ill-conditioning
4. **Robustness:** Handles high polynomial degrees (p > 6) more reliably

**SVD truncation parameter:** $\text{rcond} = 10^{-10}$
- Singular values $\sigma_i < 10^{-10} \cdot \sigma_{\max}$ are treated as zero
- Provides implicit Tikhonov-style regularization
- Conservative threshold ensures stability without over-regularization

### Derivative Computation

Derivatives are computed analytically from polynomial coefficients:

$$f(x_i) = a_0, \qquad f'(x_i) = a_1, \qquad f''(x_i) = 2a_2$$

### Edge Handling

At edges, asymmetric windows are used with extrapolation of the fitted polynomial:

$$f(x_k) = \sum_{m=0}^{p} a_m \cdot (x_k - x_{n/2})^m$$

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
- On first window, reports condition number: $\kappa = \sigma_{\max} / \sigma_{\min}$
- If $\kappa > 10^8$, issues warning about potential numerical issues
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

**Fundamental principle:** For given parameters (window size, polynomial degree, derivative order), there exist universal coefficients $c_k$ such that:

$$f^{(d)}(x_i) = \sum_{k=-n_L}^{n_R} c_k \cdot y_{i+k}$$

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

$$CV = \frac{\sigma(h)}{h_{\text{avg}}}$$

| CV range | Decision |
|----------|----------|
| $CV > 0.05$ | REJECT — Grid too non-uniform for SG |
| $CV > 0.01$ | WARNING — Nearly uniform, proceed with caution |
| $CV \le 0.01$ | OK — Grid sufficiently uniform |

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

$$\hat{y}_0 = c_{-n_L} \cdot y_{-n_L} + \cdots + c_0 \cdot y_0 + \cdots + c_{n_R} \cdot y_{n_R}$$

where the coefficients $c_k$ are "undetermined" and must satisfy the condition that the filter exactly reproduces polynomials up to degree $p$.

**The key insight:** For a given window configuration and polynomial degree, these coefficients can be determined once and applied universally - but only on uniform grids!

### Coefficient Derivation

Coefficients are derived from the condition that the filter must exactly reproduce polynomials up to degree $p$.

**Moment conditions:**

$$\sum_{j=-n_L}^{n_R} c_j \cdot j^m = \delta_{m,d} \cdot d! \qquad \text{for } m = 0, 1, \ldots, p$$

where $\delta_{m,d}$ is the Kronecker delta, $d$ is the derivative order, and $d!$ is factorial.

This leads to a system of linear equations where the unknowns are the filter coefficients $c_j$.

### Matrix Formulation

The coefficients are found by solving a **normal equations system** (not a Vandermonde system):

$$A \cdot \boldsymbol{\beta} = \mathbf{b}$$

where $A$ is a symmetric $(p+1) \times (p+1)$ moment matrix with $A_{i,j} = \sum_{k=-n_L}^{n_R} k^{i+j}$:

$$A = \begin{pmatrix} \sum k^0 & \sum k^1 & \sum k^2 & \cdots & \sum k^p \\ \sum k^1 & \sum k^2 & \sum k^3 & \cdots & \sum k^{p+1} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ \sum k^p & \sum k^{p+1} & \cdots & \cdots & \sum k^{2p} \end{pmatrix}$$

and the right-hand side vector: $b_j = \delta_{j,d} \cdot d!$

This results in a symmetric positive definite $(p+1) \times (p+1)$ matrix. The filter coefficients are then:

$$c_k = \sum_{j=0}^{p} \beta_j \cdot k^j \qquad \text{for } k = -n_L, \ldots, n_R$$

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

$$\frac{dy}{dx}\bigg|_{\text{physical}} = \frac{1}{h_{\text{avg}}} \cdot \frac{dy}{dx}\bigg|_{\text{normalized}}$$

where $h_{\text{avg}}$ is the average grid spacing.

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

$$J[u] = \underbrace{\int (y(x) - u(x))^2 \, dx}_{\text{Data fidelity}} + \lambda \underbrace{\int (u''(x))^2 \, dx}_{\text{Smoothness penalty}}$$

**Discrete formulation:**

$$J[\mathbf{u}] = \|\mathbf{y} - \mathbf{u}\|^2 + \lambda \|D^2 \mathbf{u}\|^2$$

where:
- $\|\mathbf{y} - \mathbf{u}\|^2 = \sum (y_i - u_i)^2$ is the **data fidelity term**
- $\|D^2 \mathbf{u}\|^2 = \sum (D^2 u_i)^2$ is the **regularization term** (smoothness penalty)
- $\lambda$ is the **regularization parameter** controlling the balance
- $D^2$ is the discrete second derivative operator

### The Regularization Parameter λ

The parameter λ is the **heart of Tikhonov regularization** - it controls the balance between fitting the data and smoothing the result.

#### Physical Interpretation

| $\lambda$ | Effect | Functional |
|-----------|--------|------------|
| $\lambda = 0$ | No smoothing, $u = y$ | $J[u] = \|\mathbf{y} - \mathbf{u}\|^2$ only |
| $\lambda \to \infty$ | Maximum smoothing, $u \to$ straight line | $J[u] \approx \lambda\|D^2 \mathbf{u}\|^2$ dominates |
| $\lambda$ optimal | Balanced data fit and smoothness | Both terms contribute meaningfully |

#### Mathematical Role

The minimization of $J[\mathbf{u}]$ leads to:

$$(I + \lambda (D^2)^T W D^2) \, \mathbf{u} = \mathbf{y}$$

**Effect of $\lambda$ on the solution:**
- **Small $\lambda$ (< 0.01):** Matrix $\approx I$ → solution $\mathbf{u} \approx \mathbf{y}$ (minimal smoothing)
- **Large $\lambda$ (> 1.0):** Matrix $\approx \lambda (D^2)^T W D^2$ → strong curvature penalty (heavy smoothing)
- **Optimal λ:** Matrix components balanced → noise removed, signal preserved

#### Frequency Domain Interpretation

In Fourier space, Tikhonov acts as a low-pass filter:

$$\hat{H}(\omega) = \frac{1}{1 + \lambda \omega^4}$$

where $\omega$ is spatial frequency.

**Effect:**
- **Low frequencies (slow variations):** $\hat{H} \approx 1$ → preserved
- **High frequencies (noise, rapid variations):** $\hat{H} \approx 1/(\lambda \omega^4)$ → attenuated
- **Cutoff frequency:** $\omega_c \sim \lambda^{-1/4}$

**This means:** Larger $\lambda$ → lower cutoff → more aggressive low-pass filtering → smoother result. Smaller $\lambda$ → higher cutoff → less filtering → result closer to data.

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

$$\text{Effective strength} \sim \frac{\lambda}{h_{\text{avg}}^4}$$

Same $\lambda$ on finer grid → weaker smoothing; same $\lambda$ on coarser grid → stronger smoothing.

For dimensional consistency, $\lambda$ has units $[\text{Length}^4]$ (since it scales the squared second derivative).

### Second Derivative Penalty: (D²)ᵀWD² Gram Matrix

The regularization term $\|D^2 \mathbf{u}\|^2$ is discretized using the **Gram matrix** of the second derivative operator $D^2$. The operator $D^2$ is an $(n-2) \times n$ matrix — it has rows only for interior points ($k = 1, \ldots, n-2$). This means natural boundary conditions ($D^2 u = 0$ at endpoints) are **implicit**: there are simply no $D^2$ rows for boundary points.

The penalty matrix is constructed as a sum of rank-1 contributions:

$$(D^2)^T W D^2 = \sum_{k=1}^{n-2} w_k \cdot \mathbf{d}_k^T \cdot \mathbf{d}_k$$

where $\mathbf{d}_k$ is the $k$-th row of $D^2$ (a 3-element stencil at positions $k{-}1, k, k{+}1$) and $w_k$ is the integration weight. Since each $\mathbf{d}_k$ touches 3 consecutive points, the Gram matrix is **pentadiagonal** (bandwidth $kd = 2$).

Automatic selection between two discretization schemes based on grid uniformity.

#### Grid Uniformity Detection

$$CV = \frac{h_{\text{std}}}{h_{\text{avg}}}$$

| CV range | Discretization |
|----------|---------------|
| $CV < 0.15$ | Average Coefficient Method (nearly uniform) |
| $CV \ge 0.15$ | Local Spacing Method (non-uniform) |

#### Method 1: Average Coefficient (for CV < 0.15)

Used for uniform and mildly non-uniform grids. More robust numerically.

**$D^2$ stencil:** Each interior row $k$ uses uniform spacing $h_{\text{avg}}$:

$$\mathbf{d}_k = \frac{1}{h_{\text{avg}}^2} [1, \; -2, \; 1] \qquad \text{at positions } (k{-}1, \; k, \; k{+}1)$$

**Gram matrix construction:**

$$(D^2)^T D^2 = \sum_{k=1}^{n-2} \mathbf{d}_k^T \mathbf{d}_k = \frac{1}{h_{\text{avg}}^4} \sum_{k=1}^{n-2} \mathbf{d}_k^T \mathbf{d}_k$$

The accumulated result is a pentadiagonal matrix with the classical 4th-difference stencil $[1, -4, 6, -4, 1] / h_{\text{avg}}^4$.

**Example** for $n = 6$ points with $c = \lambda / h^4$:

$$A = I + \lambda (D^2)^T D^2 = \begin{pmatrix} 1+c & -2c & c & 0 & 0 & 0 \\ -2c & 1+5c & -4c & c & 0 & 0 \\ c & -4c & 1+6c & -4c & c & 0 \\ 0 & c & -4c & 1+6c & -4c & c \\ 0 & 0 & c & -4c & 1+5c & -2c \\ 0 & 0 & 0 & c & -2c & 1+c \end{pmatrix}$$

Note: Boundary points (first/last rows) receive less regularization because fewer $D^2$ stencils overlap them.

#### Method 2: Local Spacing (for CV ≥ 0.15)

Used for highly non-uniform grids. More accurate for variable spacing.

**$D^2$ stencil:** Each interior row $k$ uses local spacings $h_l = x_k - x_{k-1}$, $h_r = x_{k+1} - x_k$:

$$\mathbf{d}_k = [a_k, \; b_k, \; c_k] \qquad \text{at positions } (k{-}1, \; k, \; k{+}1)$$

where:

$$a_k = \frac{2}{(h_l + h_r) \cdot h_l}, \qquad b_k = \frac{-2}{h_l \cdot h_r}, \qquad c_k = \frac{2}{(h_l + h_r) \cdot h_r}$$

This is the **correct second derivative formula** for non-uniform grids derived from Taylor expansion.

**Integration weight:** $w_k = (h_l + h_r) / 2$

**Gram matrix construction:**

$$(D^2)^T W D^2 = \sum_{k=1}^{n-2} w_k \cdot \mathbf{d}_k^T \mathbf{d}_k$$

For each interior point $k$, the rank-1 outer product $w_k \cdot \mathbf{d}_k^T \mathbf{d}_k$ is accumulated into the pentadiagonal matrix (upper triangle only): diagonal, 1st superdiagonal, and 2nd superdiagonal.

**The resulting matrix is:**
- Symmetric (Gram matrix structure guarantees this automatically)
- Positive definite
- Pentadiagonal (bandwidth $kd = 2$)

#### Natural Boundary Conditions

Natural boundary conditions ($D^2 u = 0$ at endpoints) are **implicit** in the Gram matrix construction:

- $D^2$ is an $(n{-}2) \times n$ matrix with rows only for interior points $k = 1, \ldots, n{-}2$
- No boundary rows exist in $D^2$, so no explicit boundary penalty is applied
- Boundary points ($i = 0$, $i = n{-}1$) are regularized only through their participation in nearby interior stencils
- This approach avoids the need for explicit boundary formulas and produces cleaner edge behavior

**Boundary behavior:**
- Boundary points receive less regularization pressure than interior points
- For very small λ, boundary points may track the data more closely than interior points
- This is generally desirable: edges are less constrained, avoiding artificial boundary effects

### Functional Computation

The actual value of the minimized functional is computed for diagnostic purposes:

**Data term:**

$$\|\mathbf{y} - \mathbf{u}\|^2 = \sum_{i=0}^{n-1} (y_i - u_i)^2$$

**Regularization term (interior points only):**

Consistent with the Gram matrix construction, the regularization term sums only over interior points (natural BCs are implicit — no boundary terms):

For **average coefficient method:**

$$\|D^2 \mathbf{u}\|^2 = \sum_{i=1}^{n-2} \left[\frac{u_{i-1} - 2u_i + u_{i+1}}{h_{\text{avg}}^2}\right]^2$$

For **local spacing method:**

$$\|D^2 \mathbf{u}\|_W^2 = \sum_{i=1}^{n-2} (D^2 u_i)^2 \cdot \frac{h_l + h_r}{2}$$

where $D^2 u_i = \frac{2}{h_l + h_r} \left[\frac{u_{i-1}}{h_l} - u_i \left(\frac{1}{h_l} + \frac{1}{h_r}\right) + \frac{u_{i+1}}{h_r}\right]$

**Total functional:**

$$J[\mathbf{u}] = \|\mathbf{y} - \mathbf{u}\|^2 + \lambda \|D^2 \mathbf{u}\|^2$$

### Variational Approach

The minimum of functional $J[\mathbf{u}]$ satisfies the Euler-Lagrange equation:

$$\frac{\partial J}{\partial u_i} = 0 \implies -2(y_i - u_i) + 2\lambda \left((D^2)^T W D^2 \, \mathbf{u}\right)_i = 0$$

which leads to the linear system:

$$(I + \lambda (D^2)^T W D^2) \, \mathbf{u} = \mathbf{y}$$

where $(D^2)^T W D^2$ is the Gram matrix of the second derivative operator $D^2$.

### Matrix Representation

The linear system $(I + \lambda (D^2)^T W D^2) \, \mathbf{u} = \mathbf{y}$ has matrix $A = I + \lambda (D^2)^T W D^2$ with structure:

**Properties:**
- Symmetric (Gram matrix structure guarantees this)
- Positive definite
- Pentadiagonal (banded with bandwidth $kd = 2$)

**General pentadiagonal form:**

$$A = \begin{pmatrix} d_0 & e_0 & f_0 & 0 & \cdots & 0 \\ e_0 & d_1 & e_1 & f_1 & \cdots & 0 \\ f_0 & e_1 & d_2 & e_2 & \cdots & 0 \\ 0 & f_1 & e_2 & d_3 & \cdots & 0 \\ \vdots & & & & \ddots & \vdots \\ 0 & 0 & 0 & 0 & e_{n-2} & d_{n-1} \end{pmatrix}$$

where $d_i$ = diagonal, $e_i$ = 1st off-diagonal, $f_i$ = 2nd off-diagonal.

The pentadiagonal structure arises because each row of $D^2$ touches 3 consecutive points ($k{-}1, k, k{+}1$), so the Gram matrix $(D^2)^T D^2$ couples points up to 2 positions apart.

This structure allows efficient solution using LAPACK's banded solver `dpbsv`.

### Generalized Cross Validation (GCV)

For automatic $\lambda$ selection (`-l auto`), we minimize the GCV criterion:

$$\text{GCV}(\lambda) = \frac{n \cdot \text{RSS}(\lambda)}{(n - \text{tr}(H_\lambda))^2}$$

where:
- $\text{RSS}(\lambda) = \|\mathbf{y} - \mathbf{u}_\lambda\|^2$ is the residual sum of squares
- $H_\lambda = (I + \lambda (D^2)^T W D^2)^{-1}$ is the influence matrix (smoother matrix)
- $\text{tr}(H_\lambda)$ is the trace (effective number of parameters)

**Interpretation:**
- $\text{tr}(H_\lambda)$ measures model complexity (degrees of freedom)
- Small $\lambda$: $\text{tr}(H) \approx n$ (interpolation, overfitting)
- Large $\lambda$: $\text{tr}(H) \to 2$ (linear fit, underfitting — $D^2$ null space is constants + linear)
- Optimal $\lambda$: minimizes prediction error

**Trace estimation using eigenvalues:**

For uniform grids, the eigenvalues of $(D^2)^T D^2$ are the squares of the eigenvalues of $D^T D$ (first-derivative operator):

$$\text{tr}(H_\lambda) \approx 2 + \sum_{k=1}^{n-2} \frac{1}{1 + \lambda \mu_k}$$

where the eigenvalues of $(D^2)^T D^2$ are:

$$\theta_k = \frac{\pi k}{n}, \qquad \mu_k = \left(\frac{4 \sin^2(\theta_k / 2)}{h^2}\right)^2$$

The null space of $D^2$ is 2-dimensional (constants and linear functions), so trace starts at 2.0 (these two modes are unpenalized: $1/(1+0) = 1$ each).

**Note:** This approximation is exact for uniform grids but approximate for non-uniform grids. For highly non-uniform grids (CV > 0.2), the program issues a warning.

#### Enhanced GCV

**Over-fitting penalty:** If $\text{tr}(H)/n > 0.7$:

$$\text{GCV}_{\text{modified}} = \text{GCV} \cdot \exp\!\left(10 \cdot \left(\frac{\text{tr}(H)}{n} - 0.7\right)\right)$$

This exponential penalty prevents selection of too-small $\lambda$ that would lead to overfitting.

**L-curve backup (for n > 20000):**

For very large datasets, GCV trace approximation may be inaccurate. The program also computes the L-curve (plot of $\|D^2 \mathbf{u}\|^2$ vs $\|\mathbf{y} - \mathbf{u}\|^2$) and finds the corner point with maximum curvature:

$$\kappa = \frac{|x' y'' - y' x''|}{(x'^2 + y'^2)^{3/2}}$$

where $x = \log\|\mathbf{y} - \mathbf{u}\|^2$ and $y = \log\|D^2 \mathbf{u}\|^2$.

If GCV and L-curve disagree significantly, the program uses the more conservative (larger) $\lambda$.

### Efficient Implementation

The program uses LAPACK routine `dpbsv` for solving symmetric positive definite banded systems:

```c
// Banded matrix storage (LAPACK column-major format)
// For pentadiagonal symmetric matrix A with bandwidth kd=2:
//
//     [ AB[0,j] ]  = 2nd superdiagonal elements  (a[i,i+2])
//     [ AB[1,j] ]  = 1st superdiagonal elements  (a[i,i+1])
//     [ AB[2,j] ]  = diagonal elements            (a[i,i])
//
// Storage layout for pentadiagonal matrix:
//
//         [ *    *    a02  a13  a24  ... ]   <- row 0 (2nd superdiagonal)
//         [ *    a01  a12  a23  a34  ... ]   <- row 1 (1st superdiagonal)
//    AB = [ a00  a11  a22  a33  a44  ... ]   <- row 2 (diagonal)
//
// System solution (kd=2, ldab=3)
dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, b, &n, &info);
```

**Complexity:**
- Memory: O(n) for banded storage (3n elements)
- Time: O(n) for factorization and back-substitution

This is **optimal** for pentadiagonal systems.

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
TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda,
                                GridAnalysis *grid_info);

// Automatic λ selection
double find_optimal_lambda_gcv(double *x, double *y, int n, GridAnalysis *grid_info);

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

$$|H(j\omega)|^2 = \frac{1}{1 + (\omega / \omega_c)^{2N}}$$

where:
- $N$ = filter order (4 in our implementation)
- $\omega_c$ = cutoff frequency (3dB point)
- $\omega$ = frequency

**Key Properties:**
- **Maximally flat passband:** No ripples for ω < ωc
- **Monotonic rolloff:** Smooth transition from passband to stopband
- **-3dB at cutoff:** |H(jωc)| = 1/√2 ≈ 0.707
- **Rolloff rate:** -20N dB/decade (for N=4: -80 dB/decade)

### Digital Implementation

The smooth program implements a **4th-order digital Butterworth low-pass filter** using the following algorithm:

**Step 1: Pole Calculation**

Butterworth poles lie on unit circle in s-domain at angles:

$$\theta_k = \frac{\pi}{2} + \frac{\pi(2k+1)}{2N}, \qquad k = 0, 1, \ldots, N{-}1$$

For $N = 4$:

$$s_{\text{poles}}[k] = e^{j\theta_k} \qquad \text{where } \theta \in \left\{\frac{5\pi}{8}, \frac{7\pi}{8}, \frac{9\pi}{8}, \frac{11\pi}{8}\right\}$$

**Step 2: Frequency Scaling**

Scale poles by prewarped cutoff frequency:

$$\omega_c = \tan(\pi f_c / 2), \qquad s_{\text{scaled}} = \omega_c \cdot s_{\text{poles}}$$

**Prewarping correction:** The bilinear transform introduces frequency warping. The factor $\tan(\pi f_c / 2)$ compensates for this, ensuring the digital filter's cutoff matches the desired normalized frequency $f_c$.

**Step 3: Bilinear Transform**

Convert analog poles to digital domain:

$$z_{\text{poles}} = \frac{2 + s_{\text{scaled}}}{2 - s_{\text{scaled}}}$$

The bilinear transformation maps:
- Left half of s-plane → inside unit circle in z-plane
- jω axis → unit circle in z-plane
- Preserves stability

**Step 4: Biquad Cascade**

Form two 2nd-order sections (biquads) from conjugate pole pairs:

$$H(z) = H_1(z) \cdot H_2(z), \qquad H_i(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}$$

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
- **Effective order:** $2N = 8$ (squared magnitude response)
- **Steeper rolloff:** $|H_{\text{eff}}(j\omega)|^2 = |H(j\omega)|^4$

### Initial Conditions (Biquad IC)

To minimize edge transients, we compute initial filter state for each **biquad section** using an **analytical solution**:

**Problem:** For each 2nd-order biquad section, find initial state $\mathbf{z}_i$ such that for constant input $x = c$:

$$\mathbf{z}_i = A \cdot \mathbf{z}_i + B \cdot c$$

This ensures each biquad starts in steady-state, eliminating startup transients.

**Solution for 2nd-order biquad:** Solve the $2 \times 2$ linear system analytically:

$$(I - A) \cdot \mathbf{z}_i = \mathbf{B}$$

where for Transposed Direct Form II:

$$I - A = \begin{pmatrix} 1+a_1 & -1 \\ a_2 & 1 \end{pmatrix}, \qquad \mathbf{B} = \begin{pmatrix} b_1 - a_1 b_0 \\ b_2 - a_2 b_0 \end{pmatrix}$$

**Analytical solution using Cramer's rule:**

$$\det(I - A) = 1 + a_1 + a_2$$

$$z_i[0] = \frac{B_0 + B_1}{\det}, \qquad z_i[1] = \frac{-a_2 B_0 + (1 + a_1) B_1}{\det}$$

**Implementation advantages:**
- **No LAPACK dependency** for initial conditions (purely analytical)
- **Per-biquad computation** - simple 2×2 systems instead of 4×4
- **Numerically robust** - direct formulas avoid iterative solvers
- **Efficient** - closed-form solution, no matrix decomposition needed

This approach maintains numerical stability even for challenging filter coefficients (e.g., very low cutoff frequencies where coefficients can be extremely small ~10⁻⁸), while being simpler and faster than the previous LAPACK-based companion matrix approach.

### Normalized Cutoff Frequency

The cutoff frequency `fc` is **normalized** to the sampling rate and is the **most important parameter** for Butterworth filtering.

**Simple explanation:**
- `fc` controls how much smoothing you get
- **Smaller fc (e.g., 0.05)** → heavy smoothing, only very slow trends preserved
- **Larger fc (e.g., 0.30)** → light smoothing, more detail preserved
- Valid range: `0 < fc < 1.0` (Nyquist limit)

**Technical details:**

$$f_c = \frac{f_{\text{cutoff}}}{f_{\text{Nyquist}}} = \frac{f_{\text{cutoff}}}{f_{\text{sample}} / 2}, \qquad f_{\text{sample}} = \frac{1}{h_{\text{avg}}}$$

**Nyquist Constraint:** $0 < f_c < 1.0$
- $f_c = 1.0$ corresponds to Nyquist frequency ($f_{\text{sample}}/2$) - maximum possible
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

### Biquad Cascade Implementation

The 4th-order Butterworth filter is implemented as a **cascade of 2 second-order sections (biquads)** for superior numerical stability.

**Why biquad cascade?**
- **Numerical stability:** Each 2nd-order section has well-conditioned coefficient magnitudes
- **Reduced quantization errors:** Coefficients stay in reasonable ranges (~0.1 to 10)
- **Industry standard:** Used in professional DSP applications
- **Modular:** Each biquad can be analyzed and tested independently

**Filter structure:**
```
x → [Biquad 1] → [Biquad 2] → y

Each biquad processes conjugate pole pair from Butterworth prototype
```

Each biquad section uses **Transposed Direct Form II** (TDF-II) for optimal numerical properties:

$$y[n] = b_0 \cdot x[n] + z_0$$
$$z_0 = b_1 \cdot x[n] - a_1 \cdot y[n] + z_1$$
$$z_1 = b_2 \cdot x[n] - a_2 \cdot y[n]$$

where $\mathbf{z}$ is the biquad state (2 elements per section), $[b_0, b_1, b_2]$ are numerator coefficients, and $[1, a_1, a_2]$ are denominator coefficients ($a_0$ normalized to 1).

**Complete 4th-order filtering:**
1. Apply first biquad to input → intermediate result
2. Apply second biquad to intermediate result → final output

This cascade approach is more robust than direct 4th-order implementation, especially for very low cutoff frequencies (fc < 0.05) where monolithic implementations can suffer from coefficient scaling issues.

### Modularized Implementation

```c
// butterworth.h

/* Biquad section structure (2nd-order IIR filter) */
typedef struct {
    double b[3];  // Numerator coefficients: [b0, b1, b2]
    double a[3];  // Denominator coefficients: [a0=1, a1, a2]
} BiquadSection;

/* Complete filter structure (cascade of 2 biquads) */
typedef struct {
    BiquadSection sections[2];  // Two 2nd-order sections for 4th order
} ButterworthCoeffs;

/* Result structure */
typedef struct {
    double *y_smooth;     // Smoothed values
    int n;                // Number of points
    int order;            // Filter order (BUTTERWORTH_ORDER = 4)
    double cutoff_freq;   // Normalized cutoff frequency (0 < fc < 1)
    double sample_rate;   // Effective sample rate from data spacing
} ButterworthResult;

// Main function
ButterworthResult* butterworth_filtfilt(const double *x, const double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        const GridAnalysis *grid_info);

// Automatic cutoff selection (currently returns 0.2)
double estimate_cutoff_frequency(const double *x, const double *y, int n);

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
- **Superior numerical stability** (biquad cascade implementation)
- **No LAPACK dependency** for initial conditions (analytical solution)
- **Classical DSP approach** with extensive literature and understanding
- **Predictable frequency-domain behavior** - easy to interpret cutoff frequency
- **No ringing** (unlike Chebyshev or elliptic filters)
- **Efficient implementation** - O(n) time complexity
- **Smooth monotonic rolloff** - natural attenuation curve
- **Robust for extreme cutoffs** (fc < 0.05 handled well by biquad cascade)

**Disadvantages:**
- **Requires uniform/nearly-uniform grid** (CV < 0.15 enforced, warning for CV > 0.05)
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
- Use -f auto (currently returns default fc = 0.2)
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

# Butterworth with automatic cutoff (currently returns 0.2)
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
    |--- test_main.c              # Test runner (102 tests)
    |--- test_grid_analysis.c     # Grid analysis tests (7 tests)
    |--- test_polyfit.c           # Polyfit module tests (21 tests)
    |--- test_savgol.c            # Savgol module tests (16 tests)
    |--- test_tikhonov.c          # Tikhonov module tests (26 tests)
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
- 102 unit tests using Unity testing framework
- AAA pattern (Arrange-Act-Assert) for all tests
- Memory leak detection via Valgrind integration
- Edge case coverage for robust production use
- Test modules: grid_analysis (7), polyfit (21), savgol (16), tikhonov (26), butterworth (17), timestamp (15)

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

**Document revision:** 2026-02-07
**Program version:** smooth v5.11.0
**Dependencies:** LAPACK, BLAS
**Testing framework:** Unity (included in tests/)
**License:** MIT License

**Recent changes (v5.11.0):**
- Tikhonov regularization corrected to use true 2nd-order penalty (D²)ᵀWD²
- Pentadiagonal Gram matrix (kd=2) replaces previous tridiagonal approximation
- Natural boundary conditions now implicit (no boundary rows in D²)
- GCV eigenvalues corrected: μ_k = (4sin²(θ/2)/h²)² with 2D null space
- 102 unit tests including 26 Tikhonov-specific tests
