# 4th-order Regge--Wheeler time-domain solver

Numerical solver for the **Regge--Wheeler equation**, describing perturbations of a Schwarzschild black hole and their associated quasi-normal mode ringdown.

The code evolves a scalar or gravitational wave in tortoise coordinates using:

- 4th-order Runge--Kutta time integration,
- 4th-order finite differences in space,
- smooth windowing of the Regge--Wheeler potential,
- 4th-order Lagrange interpolation for waveform extraction.

---

## Physics background

The Regge--Wheeler equation governs odd-parity gravitational perturbations of a Schwarzschild black hole of mass \(M\):

```math
\left(
-\frac{\partial^2}{\partial t^2}
+\frac{\partial^2}{\partial r_*^2}
-V_\ell(r)
\right)\Psi(t,r_*) = 0 .
```

Here, \(r_*\) is the tortoise coordinate and \(V_\ell(r)\) is the Regge--Wheeler potential.

The Schwarzschild metric in standard coordinates is

```math
ds^2 =
-\left(1-\frac{2M}{r}\right)dt^2
+\left(1-\frac{2M}{r}\right)^{-1}dr^2
+r^2 d\Omega^2 .
```

The event horizon is located at

```math
r = 2M .
```

---

## Coordinate system

### Tortoise coordinate

The tortoise coordinate \(r_*\), denoted by `s` in the code, is defined as

```math
r_* =
r + 2M\ln\left(\frac{r}{2M}-1\right).
```

This maps the exterior Schwarzschild region

```math
r \in (2M,\infty)
```

to

```math
r_* \in (-\infty,+\infty).
```

### Inverting the tortoise coordinate

Given \(r_*\), the physical Schwarzschild radius \(r\) can be recovered using the principal branch of the Lambert W function:

```math
r(r_*) =
2M\left[
1 + W_0\left(
\exp\left(\frac{r_*}{2M}-1\right)
\right)
\right].
```

This inversion is implemented in `r_of_s(s, M)`, using `scipy.special.lambertw` when available, and `mpmath.lambertw` as a fallback.

---

## Regge--Wheeler potential

The effective potential for gravitational perturbations with angular momentum number \(\ell\) is

```math
V_\ell(r) =
\left(1-\frac{2M}{r}\right)
\left[
\frac{\ell(\ell+1)}{r^2}
-\frac{6M}{r^3}
\right].
```

Its main contributions are:

- the redshift factor \((1-2M/r)\), which vanishes at the horizon;
- the centrifugal barrier \(\ell(\ell+1)/r^2\);
- the curvature correction \(-6M/r^3\), characteristic of gravitational perturbations.

For gravitational waves, the dominant mode is usually \(\ell=2\).

---

## Windowing

To suppress spurious reflections from the grid boundaries, the potential is multiplied by a smooth window function \(W(r)\):

```math
V_{\rm eff}(r) = W(r)V_\ell(r),
\qquad
W(r) = W_L(r)W_R(r).
```

### Left ramp

The left ramp suppresses the potential close to the horizon.

| Region | Value of \(W_L(r)\) |
|---|---|
| ```mathr \le r_{\min}\)``` | \(0\) |
| \(r_{\min} < r < r_{\rm ramp}\) | \(\sin^2\left[\dfrac{\pi}{2}\dfrac{r-r_{\min}}{r_{\rm ramp}-r_{\min}}\right]\) |
| \(r \ge r_{\rm ramp}\) | \(1\) |

### Right cutoff

The right cutoff suppresses the potential near the outer boundary.

| Region | Value of \(W_R(r)\) |
|---|---|
| \(r \le r_{\rm start}\) | \(1\) |
| \(r_{\rm start} < r < r_{\rm cut}\) | \(\cos^2\left[\dfrac{\pi}{2}\dfrac{r-r_{\rm start}}{r_{\rm cut}-r_{\rm start}}\right]\) |
| \(r \ge r_{\rm cut}\) | \(0\) |

Default parameters:

```math
r_{\min}=2.01M,\qquad
r_{\rm ramp}=2.2M,\qquad
r_{\rm start}=500M,\qquad
r_{\rm cut}=1000M.
```

---

## Numerical methods

### Spatial discretisation

The code uses 4th-order finite differences in the tortoise coordinate \(s=r_*\).

#### First derivative

The first derivative is used, for example, to initialise \(\Pi=\partial_t\Psi\) from \(\Psi\):

```math
u'_i =
\frac{
-u_{i+2}
+8u_{i+1}
-8u_{i-1}
+u_{i-2}
}{12\Delta s}
+
\mathcal{O}(\Delta s^4).
```

#### Spatial operator

The spatial operator is

```math
\mathcal{L}[u]_i =
\frac{
-u_{i+2}
+16u_{i+1}
-30u_i
+16u_{i-1}
-u_{i-2}
}{12\Delta s^2}
-
V_i u_i
+
\mathcal{O}(\Delta s^4).
```

One-sided stencils of equivalent order are applied at the four boundary points.

### Time integration

The wave equation is recast as the first-order system

```math
\Pi \equiv \partial_t\Psi,
```

so that

```math
\partial_t
\begin{pmatrix}
\Psi \\
\Pi
\end{pmatrix}
=
\begin{pmatrix}
\Pi \\
\mathcal{L}[\Psi]
\end{pmatrix}.
```

The system is evolved using the standard 4th-order Runge--Kutta method:

```math
\mathbf{u}^{n+1}
=
\mathbf{u}^n
+
\frac{\Delta t}{6}
\left(
\mathbf{k}_1
+2\mathbf{k}_2
+2\mathbf{k}_3
+\mathbf{k}_4
\right).
```

### CFL condition

Numerical stability requires

```math
\Delta t =
\lambda_{\rm CFL}\Delta s,
\qquad
\lambda_{\rm CFL} \leq 1 .
```

The default value is

```math
\lambda_{\rm CFL}=0.1.
```

---

## Initial data

The default initial wavepacket is a modulated Gaussian:

```math
\Psi(s,0)
=
A\exp\left[
-\frac{(s-s_0)^2}{2\sigma^2}
\right]
\cos\left[
\omega_0(s-s_0)
\right].
```

The parameters are:

- \(A\): amplitude,
- \(s_0\): centre of the pulse in tortoise coordinates,
- \(\sigma\): Gaussian width,
- \(\omega_0\): carrier frequency.

---

## Observer extraction

The waveform at an arbitrary observer location \(r_*^{\rm obs}\) is extracted using 4-point Lagrange interpolation over the nearest grid nodes.

The interpolation formula is

```math
\Psi(r_*^{\rm obs},t)
=
\sum_{j=0}^{3}
\Psi(s_{i_j},t)
\prod_{\substack{k=0 \\ k\neq j}}^{3}
\frac{
r_*^{\rm obs}-s_{i_k}
}{
s_{i_j}-s_{i_k}
}.
```

---

## API reference

### `r_of_s(s, M=1.0)`

Converts the tortoise coordinate \(r_*\) to the Schwarzschild radius \(r\).

| Parameter | Type | Description |
|---|---|---|
| `s` | `array_like` | Tortoise coordinates \(r_*\). |
| `M` | `float` | Black hole mass in geometric units. Default: `1.0`. |

Returns:

- `ndarray` containing the corresponding Schwarzschild radii.

---

### `SchwarzPot(r, l, M=1.0)`

Evaluates the Regge--Wheeler potential \(V_\ell(r)\).

| Parameter | Type | Description |
|---|---|---|
| `r` | `array_like` | Schwarzschild radii. |
| `l` | `int` | Angular momentum number \(\ell\). |
| `M` | `float` | Black hole mass. Default: `1.0`. |

Returns:

- `ndarray` containing the potential values.

---

### `W_window(r, M=1.0, rmin=2.01, rramp=2.2, rstart=500.0, rcut=1000.0)`

Computes the smooth window function \(W(r)\).

| Parameter | Type | Description |
|---|---|---|
| `r` | `array_like` | Schwarzschild radii. |
| `M` | `float` | Black hole mass. Default: `1.0`. |
| `rmin` | `float` | Radius where the left ramp starts. |
| `rramp` | `float` | Radius where the left ramp reaches one. |
| `rstart` | `float` | Radius where the right cutoff starts. |
| `rcut` | `float` | Radius where the right cutoff reaches zero. |

Returns:

- `ndarray` containing the window values.

---

### `gaussiana(A, x_min, x_max, N, w0, s0, sigma)`

Constructs the initial Gaussian pulse on a uniform tortoise-coordinate grid.

| Parameter | Type | Description |
|---|---|---|
| `A` | `float` | Amplitude. |
| `x_min` | `float` | Left boundary of the \(r_*\) grid. |
| `x_max` | `float` | Right boundary of the \(r_*\) grid. |
| `N` | `int` | Number of grid points. |
| `w0` | `float` | Carrier frequency \(\omega_0\). |
| `s0` | `float` | Pulse centre in \(r_*\). |
| `sigma` | `float` | Gaussian width \(\sigma\). |

Returns:

- `ndarray` containing the initial field profile.

---

### `simula_full(N, x_min, x_max, Psi_initial, TF=400.0, CFL=0.1, l=2, M=1.0)`

Main time-evolution routine.

| Parameter | Type | Description |
|---|---|---|
| `N` | `int` | Number of grid points. |
| `x_min` | `float` | Left boundary of the \(r_*\) grid. |
| `x_max` | `float` | Right boundary of the \(r_*\) grid. |
| `Psi_initial` | `array_like` | Initial field array. |
| `TF` | `float` | Final evolution time. Default: `400.0`. |
| `CFL` | `float` | Courant factor \(\lambda_{\rm CFL}\). Default: `0.1`. |
| `l` | `int` | Angular momentum number \(\ell\). Default: `2`. |
| `M` | `float` | Black hole mass. Default: `1.0`. |

Returns:

- `U`: field snapshots with shape `(steps, N)`,
- `T`: time array,
- `s`: tortoise-coordinate grid.

---

### `extract_observer_4th(U, s, target)`

Extracts the waveform at a fixed observer position using 4-point Lagrange interpolation.

| Parameter | Type | Description |
|---|---|---|
| `U` | `array_like` | Field snapshots with shape `(steps, N)`. |
| `s` | `array_like` | Tortoise-coordinate grid. |
| `target` | `float` | Observer location \(r_*^{\rm obs}\). |

Returns:

- 1-D array containing \(\Psi(r_*^{\rm obs},t_n)\).

---

## Quick start

```python
import numpy as np
import matplotlib.pyplot as plt

from simulator import gaussiana, simula_full, extract_observer_4th

M = 1.0
N = 2000
x_min = -50.0
x_max = 300.0

# Initial Gaussian pulse centred at r* = 10 M
Psi0 = gaussiana(
    A=1.0,
    x_min=x_min,
    x_max=x_max,
    N=N,
    w0=0.5,
    s0=10.0,
    sigma=3.0,
)

# Evolve up to t = 300 M
U, T, s = simula_full(
    N=N,
    x_min=x_min,
    x_max=x_max,
    Psi_initial=Psi0,
    TF=300.0,
    CFL=0.1,
    l=2,
    M=M,
)

# Extract waveform at r* = 100 M
obs = extract_observer_4th(U, s, target=100.0)

plt.plot(T, obs)
plt.xlabel(r"$t\,[M]$")
plt.ylabel(r"$\Psi$")
plt.title(r"Quasi-normal ringdown -- Schwarzschild BH, $\ell=2$")
plt.tight_layout()
plt.show()
```

---

## Dependencies

| Package | Role | Status |
|---|---|---|
| `numpy` | Array operations and grid generation. | Required |
| `scipy` | Lambert W function via `scipy.special.lambertw`. | Preferred |
| `mpmath` | Fallback implementation of Lambert W. | Optional fallback |
| `matplotlib` | Plotting. | Optional |

Install them with

```bash
pip install numpy scipy mpmath matplotlib
```

---

## Notes

- The coordinate called `s` in the code corresponds to the tortoise coordinate \(r_*\).
- The dominant gravitational perturbation is usually the \(\ell=2\) mode.
- The potential window is introduced only to reduce boundary reflections; it should not affect the physical region of interest.
- All quantities are written in geometric units.
