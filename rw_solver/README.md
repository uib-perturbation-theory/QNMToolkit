# 4th order Regge-Wheeler time-domain solver

Numerical solver for the **Regge–Wheeler equation** describing gravitational perturbations (quasi-normal modes) of a Schwarzschild black hole. The code evolves a scalar/gravitational wave in tortoise coordinates using a 4th-order Runge–Kutta scheme in time and 4th-order finite differences in space.

---

## Physics Background

The **Regge–Wheeler equation** governs odd-parity perturbations of a Schwarzschild black hole of mass $M$:

$$
\left( -\frac{\partial^2}{\partial t^2} + \frac{\partial^2}{\partial r_*^2} - V_\ell(r) \right) \Psi = 0
$$

where $r_*$ is the tortoise coordinate and $V_\ell(r)$ is the Regge–Wheeler potential. The Schwarzschild metric in standard coordinates is:

$$
ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2 d\Omega^2
$$

The event horizon sits at $r = 2M$.

---

## Coordinate System

### Tortoise Coordinate

The tortoise coordinate $r_*$ (labelled `s` in the code) is defined by:

$$
r_* = r + 2M \ln\left(\frac{r}{2M} - 1\right)
$$

which maps the exterior region $r \in (2M, \infty)$ to $r_* \in (-\infty, +\infty)$.

### Inverting the Tortoise Coordinate

Given $r_{*}$, the physical radius $r$ is recovered via the **Lambert W function** $W_{0}$:

$$
r(r_*) = 2M\left[1 + W_0\left(\exp\left(\frac{r_{*}}{2M} - 1\right)\right)\right]
$$

Implemented in r_of_s(s, M), using scipy.special.lambertw when available, or mpmath.lambertw as fallback.

## Potential

### Regge–Wheeler Potential

The effective potential for gravitational perturbations of angular momentum $\ell$ is:

$$
V_\ell(r) = \left(1 - \frac{2M}{r}\right)\left(\frac{\ell(\ell+1)}{r^2} - \frac{6M}{r^3}\right)
$$

- The prefactor $(1 - 2M/r)$ is the **redshift factor**, vanishing at the horizon.
- The $\ell(\ell+1)/r^2$ term is the **centrifugal barrier**.
- The $-6M/r^3$ term is the **gravitational curvature** correction (absent for scalar fields).

For physical gravitational waves the dominant mode is $\ell = 2$.

---

## Windowing

To suppress spurious reflections at the grid boundaries, the potential is multiplied by a smooth window function $W(r)$:

$$
V_{\rm eff}(r) = W(r)\cdot V_\ell(r), \qquad W = W_L \cdot W_R
$$

**Left ramp** (near-horizon suppression):

$$
W_L(r) = \begin{cases}
0 & r \leq r_{\min} \\[4pt]
\sin^2\!\left(\dfrac{\pi}{2}\,\dfrac{r - r_{\min}}{r_{\rm ramp} - r_{\min}}\right) & r_{\min} < r < r_{\rm ramp} \\[4pt]
1 & r \geq r_{\rm ramp}
\end{cases}
$$

**Right cutoff** (outer-boundary suppression):

$$
W_R(r) = \begin{cases}
1 & r \leq r_{\rm start} \\[4pt]
\cos^2\!\left(\dfrac{\pi}{2}\,\dfrac{r - r_{\rm start}}{r_{\rm cut} - r_{\rm start}}\right) & r_{\rm start} < r < r_{\rm cut} \\[4pt]
0 & r \geq r_{\rm cut}
\end{cases}
$$

Default parameters: $r_{\min} = 2.01M$, $r_{\rm ramp} = 2.2M$, $r_{\rm start} = 500M$, $r_{\rm cut} = 1000M$.

---

## Numerical Methods

### Spatial Discretisation — 4th-Order Finite Differences

**First derivative** (used to initialise $\Pi = \partial_t\Psi$ from $\Psi$):

$$
u'_i = \frac{-u_{i+2} + 8u_{i+1} - 8u_{i-1} + u_{i-2}}{12\,\Delta s} + \mathcal{O}(\Delta s^4)
$$

**Spatial operator** $\mathcal{L}$ (second derivative minus potential):

$$
\mathcal{L}[u]_i = \frac{-u_{i+2} + 16u_{i+1} - 30u_i + 16u_{i-1} - u_{i-2}}{12\,\Delta s^2} - V_i\,u_i + \mathcal{O}(\Delta s^4)
$$

One-sided stencils of equivalent order are applied at the four boundary points.

### Time Integration — 4th-Order Runge–Kutta

The wave equation is recast as the first-order system with $\Pi \equiv \partial_t \Psi$:

$$
\partial_t \begin{pmatrix} \Psi \\ \Pi \end{pmatrix} = \begin{pmatrix} \Pi \\ \mathcal{L}[\Psi] \end{pmatrix}
$$

Standard RK4 advances the state by one timestep $\Delta t$:

$$
\mathbf{u}^{n+1} = \mathbf{u}^n + \frac{\Delta t}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)
$$

### CFL Condition

Numerical stability requires:

$$
\Delta t = \lambda_{\rm CFL}\,\Delta s, \qquad \lambda_{\rm CFL} \leq 1
$$

Default: $\lambda_{\rm CFL} = 0.1$.

---

## Initial Data — Gaussian Pulse

The default initial wavepacket is a **modulated Gaussian**:

$$
\Psi(s,\,0) = A\,\exp\!\left(-\frac{(s-s_0)^2}{2\sigma^2}\right)\cos\!\left(\omega_0(s-s_0)\right)
$$

where $A$ is the amplitude, $s_0$ the centre, $\sigma$ the width, and $\omega_0$ the carrier frequency.

---

## Observer Extraction — 4th-Order Lagrange Interpolation

The waveform at an arbitrary observer location $r_*^{\rm obs}$ is extracted using 4-point Lagrange interpolation over the nearest grid nodes $\{s_{i_j}\}_{j=0}^{3}$:

$$
\Psi\!\left(r_*^{\rm obs},\,t\right) = \sum_{j=0}^{3} \Psi(s_{i_j},\,t)\prod_{\substack{k=0\\k\neq j}}^{3} \frac{r_*^{\rm obs} - s_{i_k}}{s_{i_j} - s_{i_k}}
$$

---

## API Reference

### `r_of_s(s, M=1.0)`
Converts tortoise coordinate $r_*$ to Schwarzschild radius $r$.

| Parameter | Type | Description |
|-----------|------|-------------|
| `s` | `array_like` | Tortoise coordinates $r_*$ |
| `M` | `float` | Black hole mass (geometric units, default `1.0`) |

**Returns:** `ndarray` of radii $r$.

---

### `SchwarzPot(r, l, M=1.0)`
Evaluates the Regge–Wheeler potential $V_\ell(r)$.

| Parameter | Type | Description |
|-----------|------|-------------|
| `r` | `array_like` | Schwarzschild radii |
| `l` | `int` | Angular momentum number $\ell$ |
| `M` | `float` | Black hole mass |

---

### `W_window(r, M=1.0, rmin=2.01, rramp=2.2, rstart=500.0, rcut=1000.0)`
Computes the smooth window function $W(r)$.

---

### `gaussiana(A, x_min, x_max, N, w0, s0, sigma)`
Constructs the initial Gaussian pulse on a uniform $r_*$ grid.

| Parameter | Description |
|-----------|-------------|
| `A` | Amplitude |
| `w0` | Carrier frequency $\omega_0$ |
| `s0` | Pulse centre in $r_*$ |
| `sigma` | Pulse width $\sigma$ |

---

### `simula_full(N, x_min, x_max, Psi_initial, TF=400.0, CFL=0.1, l=2, M=1.0)`
Main time-evolution routine.

| Parameter | Description |
|-----------|-------------|
| `N` | Number of grid points |
| `x_min`, `x_max` | Grid extent in $r_*$ |
| `Psi_initial` | Initial field array |
| `TF` | Final evolution time |
| `CFL` | Courant factor $\lambda_{\rm CFL}$ |
| `l` | Angular momentum $\ell$ |
| `M` | Black hole mass |

**Returns:** `(U, T, s)` — field snapshots `(steps, N)`, time array, tortoise grid.

---

### `extract_observer_4th(U, s, target)`
Extracts the waveform at a fixed observer position via Lagrange interpolation.

| Parameter | Description |
|-----------|-------------|
| `U` | Field snapshots array `(steps, N)` |
| `s` | Tortoise coordinate grid |
| `target` | Observer location $r_*^{\rm obs}$ |

**Returns:** 1-D array of $\Psi(r_*^{\rm obs}, t_n)$.

---

## Quick Start

```python
import numpy as np
import matplotlib.pyplot as plt
from simulator import gaussiana, simula_full, extract_observer_4th

M     = 1.0
N     = 2000
x_min = -50.0
x_max = 300.0

# Initial Gaussian pulse centred at r* = 10 M
Psi0 = gaussiana(A=1.0, x_min=x_min, x_max=x_max,
                 N=N, w0=0.5, s0=10.0, sigma=3.0)

# Evolve up to t = 300 M
U, T, s = simula_full(N, x_min, x_max, Psi0,
                      TF=300.0, CFL=0.1, l=2, M=M)

# Extract waveform at r* = 100 M
obs = extract_observer_4th(U, s, target=100.0)

plt.plot(T, obs)
plt.xlabel(r"$t\;[M]$")
plt.ylabel(r"$\Psi$")
plt.title("Quasi-normal ringdown — Schwarzschild BH ($\ell=2$)")
plt.tight_layout()
plt.show()
```

---

## Dependencies

| Package | Role | Status |
|---------|------|--------|
| `numpy` | Array operations, grid generation | required |
| `scipy` | Lambert W via `scipy.special.lambertw` | preferred |
| `mpmath` | Lambert W fallback | auto-fallback |
| `matplotlib` | Plotting (not imported in core) | optional |

```bash
pip install numpy scipy mpmath matplotlib
```
