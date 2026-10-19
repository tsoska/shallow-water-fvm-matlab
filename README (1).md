# Finite Volume Simulation of the 1D Shallow Water Equations

A MATLAB implementation of a **one-dimensional finite volume solver for the shallow water equations**, developed as part of a Computational Hydrodynamics course project at the National Technical University of Athens (NTUA).

The project investigates the influence of **approximate Riemann solvers**, **explicit time-integration schemes**, and the **Courant–Friedrichs–Lewy (CFL) number** on a dam-break (Riemann) problem. Numerical solutions are also compared with an analytical reference solution.

## Governing equations

The inviscid, one-dimensional shallow water equations are solved in conservative form:

$$
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}(\mathbf{U})}{\partial x} = 0,
$$

with

$$
\mathbf{U} = \begin{bmatrix} h \\ hu \end{bmatrix}, \qquad
\mathbf{F}(\mathbf{U}) = \begin{bmatrix} hu \\ hu^2 + \tfrac{1}{2}gh^2 \end{bmatrix},
$$

where $h$ is the water depth, $u$ is the depth-averaged velocity, and $g$ is gravitational acceleration. The scripts use `H` for $h$ and `Q = H*u` for the discharge per unit width.

## Numerical methods

The spatial domain is divided into uniform finite volumes. Cell-averaged states are updated using the interface numerical fluxes:

$$
\frac{d\mathbf{U}_i}{dt} = -\frac{\widehat{\mathbf{F}}_{i+1/2} - \widehat{\mathbf{F}}_{i-1/2}}{\Delta x}.
$$

**Approximate Riemann solvers implemented:**

- **Rusanov** (local Lax–Friedrichs)
- **Roe**
- **HLL** (Harten–Lax–van Leer)
- **HLL with Roe-averaged wave speeds** (`hllroe`)

**Explicit time integrators:**

- Forward Euler (FE)
- Classical fourth-order Runge–Kutta (RK4)

The time step is recalculated from the maximum estimated characteristic speed:

$$
\Delta t = C_{\mathrm{max}}\frac{\Delta x}{\max_i\left(|u_i|+\sqrt{gh_i}\right)}.
$$

The implementation uses piecewise-constant cell states at interfaces. The use of RK4 does **not** make the full space–time discretization fourth-order accurate.

## Test case

A one-dimensional dam-break problem is initialized with discontinuous water depth and zero velocity:

$$
h(x,0) = \begin{cases}
50\ \mathrm{m}, & x < 0,\\
1\ \mathrm{m}, & x > 0,
\end{cases}
\qquad u(x,0) = 0.
$$

| Parameter | Value |
| --- | --- |
| Spatial domain | $x \in [-1,1]$ m |
| Number of cells | 1,000 |
| Gravitational acceleration | $9.81$ m/s² |
| Final simulation time | $0.02$ s |
| Boundary treatment | Zero discharge at both domain boundaries, with hydrostatic momentum-flux terms |
| CFL values in stability study | 0.9, 1.2, 2.0 |
| CFL value in solver comparison | 0.9 |

## Repository files

```text
shallow-water-fvm-matlab/
├── README.md
├── stability_analysis.m
└── solver_validation.m
```

### `stability_analysis.m`

Examines the sensitivity of the numerical solution to the time integrator and CFL number, using a selected numerical flux (HLL by default).

- Runs both Forward Euler and RK4.
- Tests $C_{\mathrm{max}} = 0.9$, $1.2$, and $2.0$.
- Plots water-depth and velocity profiles at $t = 0$, $0.0067$, $0.0133$, and $0.02$ s.
- Includes a basic check for non-finite, complex, or non-positive water depths.

### `solver_validation.m`

Compares all four Riemann solvers for both time integrators at $C_{\mathrm{max}} = 0.9$.

- Computes a reference solution for the dam-break Riemann problem using the `swe_analytical` function supplied in the script.
- Overlays numerical and reference profiles for water depth and velocity at $t = 0.02$ s.
- Reports **root-mean-square error (RMSE)** for each variable and solver:

$$
\operatorname{RMSE}(q) = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left(q_i^{\mathrm{num}} - q_i^{\mathrm{ref}}\right)^2}.
$$

The console labels `H_L2` and `u_L2` refer to this discrete RMSE calculation (not a grid-weighted integral $L^2$ norm).

## How to run

1. Download or clone this repository.
2. Open MATLAB and set the **Current Folder** to the directory containing both `.m` files.
3. Run either script in the Command Window:

   ```matlab
   stability_analysis
   ```

   or

   ```matlab
   solver_validation
   ```

No external datasets are required. The scripts use standard MATLAB functions; a MATLAB version supporting **local functions in scripts** and **string arrays** is needed (R2016b or later).

### Changing simulation settings

- In `stability_analysis.m`, set `flux_type` to `'rusanov'`, `'roe'`, `'hll'`, or `'hllroe'`.
- Modify `Cmax_list` to test other Courant numbers.
- Change `nx` to alter spatial resolution.
- In `solver_validation.m`, modify `solvers`, `time_schemes`, or `Cmax` to select comparisons.

## Discussion and limitations

The experiments are intended to illustrate how flux selection and time integration affect numerical dissipation, accuracy, and stability near a discontinuity. As discussed in the accompanying coursework report, the Roe-family implementations can display nonphysical features near the rarefaction region in this test case.

**Scope notes:**

- CFL values above unity are used deliberately to explore stability limits; successful behavior for an individual run is not a general stability guarantee.
- The implementations do not include a dedicated Roe entropy fix or a general positivity-preserving treatment.
- This project compares against a reference solution at one spatial resolution; it does **not** include a formal grid-convergence study.

## Project context

Computational Hydrodynamics coursework, **School of Naval Architecture and Marine Engineering, National Technical University of Athens (NTUA)**.

**Author:** Georgios Katsos
