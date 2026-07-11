# fortran-pde

Poisson/Laplace and transport (advection-diffusion) PDE solvers, written
from scratch in modern Fortran. No external PDE libraries — only
LAPACK/BLAS for the one dense linear solve used in 1D.

Verified: the 2D Poisson solver is checked against a manufactured
solution and reproduces the expected **2nd-order** convergence rate
(see `test/check_poisson.f90` — halving the grid spacing quarters the
error, as it should for central differences).

## What's solved

**Poisson / Laplace equation**, for *any* source function `f`:

```
-Laplacian(u) = f(x)      (1D: -u'' = f(x))
           u  = g(x)      on the boundary (Dirichlet)
```

`f` and the boundary function `g` are ordinary Fortran functions you
pass in as arguments — there's no restriction on what `f` looks like
(smooth, localized, oscillatory, a sum of Gaussians, ...). See
`src/problem_definitions_module.f90` for examples, including a
two-source "custom" case, and add your own following the same pattern.

**Transport (advection-diffusion) equation**, in 1D and 2D, with time:

```
du/dt + v . grad(u) = D * Laplacian(u)
```

## Repo structure

```
fortran-pde/
├── fpm.toml                 # manifest for the Fortran Package Manager (fpm)
├── Makefile                 # plain-Makefile fallback (no fpm needed)
├── src/
│   ├── kinds_module.f90              # real64 precision, central single source of truth
│   ├── grid_module.f90               # structured Cartesian grid (1D/2D)
│   ├── io_module.f90                 # write fields to plain text for plotting
│   ├── poisson_module.f90            # -Laplacian(u)=f solver (1D direct, 2D matrix-free CG)
│   ├── transport_module.f90          # du/dt + v.grad(u) = D*Laplacian(u), explicit time-stepping
│   └── problem_definitions_module.f90  # <- EDIT THIS to change f(x), initial/boundary conditions
├── app/
│   ├── main_poisson.f90     # driver: solves 1D + 2D Poisson, writes results/*.dat
│   └── main_transport.f90   # driver: time-steps 1D + 2D transport, writes snapshot series
├── test/
│   └── check_poisson.f90    # manufactured-solution convergence test
├── scripts/
│   ├── plot_results.py      # matplotlib: contour plots, snapshot overlays
│   └── animate_transport.py # Manim scene animating the transport snapshots
└── results/                 # solver output (.dat files), created on first run
```

## Building and running

### Option A — plain `make` (no extra tools needed beyond gfortran + LAPACK/BLAS)

```bash
sudo apt install gfortran liblapack-dev libopenblas-dev   # Ubuntu/Debian
# or: conda install -c conda-forge gfortran lapack openblas

make run-poisson      # builds + runs the Poisson driver
make run-transport    # builds + runs the transport driver
make test             # builds + runs the convergence test
```

### Option B — `fpm` (once you have it installed)

```bash
fpm build
fpm run --example main_poisson     # or: fpm run main_poisson  (name depends on fpm version)
fpm run main_transport
fpm test
```

The `fpm.toml` is already set up with `link = ["lapack", "blas"]`.

## Visualizing results

```bash
pip install numpy matplotlib manim

python scripts/plot_results.py poisson2d results/poisson2d_custom_source.dat
python scripts/plot_results.py transport1d "results/transport1d_t*.dat"

# Animated version (renders an .mp4):
manim -qm -pql scripts/animate_transport.py TransportWave
```

## How the numerics work (short version)

- **Discretization**: standard second-order central finite differences
  (3-point stencil in 1D, 5-point stencil in 2D).
- **1D Poisson**: the resulting tridiagonal system is solved directly
  with LAPACK's `dgtsv` — exact (to machine precision) for a given grid.
- **2D Poisson**: the discrete Laplacian is symmetric positive definite,
  so it's solved with a **matrix-free Conjugate Gradient** method — the
  sparse matrix is never assembled or stored; we only ever apply the
  5-point stencil (`apply_laplacian_2d`). This scales to much larger
  grids than a dense solve would.
- **Transport**: explicit forward-Euler in time, first-order upwind for
  advection (stable but numerically diffusive), central differences for
  the diffusion term. `max_stable_dt_1d`/`max_stable_dt_2d` compute a
  timestep that respects both the CFL condition and the diffusion-number
  stability limit, so you don't have to hand-tune `dt`.

## Extending to 3D / your own equation

- Add a `make_grid_3d` to `grid_module.f90` (same pattern as the 2D one,
  one more coordinate array).
- Add a 7-point-stencil version of `apply_laplacian_2d` in
  `poisson_module.f90` — the CG solver itself (`solve_poisson_2d`'s loop
  logic) doesn't change, only the stencil.
- For the transport equation, add a `vz`/`z` term to
  `step_transport_2d` following the same upwind+central pattern.
- To change *what* you're solving for, you almost never need to touch
  the solver modules — just add a new function to
  `problem_definitions_module.f90` and pass it into the driver.

## Extending the numerics

- Swap first-order upwind for a higher-order scheme (e.g. QUICK, or a
  flux limiter for a proper finite-volume formulation) inside
  `transport_module.f90` — this is the natural next step if you want a
  genuine FVM code rather than FDM.
- Swap the 2D CG solver for a real FEM assembly (unstructured
  triangular mesh, weak form, element stiffness matrices) once you're
  comfortable with the structured-grid version here — the boundary
  condition handling and I/O modules can be reused as-is.
