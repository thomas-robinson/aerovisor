# Aerovisor Project

A 3D computational fluid dynamics (CFD) simulation for analyzing the aerodynamic effects of a visor on a runner during a triathlon. The solver is written in modern Fortran with OpenMP parallelization.

## Domain Configuration
- **Grid**: 60 × 45 × 150 cells (non-uniform in z for visor resolution)
- **Domain Size**: 2.0 m × 1.7 m × 2.5 m (length × width × height)
- **Grid Spacing**: Uniform dx/dy ≈ 3.3 cm; refined dz ≈ 4.2 mm near head (1.6–2.1 m)
- **Runner Position**: (0.9 m, 0.875 m, 0.0 m) — 45% from inlet, centered in width
- **Visor Thickness**: 5 mm (spans ~1.2 vertical cells in refined zone)

## Project Structure

### Fortran Source (`src/`)
| File | Purpose |
|------|---------|
| `grid_module.f90` | Non-uniform grid generation (refined z near head) |
| `geometry_module.f90` | Runner & visor geometry; solid-mask creation |
| `solver_module.f90` | FVM solver: explicit time-stepping, SIMPLE pressure correction, drag computation |
| `visor_experiment.f90` | Experiment driver: baseline + angle sweep, snapshot output, CSV results |

### Build & Run
| File | Purpose |
|------|---------|
| `Makefile` | Compiles Fortran modules and links `visor_cfd` executable |
| `run_experiment.sh` | End-to-end script: build → run → postprocess (accepts iteration count arg) |
| `visor.nml` | Namelist configuration file (iterations, velocity, intervals) |
| `visor_cfd` | Compiled executable (after `make`) |

### Postprocessing (Python)
| File | Purpose |
|------|---------|
| `plot_results.py` | Reads binary snapshots, creates animations (GIF), drag-comparison plots, and cross-section visualizations |

### Output
Output files include iteration count suffix (e.g., `_500`, `_2000`):

| Directory/File | Contents |
|----------------|----------|
| `output/snapshots/` | Binary flow-field snapshots for animation |
| `output/visor_results_*.csv` | Drag results per visor angle |
| `output/grid_info_*.csv` | Grid metadata |
| `output/z_coords_*.csv` | Vertical cell coordinates |
| `output/drag_comparison_*.*` | Drag analysis plots (PNG, EPS, SVG) |
| `output/grid_resolution_*.*` | Grid resolution plots (PNG, EPS, SVG) |
| `output/cross_sections_*.*` | Flow field cross-sections (PNG, EPS, SVG) |
| `output/animation_*.gif` | Flow field animations |

## Dependencies

### Build (Fortran)
- **gfortran** ≥ 9 (with OpenMP support)

### Postprocessing (Python)
- NumPy
- Matplotlib
- Pandas

## Building & Running

```bash
# Build the Fortran solver
make clean && make

# Quick test (500 iterations)
./run_experiment.sh 500

# Production run (2000 iterations, well converged)
./run_experiment.sh 2000

# Or run solver directly with custom iterations
OMP_NUM_THREADS=8 ./visor_cfd 2000
```

## Runtime Configuration

Edit `visor.nml` to change defaults:
```fortran
&experiment
    n_iterations = 500      ! 500=quick, 2000=production
    report_interval = 50
    snapshot_interval = 50
    inlet_velocity = 5.0    ! m/s
/
```

## Postprocessing

```bash
# Activate Python environment and run plotting
source venv/bin/activate
python plot_results.py 2000  # Process 2000-iteration run
```

Outputs (with iteration suffix, e.g., `_2000`):
- `output/drag_comparison_*.{png,eps,svg}` — Drag analysis plots
- `output/grid_resolution_*.{png,eps,svg}` — Grid resolution visualization
- `output/cross_sections_*.{png,eps,svg}` — X-Z, Y-Z, X-Y flow cross-sections with masked runner
- `output/animation_*.gif` — Flow field animations (one per visor configuration)

## Key Solver Details
- **Method**: Finite Volume Method (FVM) with explicit time-stepping
- **Convection**: First-order upwind
- **Pressure**: SIMPLE-like correction with Jacobi iterations
- **Parallelization**: OpenMP (`!$OMP PARALLEL DO` on main loops)
- **Boundary Conditions**: Uniform inlet (5.0 m/s by default), zero-gradient outlet, slip walls
- **Visor Study**: 20 configurations (baseline + 19 angles from -45° to +45° in 5° increments)
