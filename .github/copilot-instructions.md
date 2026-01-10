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
| `run_experiment.sh` | End-to-end script: build → run → postprocess |
| `visor_cfd` | Compiled executable (after `make`) |

### Postprocessing (Python)
| File | Purpose |
|------|---------|
| `plot_results.py` | Reads binary snapshots, creates animations (GIF), drag-comparison plots, and cross-section visualizations |

### Output
| Directory/File | Contents |
|----------------|----------|
| `output/snapshots/` | Binary flow-field snapshots for animation |
| `output/visor_results.csv` | Drag results per visor angle |
| `output/grid_info.csv` | Grid metadata |
| `output/z_coords.csv` | Vertical cell coordinates |
| `output/drag_comparison.*` | Drag analysis plots (PNG, EPS, SVG) |
| `output/grid_resolution.*` | Grid resolution plots (PNG, EPS, SVG) |
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

# Run the full experiment (uses 8 threads by default)
OMP_NUM_THREADS=8 ./visor_cfd

# Or use the all-in-one workflow script
./run_experiment.sh
```

## Postprocessing

```bash
# Activate Python environment and run plotting
source venv/bin/activate
python plot_results.py
```

Outputs:
- `output/drag_comparison.{png,eps,svg}` — Drag analysis plots
- `output/grid_resolution.{png,eps,svg}` — Grid resolution visualization
- `output/cross_sections_*.{png,eps,svg}` — X-Z, Y-Z, X-Y flow cross-sections with masked runner
- `output/animation_*.gif` — Flow field animations (one per visor configuration)

## Key Solver Details
- **Method**: Finite Volume Method (FVM) with explicit time-stepping
- **Convection**: First-order upwind
- **Pressure**: SIMPLE-like correction with Jacobi iterations
- **Parallelization**: OpenMP (`!$OMP PARALLEL DO` on main loops)
- **Boundary Conditions**: Uniform inlet (8 m/s), zero-gradient outlet, slip walls
