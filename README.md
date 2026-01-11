# Aerovisor

A 3D computational fluid dynamics (CFD) simulation for analyzing the aerodynamic effects of a visor on a runner during a triathlon.

## Overview

This project simulates airflow around a runner at race pace to study how different visor angles affect aerodynamic drag. The solver is written in modern Fortran with OpenMP parallelization, and includes Python scripts for postprocessing and visualization.

## Features

- **3D Incompressible Flow Solver** — Finite Volume Method with SIMPLE pressure correction
- **Parametric Visor Study** — Automatic sweep of visor angles from -45° to +45°
- **Non-Uniform Grid** — Refined resolution near the head/visor region for accuracy
- **OpenMP Parallelization** — Efficient multi-core execution
- **Configurable Runtime** — Namelist input file or command-line parameters
- **Automated Workflow** — Single script to build, run, and generate all outputs

## Quick Start

```bash
# Quick test run (500 iterations, ~5 min)
./run_experiment.sh 500

# Production run (2000 iterations, ~20 min, well converged)
./run_experiment.sh 2000

# Use default from visor.nml
./run_experiment.sh
```

This will:
1. Build the Fortran solver
2. Set up a Python virtual environment with dependencies
3. Run the CFD simulation for baseline and 13 visor angles
4. Generate drag analysis plots, cross-section visualizations, and flow animations

## Runtime Configuration

Edit `visor.nml` to change default parameters:

```fortran
&experiment
    n_iterations = 500      ! 500=quick, 2000=production, 5000=high-accuracy
    report_interval = 50    ! Progress reporting frequency
    snapshot_interval = 50  ! Animation frame frequency
    inlet_velocity = 5.0    ! m/s (5.0 = 18 km/h running pace)
/
```

Or override iterations from command line:
```bash
./visor_cfd 2000           # Run solver directly with 2000 iterations
./run_experiment.sh 2000   # Full workflow with 2000 iterations
```

Output files are automatically tagged with iteration count (e.g., `visor_results_2000.csv`).

## Requirements

### Fortran Compiler
- **gfortran** ≥ 9 with OpenMP support

### Python (for postprocessing)
- Python 3.8+
- NumPy
- Matplotlib
- Pandas

Dependencies are automatically installed by `run_experiment.sh`.

## Project Structure

```
aerovisor/
├── src/                      # Fortran source files
│   ├── grid_module.f90       # Non-uniform grid generation
│   ├── geometry_module.f90   # Runner & visor geometry
│   ├── solver_module.f90     # FVM solver with pressure correction
│   └── visor_experiment.f90  # Experiment driver
├── output/                   # Generated outputs (after running)
│   ├── snapshots/            # Binary flow-field snapshots
│   ├── visor_results_*.csv   # Drag data for all configurations
│   ├── drag_comparison_*.*   # Drag analysis plots (PNG/EPS/SVG)
│   ├── grid_resolution_*.*   # Grid visualization (PNG/EPS/SVG)
│   ├── cross_sections_*.*    # Flow cross-sections (PNG/EPS/SVG)
│   └── animation_*.gif       # Flow field animations
├── visor.nml                 # Runtime configuration (namelist)
├── Makefile                  # Build configuration
├── run_experiment.sh         # Automated workflow script
├── plot_results.py           # Postprocessing & visualization
└── requirements.txt          # Python dependencies
```

## Simulation Details

### Domain Configuration
| Parameter | Value |
|-----------|-------|
| Domain Size | 2.0 m × 1.7 m × 2.5 m (L × W × H) |
| Grid Cells | 60 × 45 × 150 |
| Inlet Velocity | 8 m/s (≈ 5:00/km pace) |
| Visor Thickness | 5 mm |
| Grid Refinement | 4.2 mm vertical resolution near head |

### Solver Parameters
- **Method**: Finite Volume Method (FVM)
- **Time Stepping**: Explicit with CFL = 0.2
- **Pressure**: SIMPLE-like correction with Jacobi iterations
- **Convection**: First-order upwind
- **Boundary Conditions**: Uniform inlet, zero-gradient outlet, slip walls

### Visor Configurations
The experiment tests 20 configurations:
- 1 baseline (no visor)
- 19 visor angles: -45°, -40°, -35°, ..., 0°, ..., +35°, +40°, +45° (5° increments)

## Manual Build & Run

```bash
# Build
make clean && make

# Run with custom thread count and iterations
OMP_NUM_THREADS=4 ./visor_cfd 2000

# Postprocess specific iteration run
source venv/bin/activate
python plot_results.py 2000
```

## Output Files

Output files include iteration count in filename (e.g., `_500` or `_2000`):

### Data
- `visor_results_*.csv` — Drag force for each configuration
- `grid_info_*.csv` — Grid dimensions and domain size
- `z_coords_*.csv` — Vertical cell coordinates

### Visualizations
- `drag_comparison_*.*` — Multi-panel drag analysis (absolute, relative, polar)
- `grid_resolution_*.*` — Non-uniform grid spacing visualization
- `cross_sections_*.*` — X-Z (side), Y-Z (front), X-Y (top) flow views
- `animation_*.gif` — Time evolution of flow field

## Iteration Guidelines

| Iterations | Time | Use Case |
|------------|------|----------|
| 100-500 | ~2-5 min | Quick tests, debugging |
| 500 | ~5 min | Fast comparison (may not be fully converged) |
| 2000 | ~20 min | Production runs (well converged) |
| 5000 | ~50 min | High-accuracy validation |

## License

MIT License

## Author

@thomas-robinson Tom Robinson 
