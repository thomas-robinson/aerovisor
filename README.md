# Aerovisor

A 3D computational fluid dynamics (CFD) simulation for analyzing the aerodynamic effects of a visor on a runner during a triathlon.

## Overview

This project simulates airflow around a runner at race pace (8 m/s ≈ 5:00/km) to study how different visor angles affect aerodynamic drag. The solver is written in modern Fortran with OpenMP parallelization, and includes Python scripts for postprocessing and visualization.

## Features

- **3D Incompressible Flow Solver** — Finite Volume Method with SIMPLE pressure correction
- **Parametric Visor Study** — Automatic sweep of visor angles from -45° to +45°
- **Non-Uniform Grid** — Refined resolution near the head/visor region for accuracy
- **OpenMP Parallelization** — Efficient multi-core execution
- **Automated Workflow** — Single script to build, run, and generate all outputs

## Quick Start

```bash
# Run the complete workflow (build + simulate + postprocess)
./run_experiment.sh
```

This will:
1. Build the Fortran solver
2. Set up a Python virtual environment with dependencies
3. Run the CFD simulation for baseline and 13 visor angles
4. Generate drag analysis plots, cross-section visualizations, and flow animations

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
│   ├── visor_results.csv     # Drag data for all configurations
│   ├── drag_comparison.*     # Drag analysis plots (PNG/EPS/SVG)
│   ├── grid_resolution.*     # Grid visualization (PNG/EPS/SVG)
│   ├── cross_sections_*.*    # Flow cross-sections (PNG/EPS/SVG)
│   └── animation_*.gif       # Flow field animations
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
The experiment tests 14 configurations:
- 1 baseline (no visor)
- 13 visor angles: -45°, -40°, ..., 0°, ..., +40°, +45° (by 7.5° increments)

## Manual Build & Run

```bash
# Build
make clean && make

# Run with custom thread count
OMP_NUM_THREADS=4 ./visor_cfd

# Postprocess
source venv/bin/activate
python plot_results.py
```

## Output Files

### Data
- `visor_results.csv` — Drag force for each configuration
- `grid_info.csv` — Grid dimensions and domain size
- `z_coords.csv` — Vertical cell coordinates

### Visualizations
- `drag_comparison.*` — Multi-panel drag analysis (absolute, relative, polar)
- `grid_resolution.*` — Non-uniform grid spacing visualization
- `cross_sections_*.*` — X-Z (side), Y-Z (front), X-Y (top) flow views
- `animation_*.gif` — Time evolution of flow field

## License

MIT License

## Author

Tom
