#!/bin/bash
#===============================================================================
# run_experiment.sh
# Complete workflow for visor aerodynamics experiment
# Builds code, runs simulation, generates animations and plots
#
# Usage:
#   ./run_experiment.sh           # Use default from visor.nml (500 iterations)
#   ./run_experiment.sh 500       # Quick test run
#   ./run_experiment.sh 2000      # Production run (well converged)
#   ./run_experiment.sh 5000      # High-accuracy run
#===============================================================================

set -e  # Exit on error

# Parse command-line argument for iteration count
N_ITERATIONS=${1:-}

echo "========================================================================"
echo "VISOR AERODYNAMICS EXPERIMENT - COMPLETE WORKFLOW"
echo "========================================================================"
echo ""

# Configuration
export OMP_NUM_THREADS=8  # Use all 8 logical threads
echo "Using $OMP_NUM_THREADS OpenMP threads"
if [ -n "$N_ITERATIONS" ]; then
    echo "Iterations per configuration: $N_ITERATIONS (from command line)"
else
    echo "Iterations per configuration: from visor.nml or default"
fi
echo ""

# Step 1: Clean and build
echo "Step 1: Building Fortran code..."
echo "------------------------------------------------------------------------"
make clean
make
echo ""

# Step 2: Create output directory
echo "Step 2: Creating output directories..."
echo "------------------------------------------------------------------------"
mkdir -p output/snapshots
echo "  Created: output/"
echo "  Created: output/snapshots/"
echo ""

# Step 3: Set up Python environment
echo "Step 3: Setting up Python environment..."
echo "------------------------------------------------------------------------"
if [ ! -d "venv" ]; then
    echo "  Creating virtual environment..."
    python3 -m venv venv
fi
source venv/bin/activate
echo "  Installing/updating dependencies from requirements.txt..."
pip install -q -r requirements.txt
echo "  Python environment ready."
echo ""

# Step 4: Run the simulation
echo "Step 4: Running CFD simulation..."
echo "------------------------------------------------------------------------"
echo "  This may take several minutes..."
echo ""
if [ -n "$N_ITERATIONS" ]; then
    time ./visor_cfd "$N_ITERATIONS"
else
    time ./visor_cfd
fi
echo ""

# Step 5: Generate plots and animations
echo "Step 5: Generating plots and animations..."
echo "------------------------------------------------------------------------"
if [ -n "$N_ITERATIONS" ]; then
    python plot_results.py "$N_ITERATIONS"
else
    python plot_results.py
fi
echo ""

# Step 6: Summary
echo "========================================================================"
echo "WORKFLOW COMPLETE"
echo "========================================================================"
echo ""
if [ -n "$N_ITERATIONS" ]; then
    SUFFIX="_${N_ITERATIONS}"
else
    SUFFIX=""
fi
echo "Output files:"
echo "  - output/visor_results${SUFFIX}.csv: Raw drag data"
echo "  - output/grid_info${SUFFIX}.csv: Grid configuration"
echo "  - output/z_coords${SUFFIX}.csv: Vertical level coordinates"
echo "  - output/drag_comparison${SUFFIX}.png: Drag analysis plot"
echo "  - output/grid_resolution${SUFFIX}.png: Grid resolution plot"
echo "  - output/cross_sections_*${SUFFIX}.png: Flow cross-sections"
echo "  - output/animation_*${SUFFIX}.gif: Flow field animations"
echo ""
echo "To view results:"
echo "  open output/drag_comparison${SUFFIX}.png"
echo "  open output/animation_baseline${SUFFIX}.gif"
