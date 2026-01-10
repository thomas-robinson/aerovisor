#!/bin/bash
#===============================================================================
# run_experiment.sh
# Complete workflow for visor aerodynamics experiment
# Builds code, runs simulation, generates animations and plots
#===============================================================================

set -e  # Exit on error

echo "========================================================================"
echo "VISOR AERODYNAMICS EXPERIMENT - COMPLETE WORKFLOW"
echo "========================================================================"
echo ""

# Configuration
export OMP_NUM_THREADS=8  # Use all 8 logical threads
echo "Using $OMP_NUM_THREADS OpenMP threads"
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
time ./visor_cfd
echo ""

# Step 5: Generate plots and animations
echo "Step 5: Generating plots and animations..."
echo "------------------------------------------------------------------------"
python plot_results.py
echo ""

# Step 6: Summary
echo "========================================================================"
echo "WORKFLOW COMPLETE"
echo "========================================================================"
echo ""
echo "Output files:"
echo "  - output/visor_results.csv: Raw drag data"
echo "  - output/grid_info.csv: Grid configuration"
echo "  - output/z_coords.csv: Vertical level coordinates"
echo "  - output/drag_comparison.png: Drag analysis plot"
echo "  - output/grid_resolution.png: Grid resolution plot"
echo "  - output/animation_*.gif: Flow field animations"
echo ""
echo "To view results:"
echo "  open output/drag_comparison.png"
echo "  open output/animation_baseline.gif"
