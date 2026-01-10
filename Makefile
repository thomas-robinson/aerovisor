# Makefile for Visor Aerodynamics CFD Simulation
# Using gfortran with OpenMP support

# Compiler
FC = gfortran

# Compiler flags
# -O3: High optimization
# -fopenmp: OpenMP parallelization
# -march=native: Optimize for current CPU
# -ffast-math: Aggressive floating-point optimizations
FFLAGS = -O3 -fopenmp -march=native -ffast-math

# Debug flags (uncomment for debugging)
# FFLAGS = -g -Wall -Wextra -fcheck=all -fbacktrace -fopenmp

# Source directory
SRCDIR = src

# Object files (in dependency order)
OBJS = grid_module.o geometry_module.o solver_module.o visor_experiment.o

# Executable name
TARGET = visor_cfd

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo ""
	@echo "Build complete: $(TARGET)"
	@echo "Run with: ./$(TARGET)"
	@echo ""
	@echo "For parallel execution, set OMP_NUM_THREADS:"
	@echo "  export OMP_NUM_THREADS=8"
	@echo "  ./$(TARGET)"

# Compile rules (order matters for module dependencies)
grid_module.o: $(SRCDIR)/grid_module.f90
	$(FC) $(FFLAGS) -c $<

geometry_module.o: $(SRCDIR)/geometry_module.f90 grid_module.o
	$(FC) $(FFLAGS) -c $<

solver_module.o: $(SRCDIR)/solver_module.f90 grid_module.o
	$(FC) $(FFLAGS) -c $<

visor_experiment.o: $(SRCDIR)/visor_experiment.f90 grid_module.o geometry_module.o solver_module.o
	$(FC) $(FFLAGS) -c $<

# Clean
clean:
	rm -f *.o *.mod $(TARGET)

# Clean all including results
cleanall: clean
	rm -f visor_results.csv

# Quick test (fewer iterations)
test: $(TARGET)
	@echo "Running quick test..."
	./$(TARGET)

# Help
help:
	@echo "Visor Aerodynamics CFD Simulation"
	@echo ""
	@echo "Targets:"
	@echo "  all      - Build the simulation (default)"
	@echo "  clean    - Remove object files and executable"
	@echo "  cleanall - Remove all generated files including results"
	@echo "  test     - Build and run simulation"
	@echo "  help     - Show this help"
	@echo ""
	@echo "Configuration:"
	@echo "  Grid: 60 x 45 x 150 cells"
	@echo "  Domain: 2.0m x 1.7m x 2.5m (reduced height)"
	@echo "  Non-uniform vertical grid: 3mm near head, coarser elsewhere"
	@echo "  Visor thickness: 5mm (realistic)"
	@echo ""
	@echo "Environment variables:"
	@echo "  OMP_NUM_THREADS - Number of OpenMP threads (default: all cores)"

.PHONY: all clean cleanall test help
