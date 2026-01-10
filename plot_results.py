#!/usr/bin/env python3
"""
plot_results.py
Generate animations and comparison plots from visor aerodynamics experiment
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
import pandas as pd
import glob
import struct


def read_snapshot(filename):
    """Read binary snapshot file from Fortran simulation."""
    with open(filename, 'rb') as f:
        # Read dimensions (3 integers)
        nx = struct.unpack('i', f.read(4))[0]
        ny = struct.unpack('i', f.read(4))[0]
        nz = struct.unpack('i', f.read(4))[0]
        
        n_cells = nx * ny * nz
        
        # Read velocity fields (double precision)
        u = np.frombuffer(f.read(n_cells * 8), dtype=np.float64).reshape((nx, ny, nz), order='F')
        v = np.frombuffer(f.read(n_cells * 8), dtype=np.float64).reshape((nx, ny, nz), order='F')
        w = np.frombuffer(f.read(n_cells * 8), dtype=np.float64).reshape((nx, ny, nz), order='F')
        
        # Read pressure
        p = np.frombuffer(f.read(n_cells * 8), dtype=np.float64).reshape((nx, ny, nz), order='F')
        
        # Read solid mask (logical, 1 byte per element in gfortran)
        solid_bytes = f.read(n_cells)
        solid = np.frombuffer(solid_bytes, dtype=np.uint8).reshape((nx, ny, nz), order='F')
        solid = solid.astype(bool)
        
    return {
        'nx': nx, 'ny': ny, 'nz': nz,
        'u': u, 'v': v, 'w': w,
        'p': p, 'solid': solid
    }


def read_grid_info():
    """Read grid information from CSV files."""
    grid_info = pd.read_csv('output/grid_info.csv')
    z_coords = pd.read_csv('output/z_coords.csv')
    
    nx = int(grid_info['nx'].iloc[0])
    ny = int(grid_info['ny'].iloc[0])
    nz = int(grid_info['nz'].iloc[0])
    length = float(grid_info['length'].iloc[0])
    width = float(grid_info['width'].iloc[0])
    height = float(grid_info['height'].iloc[0])
    
    z = z_coords['z'].values
    dz = z_coords['dz'].values
    
    # Create x and y coordinates (uniform)
    x = np.linspace(length / (2*nx), length - length / (2*nx), nx)
    y = np.linspace(width / (2*ny), width - width / (2*ny), ny)
    
    return {
        'nx': nx, 'ny': ny, 'nz': nz,
        'length': length, 'width': width, 'height': height,
        'x': x, 'y': y, 'z': z, 'dz': dz
    }


def create_flow_animation(config_name, output_name, grid):
    """Create animation of flow field evolution for a specific configuration."""
    
    # Find all snapshots for this configuration
    pattern = f'output/snapshots/{config_name}_*.bin'
    files = sorted(glob.glob(pattern))
    
    if len(files) == 0:
        print(f"  No snapshots found for {config_name}")
        return False
    
    print(f"  Creating animation from {len(files)} snapshots...")
    
    # Read first snapshot to set up plot
    data = read_snapshot(files[0])
    nx, ny, nz = data['nx'], data['ny'], data['nz']
    
    # Calculate speed
    speed = np.sqrt(data['u']**2 + data['v']**2 + data['w']**2)
    
    # Find runner center for slice
    solid_idx = np.where(data['solid'])
    if len(solid_idx[0]) > 0:
        runner_j = int(np.mean(solid_idx[1]))
        runner_i = int(np.mean(solid_idx[0]))
    else:
        runner_j = ny // 2
        runner_i = nx // 2
    
    # Find head height index
    head_k = np.searchsorted(grid['z'], 1.84)
    head_k = min(head_k, nz - 1)
    
    # Create figure with 2 views
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Side view (X-Z plane)
    ax1 = axes[0]
    speed_xz = speed[:, runner_j, :].T
    X_xz, Z_xz = np.meshgrid(grid['x'], grid['z'])
    
    vmax = 10.0  # Fixed scale for consistency
    im1 = ax1.pcolormesh(X_xz, Z_xz, speed_xz, cmap='viridis', vmin=0, vmax=vmax, shading='auto')
    
    # Add geometry outline
    mask_xz = data['solid'][:, runner_j, :].T
    ax1.contour(X_xz, Z_xz, mask_xz.astype(float), levels=[0.5], colors='white', linewidths=2)
    
    ax1.set_xlabel('X (m) - Flow direction →')
    ax1.set_ylabel('Z (m) - Height')
    ax1.set_title('Side View')
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1, label='Speed (m/s)')
    
    # Top view at head height (X-Y plane)
    ax2 = axes[1]
    speed_xy = speed[:, :, head_k].T
    X_xy, Y_xy = np.meshgrid(grid['x'], grid['y'])
    
    im2 = ax2.pcolormesh(X_xy, Y_xy, speed_xy, cmap='viridis', vmin=0, vmax=vmax, shading='auto')
    
    mask_xy = data['solid'][:, :, head_k].T
    ax2.contour(X_xy, Y_xy, mask_xy.astype(float), levels=[0.5], colors='white', linewidths=2)
    
    ax2.set_xlabel('X (m) - Flow direction →')
    ax2.set_ylabel('Y (m) - Width')
    ax2.set_title(f'Top View (Z={grid["z"][head_k]:.2f}m)')
    ax2.set_aspect('equal')
    plt.colorbar(im2, ax=ax2, label='Speed (m/s)')
    
    title = fig.suptitle(f'{config_name.replace("_", " ").title()} - Frame 1/{len(files)}', fontsize=14)
    
    plt.tight_layout()
    
    def update(frame):
        data = read_snapshot(files[frame])
        speed = np.sqrt(data['u']**2 + data['v']**2 + data['w']**2)
        
        # Update side view
        speed_xz = speed[:, runner_j, :].T
        im1.set_array(speed_xz.ravel())
        
        # Update top view
        speed_xy = speed[:, :, head_k].T
        im2.set_array(speed_xy.ravel())
        
        title.set_text(f'{config_name.replace("_", " ").title()} - Frame {frame+1}/{len(files)}')
        
        return [im1, im2, title]
    
    anim = FuncAnimation(fig, update, frames=len(files), interval=200, blit=False)
    
    # Save animation
    writer = PillowWriter(fps=5)
    anim.save(f'output/{output_name}.gif', writer=writer, dpi=100)
    plt.close(fig)
    
    print(f"  Saved: output/{output_name}.gif")
    return True


def create_drag_comparison_plot():
    """Create comprehensive drag comparison plot."""
    
    print("Creating drag comparison plot...")
    
    # Read results
    results = pd.read_csv('output/visor_results.csv')
    
    # Separate baseline and visor data
    baseline_row = results[results['angle_degrees'] == 'baseline']
    baseline_drag = float(baseline_row['drag_N'].iloc[0])
    
    visor_data = results[results['angle_degrees'] != 'baseline'].copy()
    visor_data['angle_degrees'] = visor_data['angle_degrees'].astype(float)
    visor_data = visor_data.sort_values('angle_degrees')
    
    angles = visor_data['angle_degrees'].values
    drags = visor_data['drag_N'].values
    drag_diffs = visor_data['drag_diff_N'].values
    drag_pcts = visor_data['drag_diff_percent'].values
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 10))
    
    # Plot 1: Absolute drag vs angle
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.axhline(y=baseline_drag, color='gray', linestyle='--', linewidth=2, 
                label=f'No visor: {baseline_drag:.2f} N')
    ax1.plot(angles, drags, 'bo-', markersize=8, linewidth=2, label='With visor')
    ax1.axvline(x=0, color='lightgray', linestyle=':', alpha=0.5)
    ax1.set_xlabel('Visor Angle (degrees)', fontsize=12)
    ax1.set_ylabel('Drag Force (N)', fontsize=12)
    ax1.set_title('Drag Force vs Visor Angle', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Highlight optimal angle
    opt_idx = np.argmin(np.abs(drags))
    ax1.scatter([angles[opt_idx]], [drags[opt_idx]], color='green', s=200, 
                zorder=5, marker='*', label=f'Optimal: {angles[opt_idx]:.0f}°')
    
    # Plot 2: Drag difference from baseline
    ax2 = fig.add_subplot(2, 2, 2)
    colors = ['green' if d < 0 else 'red' for d in drag_diffs]
    bars = ax2.bar(angles, drag_diffs, color=colors, alpha=0.7, width=4)
    ax2.axhline(y=0, color='black', linewidth=1)
    ax2.axvline(x=0, color='lightgray', linestyle=':', alpha=0.5)
    ax2.set_xlabel('Visor Angle (degrees)', fontsize=12)
    ax2.set_ylabel('Drag Change from Baseline (N)', fontsize=12)
    ax2.set_title('Visor Drag Impact', fontsize=14)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Percent change from baseline
    ax3 = fig.add_subplot(2, 2, 3)
    colors = ['green' if d < 0 else 'red' for d in drag_pcts]
    ax3.bar(angles, drag_pcts, color=colors, alpha=0.7, width=4)
    ax3.axhline(y=0, color='black', linewidth=1)
    ax3.axvline(x=0, color='lightgray', linestyle=':', alpha=0.5)
    ax3.set_xlabel('Visor Angle (degrees)', fontsize=12)
    ax3.set_ylabel('Drag Change (%)', fontsize=12)
    ax3.set_title('Percentage Drag Change from Baseline', fontsize=14)
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Polar plot showing drag coefficient variation
    ax4 = fig.add_subplot(2, 2, 4, projection='polar')
    # Convert angles to radians, shift so 0° is at top
    theta = np.radians(angles + 90)  # Shift so 0° points up
    
    # Normalize drag for polar plot (use absolute values)
    r = np.abs(drags) / np.max(np.abs(drags))
    
    ax4.plot(theta, r, 'b-', linewidth=2)
    ax4.scatter(theta, r, c=drags, cmap='RdYlGn_r', s=100, zorder=5)
    ax4.set_title('Drag Magnitude by Visor Angle\n(0° = horizontal visor)', fontsize=12)
    
    # Add angle labels
    ax4.set_xticks(np.radians([45, 90, 135, 180, 225, 270, 315, 360]))
    ax4.set_xticklabels(['-45°', '0°', '45°', '', '-45°', '0°', '45°', ''])
    
    fig.suptitle('Visor Aerodynamics Analysis', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    # Save plot in multiple formats
    for ext in ['png', 'eps', 'svg']:
        fig.savefig(f'output/drag_comparison.{ext}', dpi=150, bbox_inches='tight')
    
    plt.close(fig)
    print("  Saved: output/drag_comparison.{png,eps,svg}")
    
    # Print summary
    print("\n" + "="*60)
    print("DRAG ANALYSIS SUMMARY")
    print("="*60)
    print(f"Baseline (no visor): {baseline_drag:.4f} N")
    print(f"\nOptimal visor angle: {angles[opt_idx]:.0f}°")
    print(f"Drag at optimal: {drags[opt_idx]:.4f} N")
    print(f"Change from baseline: {drag_diffs[opt_idx]:.4f} N ({drag_pcts[opt_idx]:.2f}%)")
    print("\nAll configurations:")
    print("-"*60)
    for i, angle in enumerate(angles):
        sign = '+' if drag_diffs[i] > 0 else ''
        print(f"  {angle:+6.0f}°: {drags[i]:12.4f} N  ({sign}{drag_pcts[i]:.2f}%)")


def create_grid_resolution_plot():
    """Create plot showing non-uniform grid resolution."""
    
    print("Creating grid resolution plot...")
    
    grid = read_grid_info()
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: dz vs height
    ax1 = axes[0]
    ax1.plot(grid['dz'] * 1000, grid['z'], 'b-', linewidth=2)
    ax1.axhline(y=1.84, color='red', linestyle='--', label='Head height (1.84m)')
    ax1.axhspan(1.6, 2.1, alpha=0.2, color='green', label='Refinement zone')
    ax1.set_xlabel('Cell Size dz (mm)', fontsize=12)
    ax1.set_ylabel('Height z (m)', fontsize=12)
    ax1.set_title('Vertical Grid Resolution', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add visor resolution annotation
    head_idx = np.searchsorted(grid['z'], 1.84)
    dz_head = grid['dz'][head_idx] * 1000
    visor_cells = 5.0 / dz_head  # 5mm visor
    ax1.annotate(f'dz={dz_head:.2f}mm\nVisor spans {visor_cells:.1f} cells',
                xy=(dz_head, 1.84), xytext=(dz_head + 5, 1.5),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red')
    
    # Plot 2: z-coordinate distribution
    ax2 = axes[1]
    ax2.plot(range(1, len(grid['z'])+1), grid['z'], 'b-', linewidth=2)
    ax2.axhline(y=1.84, color='red', linestyle='--', label='Head height')
    ax2.axhspan(1.6, 2.1, alpha=0.2, color='green', label='Refinement zone')
    ax2.set_xlabel('Vertical Level Index', fontsize=12)
    ax2.set_ylabel('Height z (m)', fontsize=12)
    ax2.set_title('Vertical Level Distribution', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    fig.suptitle(f'Non-Uniform Grid: {grid["nx"]}×{grid["ny"]}×{grid["nz"]} cells', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    for ext in ['png', 'eps', 'svg']:
        fig.savefig(f'output/grid_resolution.{ext}', dpi=150, bbox_inches='tight')
    
    plt.close(fig)
    print("  Saved: output/grid_resolution.{png,eps,svg}")


def create_cross_section_plots(config_name, grid):
    """Create cross-section plots showing flow field with masked runner.
    
    Creates three views:
    - X-Z plane (side view, through runner center)
    - Y-Z plane (front view, through runner center)
    - X-Y plane (top view, at visor/head level)
    """
    
    # Find the final snapshot for this configuration
    pattern = f'output/snapshots/{config_name}_*.bin'
    files = sorted(glob.glob(pattern))
    
    if len(files) == 0:
        print(f"  No snapshots found for {config_name}")
        return False
    
    # Use the last snapshot (final state)
    final_file = files[-1]
    print(f"  Creating cross-sections for {config_name} from {Path(final_file).name}...")
    
    data = read_snapshot(final_file)
    nx, ny, nz = data['nx'], data['ny'], data['nz']
    
    # Calculate velocity magnitude
    speed = np.sqrt(data['u']**2 + data['v']**2 + data['w']**2)
    
    # Find runner center indices
    solid_idx = np.where(data['solid'])
    if len(solid_idx[0]) > 0:
        runner_i = int(np.mean(solid_idx[0]))  # X index (center of runner)
        runner_j = int(np.mean(solid_idx[1]))  # Y index (center of runner)
    else:
        runner_i = nx // 2
        runner_j = ny // 2
    
    # Find head/visor height index (~1.84m)
    head_k = np.searchsorted(grid['z'], 1.84)
    head_k = min(head_k, nz - 1)
    
    # Color settings
    vmax = 12.0  # Max speed for colormap (slightly above inlet 8 m/s)
    cmap = 'viridis'
    
    # Create figure with 3 subplots
    fig = plt.figure(figsize=(18, 6))
    
    # =========================================================================
    # Plot 1: X-Z plane (side view) through runner center (y = runner_j)
    # =========================================================================
    ax1 = fig.add_subplot(1, 3, 1)
    
    # Extract slice
    speed_xz = speed[:, runner_j, :].T
    solid_xz = data['solid'][:, runner_j, :].T
    u_xz = data['u'][:, runner_j, :].T
    w_xz = data['w'][:, runner_j, :].T
    
    # Mask solid regions
    speed_xz_masked = np.ma.array(speed_xz, mask=solid_xz)
    
    # Create mesh
    X_xz, Z_xz = np.meshgrid(grid['x'], grid['z'])
    
    # Plot speed field
    im1 = ax1.pcolormesh(X_xz, Z_xz, speed_xz_masked, cmap=cmap, vmin=0, vmax=vmax, shading='auto')
    
    # Overlay solid geometry in gray
    ax1.contourf(X_xz, Z_xz, solid_xz.astype(float), levels=[0.5, 1.5], colors=['dimgray'], alpha=1.0)
    ax1.contour(X_xz, Z_xz, solid_xz.astype(float), levels=[0.5], colors='black', linewidths=1.5)
    
    # Add velocity vectors (subsampled)
    skip_x, skip_z = max(1, nx // 20), max(1, nz // 25)
    X_sub = X_xz[::skip_z, ::skip_x]
    Z_sub = Z_xz[::skip_z, ::skip_x]
    U_sub = u_xz[::skip_z, ::skip_x]
    W_sub = w_xz[::skip_z, ::skip_x]
    solid_sub = solid_xz[::skip_z, ::skip_x]
    
    # Mask arrows in solid regions
    U_sub = np.ma.array(U_sub, mask=solid_sub)
    W_sub = np.ma.array(W_sub, mask=solid_sub)
    
    ax1.quiver(X_sub, Z_sub, U_sub, W_sub, color='white', alpha=0.6, scale=100, width=0.003)
    
    ax1.set_xlabel('X (m) - Flow direction →', fontsize=11)
    ax1.set_ylabel('Z (m) - Height', fontsize=11)
    ax1.set_title(f'Side View (X-Z plane, Y={grid["y"][runner_j]:.2f}m)', fontsize=12, fontweight='bold')
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1, label='Speed (m/s)', shrink=0.8)
    
    # =========================================================================
    # Plot 2: Y-Z plane (front view) through runner center (x = runner_i)
    # =========================================================================
    ax2 = fig.add_subplot(1, 3, 2)
    
    # Extract slice
    speed_yz = speed[runner_i, :, :].T
    solid_yz = data['solid'][runner_i, :, :].T
    v_yz = data['v'][runner_i, :, :].T
    w_yz = data['w'][runner_i, :, :].T
    
    # Mask solid regions
    speed_yz_masked = np.ma.array(speed_yz, mask=solid_yz)
    
    # Create mesh
    Y_yz, Z_yz = np.meshgrid(grid['y'], grid['z'])
    
    # Plot speed field
    im2 = ax2.pcolormesh(Y_yz, Z_yz, speed_yz_masked, cmap=cmap, vmin=0, vmax=vmax, shading='auto')
    
    # Overlay solid geometry
    ax2.contourf(Y_yz, Z_yz, solid_yz.astype(float), levels=[0.5, 1.5], colors=['dimgray'], alpha=1.0)
    ax2.contour(Y_yz, Z_yz, solid_yz.astype(float), levels=[0.5], colors='black', linewidths=1.5)
    
    # Add velocity vectors (subsampled)
    skip_y, skip_z = max(1, ny // 15), max(1, nz // 25)
    Y_sub = Y_yz[::skip_z, ::skip_y]
    Z_sub = Z_yz[::skip_z, ::skip_y]
    V_sub = v_yz[::skip_z, ::skip_y]
    W_sub = w_yz[::skip_z, ::skip_y]
    solid_sub = solid_yz[::skip_z, ::skip_y]
    
    V_sub = np.ma.array(V_sub, mask=solid_sub)
    W_sub = np.ma.array(W_sub, mask=solid_sub)
    
    ax2.quiver(Y_sub, Z_sub, V_sub, W_sub, color='white', alpha=0.6, scale=50, width=0.003)
    
    ax2.set_xlabel('Y (m) - Width', fontsize=11)
    ax2.set_ylabel('Z (m) - Height', fontsize=11)
    ax2.set_title(f'Front View (Y-Z plane, X={grid["x"][runner_i]:.2f}m)', fontsize=12, fontweight='bold')
    ax2.set_aspect('equal')
    plt.colorbar(im2, ax=ax2, label='Speed (m/s)', shrink=0.8)
    
    # =========================================================================
    # Plot 3: X-Y plane (top view) at visor/head level (z = head_k)
    # =========================================================================
    ax3 = fig.add_subplot(1, 3, 3)
    
    # Extract slice
    speed_xy = speed[:, :, head_k].T
    solid_xy = data['solid'][:, :, head_k].T
    u_xy = data['u'][:, :, head_k].T
    v_xy = data['v'][:, :, head_k].T
    
    # Mask solid regions
    speed_xy_masked = np.ma.array(speed_xy, mask=solid_xy)
    
    # Create mesh
    X_xy, Y_xy = np.meshgrid(grid['x'], grid['y'])
    
    # Plot speed field
    im3 = ax3.pcolormesh(X_xy, Y_xy, speed_xy_masked, cmap=cmap, vmin=0, vmax=vmax, shading='auto')
    
    # Overlay solid geometry
    ax3.contourf(X_xy, Y_xy, solid_xy.astype(float), levels=[0.5, 1.5], colors=['dimgray'], alpha=1.0)
    ax3.contour(X_xy, Y_xy, solid_xy.astype(float), levels=[0.5], colors='black', linewidths=1.5)
    
    # Add velocity vectors (subsampled)
    skip_x, skip_y = max(1, nx // 20), max(1, ny // 15)
    X_sub = X_xy[::skip_y, ::skip_x]
    Y_sub = Y_xy[::skip_y, ::skip_x]
    U_sub = u_xy[::skip_y, ::skip_x]
    V_sub = v_xy[::skip_y, ::skip_x]
    solid_sub = solid_xy[::skip_y, ::skip_x]
    
    U_sub = np.ma.array(U_sub, mask=solid_sub)
    V_sub = np.ma.array(V_sub, mask=solid_sub)
    
    ax3.quiver(X_sub, Y_sub, U_sub, V_sub, color='white', alpha=0.6, scale=100, width=0.003)
    
    ax3.set_xlabel('X (m) - Flow direction →', fontsize=11)
    ax3.set_ylabel('Y (m) - Width', fontsize=11)
    ax3.set_title(f'Top View (X-Y plane, Z={grid["z"][head_k]:.2f}m - Visor Level)', fontsize=12, fontweight='bold')
    ax3.set_aspect('equal')
    plt.colorbar(im3, ax=ax3, label='Speed (m/s)', shrink=0.8)
    
    # Main title
    title_name = config_name.replace('_', ' ').title()
    if 'visor' in config_name.lower():
        try:
            angle_code = int(config_name.split('_')[1])
            angle = angle_code - 100
            title_name = f'Visor at {angle:+d}°'
        except:
            pass
    
    fig.suptitle(f'Flow Field Cross-Sections: {title_name}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save in multiple formats
    output_base = f'output/cross_sections_{config_name}'
    for ext in ['png', 'eps', 'svg']:
        fig.savefig(f'{output_base}.{ext}', dpi=150, bbox_inches='tight')
    
    plt.close(fig)
    print(f"    Saved: {output_base}.{{png,eps,svg}}")
    return True


def main():
    """Main function to generate all plots and animations."""
    
    print("="*60)
    print("VISOR AERODYNAMICS - POST-PROCESSING")
    print("="*60)
    print()
    
    # Check for output directory
    if not Path('output').exists():
        print("ERROR: output/ directory not found. Run the Fortran simulation first.")
        return
    
    # Read grid info
    try:
        grid = read_grid_info()
        print(f"Grid: {grid['nx']}x{grid['ny']}x{grid['nz']}")
        print(f"Domain: {grid['length']}m x {grid['width']}m x {grid['height']}m")
        print()
    except Exception as e:
        print(f"ERROR reading grid info: {e}")
        return
    
    # Create grid resolution plot
    create_grid_resolution_plot()
    print()
    
    # Create drag comparison plot
    try:
        create_drag_comparison_plot()
    except Exception as e:
        print(f"ERROR creating drag plot: {e}")
    print()
    
    # Create animations for each configuration
    print("Creating flow field animations...")
    
    # Baseline animation
    create_flow_animation('baseline', 'animation_baseline', grid)
    
    # Find all visor configurations
    snapshot_files = glob.glob('output/snapshots/visor_*_0001.bin')
    configs = set()
    for f in snapshot_files:
        # Extract config name (e.g., visor_0055 for -45 degrees, visor_0145 for +45)
        parts = Path(f).stem.split('_')
        if len(parts) >= 2:
            config = f"{parts[0]}_{parts[1]}"
            configs.add(config)
    
    for config in sorted(configs):
        # Decode angle from filename
        try:
            angle_code = int(config.split('_')[1])
            angle = angle_code - 100  # Decode: 55 -> -45, 145 -> +45
            output_name = f'animation_visor_{angle:+d}deg'
            create_flow_animation(config, output_name, grid)
        except Exception as e:
            print(f"  Error processing {config}: {e}")
    
    print()
    
    # Create cross-section plots for final flow states
    print("Creating cross-section plots...")
    
    # Baseline cross-sections
    create_cross_section_plots('baseline', grid)
    
    # Cross-sections for each visor configuration
    for config in sorted(configs):
        try:
            create_cross_section_plots(config, grid)
        except Exception as e:
            print(f"  Error processing {config}: {e}")
    
    print()
    print("="*60)
    print("POST-PROCESSING COMPLETE")
    print("="*60)
    print("\nGenerated files in output/:")
    print("  - drag_comparison.{png,eps,svg}: Drag analysis plots")
    print("  - grid_resolution.{png,eps,svg}: Grid resolution visualization")
    print("  - cross_sections_*.{png,eps,svg}: Flow field cross-sections")
    print("  - animation_*.gif: Flow field animations")


if __name__ == "__main__":
    main()
