#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import argparse
from tqdm import tqdm
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle

def create_chip_layout(ref_file, chip_type='small'):
    """
    Read reference file and create chip layout matrix
    
    Parameters:
    ref_file: Reference file for synthetic chip
    chip_type: Chip type 'small'(540x635) or 'large'(1080x636)
    
    Returns:
    sequences: List of all synthesized sequences
    seq_to_pos: Mapping from sequence to chip position
    """
    print(f"Reading reference file: {ref_file}")
    
    # Set chip dimensions
    if chip_type == 'small':
        width, height = 540, 635  # width=X轴方向, height=Y轴方向
        expected_lines = 342583
    else:  # large
        width, height = 1080, 636  # width=X轴方向, height=Y轴方向
        expected_lines = 686880
    
    # Read reference file
    with open(ref_file, 'r') as f:
        ref_sequences = [line.strip() for line in f]
    
    print(f"Read {len(ref_sequences)} sequence lines")
    
    # Ensure sequence count meets expectations
    if len(ref_sequences) < expected_lines:
        print(f"Warning: Sequence count in file is less than expected ({len(ref_sequences)} < {expected_lines})")
        # Fill with empty sequences to expected lines
        ref_sequences.extend(["0" * 150] * (expected_lines - len(ref_sequences)))
    
    # Extract substring from positions 16-120 for each sequence
    sequences = [seq[15:120] if len(seq) >= 120 else seq for seq in ref_sequences]
    
    # Create mapping from sequence to chip position
    seq_to_pos = {}
    
    index = 0
    # 先遍历X轴（从左到右），再遍历Y轴（从下到上）
    for x in range(width):
        for y in range(height):
            if index < len(sequences):
                if sequences[index] != "0" * 105:  # Exclude filling sequences
                    # 注意：排列是从左到右，每列从下到上
                    seq_to_pos[sequences[index]] = (x, y)
                index += 1
    
    print(f"Created {len(seq_to_pos)} sequence-to-position mappings")
    return sequences, seq_to_pos

def map_depth_to_chip(ngs_file, seq_to_pos, output_dir):
    """
    Map NGS depth data to chip position and generate heatmap
    
    Parameters:
    ngs_file: NGS sequencing result file
    seq_to_pos: Mapping from sequence to chip position
    output_dir: Output directory
    
    Returns:
    chip_matrix: Chip depth matrix
    mapped_count: Count of successfully mapped sequences
    """
    print(f"Reading NGS file: {ngs_file}")
    
    # Read NGS sequencing results
    try:
        # First try to read as TSV format
        df = pd.read_csv(ngs_file, sep='\t')
        print("Successfully read file in TSV format")
    except Exception as e1:
        try:
            # Then try as Excel file
            if ngs_file.endswith('.xls') or ngs_file.endswith('.xlsx'):
                df = pd.read_excel(ngs_file)
            else:
                raise Exception("Cannot identify file format")
        except Exception as e2:
            print(f"Error reading file: {e1}\nSecond attempt error: {e2}")
            return None, 0
    
    # Check required columns
    if 'Oligo seq' not in df.columns or 'Depth' not in df.columns:
        print(f"Error: Required columns not found. Available columns: {', '.join(df.columns)}")
        return None, 0
    
    # Get chip dimensions
    max_x = max([pos[0] for pos in seq_to_pos.values()]) + 1
    max_y = max([pos[1] for pos in seq_to_pos.values()]) + 1
    
    # Create chip depth matrix, initial value is NaN (indicating no data)
    chip_matrix = np.full((max_y, max_x), np.nan)
    
    # Map depth to chip position
    mapped_count = 0
    not_found = 0
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Mapping depth data"):
        oligo_seq = row['Oligo seq']
        depth = row['Depth']
        
        if oligo_seq in seq_to_pos:
            x, y = seq_to_pos[oligo_seq]
            # Note: Matrix index is (row, column), corresponding to (y, x)
            chip_matrix[y, x] = depth
            mapped_count += 1
        else:
            not_found += 1
    
    mapping_rate = mapped_count / len(df) * 100 if len(df) > 0 else 0
    print(f"Total sequences: {len(df)}")
    print(f"Successfully mapped: {mapped_count} ({mapping_rate:.2f}%)")
    print(f"Mapping not found: {not_found}")
    
    # Save chip depth matrix to CSV
    df_matrix = pd.DataFrame(chip_matrix)
    matrix_file = os.path.join(output_dir, 'chip_depth_matrix.csv')
    df_matrix.to_csv(matrix_file)
    print(f"Chip depth matrix saved to: {matrix_file}")
    
    return chip_matrix, mapped_count

def generate_heatmaps(chip_matrix, ngs_file, output_dir):
    """
    Generate depth distribution heatmaps
    
    Parameters:
    chip_matrix: Chip depth matrix
    ngs_file: NGS file name, used to generate output file name
    output_dir: Output directory
    """
    if chip_matrix is None:
        print("Error: Chip depth matrix is empty, unable to generate heatmap")
        return
    
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract file name (without extension) for output file naming
    file_name = os.path.splitext(os.path.basename(ngs_file))[0]
    
    # Create mask for valid data points
    mask = np.isnan(chip_matrix)
    
    # Calculate mean depth for thresholds
    mean_depth = np.nanmean(chip_matrix)
    print(f"Mean sequencing depth: {mean_depth:.2f}")
    
    # 1. Standard heatmap showing all depth values
    plt.figure(figsize=(14, 12))
    
    # Create custom colormap that highlights low values
    colors = plt.cm.viridis(np.linspace(0, 1, 256))
    # Make low values more visible with brighter colors
    colors[:64] = plt.cm.Reds(np.linspace(0.6, 1.0, 64))
    custom_cmap = mcolors.ListedColormap(colors)
    
    # Plot heatmap with no spaces between cells and clear edge color
    ax = sns.heatmap(chip_matrix, mask=mask, cmap=custom_cmap, 
                    cbar_kws={'label': 'Depth'}, linewidths=0, 
                    square=True, xticklabels=False, yticklabels=False)
    
    # Add a rectangle border around the chip
    height, width = chip_matrix.shape
    rect = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    ax.set_title('Sequencing Depth Distribution on Chip', fontsize=14)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Save standard heatmap
    standard_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_standard.png')
    plt.savefig(standard_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Standard heatmap saved to: {standard_heatmap_path}")
    
    # 2. Low coverage heatmap (points <0.5x mean depth highlighted)
    plt.figure(figsize=(14, 12))
    
    # Create a custom colormap for low coverage highlighting
    # Use a diverging colormap to emphasize low coverage areas
    low_coverage_cmap = plt.cm.RdYlBu_r
    
    # Set vmax to limit the color range
    vmax = mean_depth * 2
    vmin = 0
    
    # Create normalized colorscale
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Plot heatmap with custom color mapping and show ticks on top and right edge
    ax = sns.heatmap(chip_matrix, mask=mask, cmap=low_coverage_cmap, 
                    norm=norm, cbar_kws={'label': 'Depth'}, 
                    linewidths=0, square=True)
    
    # Configure tick spacing and labels
    tick_interval = max(width, height) // 10  # Add a tick every ~10% of dimension
    
    # Set x-axis ticks with even spacing
    x_ticks = np.arange(0, width, tick_interval)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks)
    
    # Set y-axis ticks with even spacing
    y_ticks = np.arange(0, height, tick_interval)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks)
    
    # Add the top and right spines with ticks
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    
    # Add a rectangle border around the chip
    rect = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    # Set threshold for low depth as 0.5x mean depth
    low_threshold = mean_depth * 0.5
    
    # Create a binary mask of low depth areas
    low_depth_mask = np.zeros_like(chip_matrix)
    low_depth_mask[~mask & (chip_matrix < low_threshold)] = 1
    
    # Plot contours around low depth regions
    if np.any(low_depth_mask):
        try:
            cs = ax.contour(low_depth_mask, levels=[0.5], colors=['black'], linewidths=1.5, 
                          extent=(-0.5, width-0.5, -0.5, height-0.5))
            # Add contour labels
            plt.clabel(cs, inline=1, fontsize=10, fmt='Low depth region')
        except:
            print("Unable to generate contours for low depth regions")
    
    ax.set_title('Depth Distribution with Emphasis on Low Coverage', fontsize=14)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add annotation for low depth threshold
    plt.annotate(f'Low depth threshold: {low_threshold:.2f} (0.5x mean)', 
                xy=(0.05, 0.95), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    
    # Save low coverage emphasis heatmap
    low_emph_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_low_emph.png')
    plt.savefig(low_emph_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Low coverage emphasis heatmap saved to: {low_emph_heatmap_path}")
    
    # 3. Focused low coverage heatmap (highlighting only low depth regions)
    plt.figure(figsize=(14, 12))
    
    # Create a binary mask: True for non-low depths and NaN values
    focused_mask = np.logical_or(mask, chip_matrix >= low_threshold)
    
    # Create a custom colormap for low coverage only - Reds
    low_only_cmap = plt.cm.Reds_r
    
    # Plot heatmap showing only low coverage areas
    ax = sns.heatmap(chip_matrix, mask=focused_mask, cmap=low_only_cmap,
                    vmin=0, vmax=low_threshold,
                    cbar_kws={'label': f'Depth < {low_threshold:.2f}'}, 
                    linewidths=0, square=True, xticklabels=False, yticklabels=False)
    
    # Add a rectangle border around the chip
    height, width = chip_matrix.shape
    rect = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    # Add grid lines to show chip coordinates (sparse grid)
    grid_spacing = max(width, height) // 10  # Add grid line every ~10% of dimension
    for x in range(0, width, grid_spacing):
        ax.axvline(x, color='grey', linestyle='-', linewidth=0.5, alpha=0.3)
    for y in range(0, height, grid_spacing):
        ax.axhline(y, color='grey', linestyle='-', linewidth=0.5, alpha=0.3)
    
    ax.set_title('Low Coverage Regions (<50% Mean Depth)', fontsize=14)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Calculate percentage of low depth regions
    low_depth_count = np.sum(~focused_mask)
    total_valid_points = np.sum(~mask)
    low_depth_percentage = (low_depth_count / total_valid_points) * 100 if total_valid_points > 0 else 0
    
    # Add annotation for low depth percentage
    plt.annotate(f'Low depth regions: {low_depth_percentage:.2f}% of chip', 
                xy=(0.05, 0.05), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    
    # Save focused low coverage heatmap
    low_focus_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_low_focus.png')
    plt.savefig(low_focus_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Focused low coverage heatmap saved to: {low_focus_heatmap_path}")
    
    # 4. High coverage heatmap (points >2x mean depth in blue)
    plt.figure(figsize=(14, 12))
    
    # Create a copy of the matrix for high coverage
    high_matrix = np.copy(chip_matrix)
    # Mask all values ≤2x mean
    high_mask = np.logical_or(mask, high_matrix <= 2.0 * mean_depth)
    
    # Create a custom colormap with only blue
    colors = [(0, 0, 1, 1)]  # Blue with full opacity
    high_cmap = mcolors.ListedColormap(colors)
    
    ax = sns.heatmap(high_matrix, mask=high_mask, cmap=high_cmap, 
                    cbar_kws={'label': f'Depth > 2x Mean ({mean_depth*2.0:.2f})'},
                    linewidths=0, square=True, xticklabels=False, yticklabels=False)
    
    # Add a rectangle border around the chip
    height, width = chip_matrix.shape
    rect = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    ax.set_title('High Coverage Regions (> 2x Mean Depth)', fontsize=14)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Save high coverage heatmap
    high_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_high_coverage.png')
    plt.savefig(high_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"High coverage heatmap saved to: {high_heatmap_path}")

def analyze_position_bias(chip_matrix, ngs_file, output_dir):
    """
    Analyze position bias, examine relationship between depth and chip position
    
    Parameters:
    chip_matrix: Chip depth matrix
    ngs_file: NGS file name, used to generate output file name
    output_dir: Output directory
    """
    if chip_matrix is None:
        print("Error: Chip depth matrix is empty, unable to analyze position bias")
        return
    
    # Extract file name (without extension) for output file naming
    file_name = os.path.splitext(os.path.basename(ngs_file))[0]
    
    # Create position data
    positions = []
    depths = []
    
    for y in range(chip_matrix.shape[0]):
        for x in range(chip_matrix.shape[1]):
            if not np.isnan(chip_matrix[y, x]):
                positions.append((x, y))
                depths.append(chip_matrix[y, x])
    
    if not positions:
        print("Not enough data for position bias analysis")
        return
    
    # Extract X and Y coordinates
    x_coords = [p[0] for p in positions]
    y_coords = [p[1] for p in positions]
    
    # Create X position-based analysis
    plt.figure(figsize=(12, 6))
    plt.scatter(x_coords, depths, alpha=0.5, s=10)
    plt.title('Depth vs X Position')
    plt.xlabel('X Position')
    plt.ylabel('Depth')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add trend line
    try:
        z = np.polyfit(x_coords, depths, 1)
        p = np.poly1d(z)
        x_range = np.linspace(min(x_coords), max(x_coords), 100)
        plt.plot(x_range, p(x_range), 'r--')
        corr = np.corrcoef(x_coords, depths)[0, 1]
        plt.annotate(f'Correlation: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"Unable to calculate X trend line: {e}")
    
    # Save X position analysis figure
    x_analysis_path = os.path.join(output_dir, f'{file_name}_x_position_analysis.png')
    plt.savefig(x_analysis_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create Y position-based analysis
    plt.figure(figsize=(12, 6))
    plt.scatter(y_coords, depths, alpha=0.5, s=10)
    plt.title('Depth vs Y Position')
    plt.xlabel('Y Position')
    plt.ylabel('Depth')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add trend line
    try:
        z = np.polyfit(y_coords, depths, 1)
        p = np.poly1d(z)
        y_range = np.linspace(min(y_coords), max(y_coords), 100)
        plt.plot(y_range, p(y_range), 'r--')
        corr = np.corrcoef(y_coords, depths)[0, 1]
        plt.annotate(f'Correlation: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"Unable to calculate Y trend line: {e}")
    
    # Save Y position analysis figure
    y_analysis_path = os.path.join(output_dir, f'{file_name}_y_position_analysis.png')
    plt.savefig(y_analysis_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Perform quartile and central region statistical analysis
    max_x = max(x_coords)
    max_y = max(y_coords)
    mid_x = max_x / 2
    mid_y = max_y / 2
    
    # Define regions
    regions = {
        "Upper Left": [(x, y) for x, y in positions if x < mid_x and y > mid_y],
        "Upper Right": [(x, y) for x, y in positions if x >= mid_x and y > mid_y],
        "Lower Left": [(x, y) for x, y in positions if x < mid_x and y <= mid_y],
        "Lower Right": [(x, y) for x, y in positions if x >= mid_x and y <= mid_y],
        "Center": [(x, y) for x, y in positions if mid_x/2 <= x <= mid_x*1.5 and mid_y/2 <= y <= mid_y*1.5]
    }
    
    # Calculate statistics for each region
    region_stats = {}
    for region, pos in regions.items():
        if pos:
            region_depths = [depths[positions.index(p)] for p in pos if p in positions]
            if region_depths:
                region_stats[region] = {
                    "Count": len(region_depths),
                    "Mean Depth": np.mean(region_depths),
                    "Median Depth": np.median(region_depths),
                    "Std Dev": np.std(region_depths),
                    "Min": np.min(region_depths),
                    "Max": np.max(region_depths)
                }
    
    # Save region statistics
    region_stats_df = pd.DataFrame(region_stats).T
    region_stats_path = os.path.join(output_dir, f'{file_name}_region_statistics.csv')
    region_stats_df.to_csv(region_stats_path)
    print(f"Region statistics saved to: {region_stats_path}")
    
    # Create region statistics bar chart
    if region_stats:
        regions = list(region_stats.keys())
        means = [region_stats[r]["Mean Depth"] for r in regions]
        
        plt.figure(figsize=(10, 6))
        bars = plt.bar(regions, means)
        
        # Add data labels
        for bar, mean in zip(bars, means):
            plt.text(bar.get_x() + bar.get_width()/2, mean + 5, 
                    f'{mean:.1f}', 
                    ha='center', va='bottom')
        
        plt.title('Average Sequencing Depth in Different Regions')
        plt.ylabel('Mean Depth')
        plt.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Save region statistics figure
        region_plot_path = os.path.join(output_dir, f'{file_name}_region_depth_analysis.png')
        plt.savefig(region_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Region analysis plot saved to: {region_plot_path}")

def main():
    parser = argparse.ArgumentParser(description='Map NGS depth data to chip physical location and generate heatmaps')
    parser.add_argument('ref_file', help='Path to the reference file for synthetic chip')
    parser.add_argument('ngs_file', help='Path to the NGS sequencing result file')
    parser.add_argument('--output_dir', '-o', default='chip_analysis_results', 
                        help='Output directory name')
    parser.add_argument('--chip_type', '-t', choices=['small', 'large'], default='small',
                        help='Chip type: small(540x635) or large(1080x636)')
    
    args = parser.parse_args()
    
    # Create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Create chip layout
    sequences, seq_to_pos = create_chip_layout(args.ref_file, args.chip_type)
    
    # Map depth to chip position
    chip_matrix, mapped_count = map_depth_to_chip(args.ngs_file, seq_to_pos, args.output_dir)
    
    if mapped_count > 0:
        # Generate heatmaps
        generate_heatmaps(chip_matrix, args.ngs_file, args.output_dir)
        
        # Analyze position bias
        analyze_position_bias(chip_matrix, args.ngs_file, args.output_dir)
        
        print(f"\nAnalysis complete! Results saved to '{args.output_dir}' directory")
    else:
        print("No sequences were successfully mapped, cannot generate heatmaps")

if __name__ == "__main__":
    main() 