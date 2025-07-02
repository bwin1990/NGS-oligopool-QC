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

def map_indel_to_chip(ngs_file, seq_to_pos, output_dir):
    """
    Map NGS indel data to chip position and generate heatmap
    
    Parameters:
    ngs_file: NGS sequencing result file
    seq_to_pos: Mapping from sequence to chip position
    output_dir: Output directory
    
    Returns:
    small_indel_matrix: Chip matrix for small indels (<3bp del)
    large_indel_matrix: Chip matrix for large indels (>=3bp del)
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
            return None, None, 0
    
    # Check required columns
    required_columns = ['Oligo seq', 'Mean indel(<3bp del)', 'Mean indel']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Error: Required columns not found: {', '.join(missing_columns)}")
        print(f"Available columns: {', '.join(df.columns)}")
        return None, None, 0
    
    # Get chip dimensions
    max_x = max([pos[0] for pos in seq_to_pos.values()]) + 1
    max_y = max([pos[1] for pos in seq_to_pos.values()]) + 1
    
    # Create chip indel matrices, initial value is NaN (indicating no data)
    small_indel_matrix = np.full((max_y, max_x), np.nan)
    large_indel_matrix = np.full((max_y, max_x), np.nan)
    
    # Map indel data to chip position
    mapped_count = 0
    not_found = 0
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Mapping indel data"):
        oligo_seq = row['Oligo seq']
        small_indel = row['Mean indel(<3bp del)']
        mean_indel = row['Mean indel']
        
        # Calculate large indel (>=3bp del)
        large_indel = mean_indel - small_indel
        
        if oligo_seq in seq_to_pos:
            x, y = seq_to_pos[oligo_seq]
            # Note: Matrix index is (row, column), corresponding to (y, x)
            small_indel_matrix[y, x] = small_indel
            large_indel_matrix[y, x] = large_indel
            mapped_count += 1
        else:
            not_found += 1
    
    mapping_rate = mapped_count / len(df) * 100 if len(df) > 0 else 0
    print(f"Total sequences: {len(df)}")
    print(f"Successfully mapped: {mapped_count} ({mapping_rate:.2f}%)")
    print(f"Mapping not found: {not_found}")
    
    # Save chip indel matrices to CSV
    small_indel_df = pd.DataFrame(small_indel_matrix)
    large_indel_df = pd.DataFrame(large_indel_matrix)
    
    small_indel_file = os.path.join(output_dir, 'chip_small_indel_matrix.csv')
    large_indel_file = os.path.join(output_dir, 'chip_large_indel_matrix.csv')
    
    small_indel_df.to_csv(small_indel_file)
    large_indel_df.to_csv(large_indel_file)
    
    print(f"Small indel matrix saved to: {small_indel_file}")
    print(f"Large indel matrix saved to: {large_indel_file}")
    
    return small_indel_matrix, large_indel_matrix, mapped_count

def generate_indel_heatmaps(small_indel_matrix, large_indel_matrix, ngs_file, output_dir):
    """
    Generate indel distribution heatmaps
    
    Parameters:
    small_indel_matrix: Chip matrix for small indels (<3bp del)
    large_indel_matrix: Chip matrix for large indels (>=3bp del)
    ngs_file: NGS file name, used to generate output file name
    output_dir: Output directory
    """
    if small_indel_matrix is None or large_indel_matrix is None:
        print("Error: Indel matrices are empty, unable to generate heatmap")
        return
    
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract file name (without extension) for output file naming
    file_name = os.path.splitext(os.path.basename(ngs_file))[0]
    
    # Create masks for valid data points
    small_mask = np.isnan(small_indel_matrix)
    large_mask = np.isnan(large_indel_matrix)
    
    # Calculate statistics
    small_mean = np.nanmean(small_indel_matrix)
    large_mean = np.nanmean(large_indel_matrix)
    small_std = np.nanstd(small_indel_matrix)
    large_std = np.nanstd(large_indel_matrix)
    
    print(f"Small indel (<3bp del) statistics:")
    print(f"  Mean: {small_mean:.6f}")
    print(f"  Std: {small_std:.6f}")
    print(f"Large indel (>=3bp del) statistics:")
    print(f"  Mean: {large_mean:.6f}")
    print(f"  Std: {large_std:.6f}")
    
    # 1. Small indel heatmap (<3bp del)
    plt.figure(figsize=(14, 12))
    
    # Use YlOrRd colormap for small indels
    small_cmap = plt.cm.YlOrRd
    
    ax = sns.heatmap(small_indel_matrix, mask=small_mask, cmap=small_cmap, 
                    cbar_kws={'label': 'Small Indel Rate (<3bp del)'}, 
                    linewidths=0, square=True, xticklabels=False, yticklabels=False)
    
    # Add a rectangle border around the chip
    height, width = small_indel_matrix.shape
    rect = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    ax.set_title('Small Indel Distribution (<3bp del) on Chip', fontsize=14)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add statistics annotation
    plt.annotate(f'Mean: {small_mean:.6f}\nStd: {small_std:.6f}', 
                xy=(0.05, 0.95), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8),
                verticalalignment='top')
    
    # Save small indel heatmap
    small_heatmap_path = os.path.join(output_dir, f'{file_name}_small_indel_heatmap.png')
    plt.savefig(small_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Small indel heatmap saved to: {small_heatmap_path}")
    
    # 2. Large indel heatmap (>=3bp del)
    plt.figure(figsize=(14, 12))
    
    # Use Reds colormap for large indels
    large_cmap = plt.cm.Reds
    
    ax = sns.heatmap(large_indel_matrix, mask=large_mask, cmap=large_cmap, 
                    cbar_kws={'label': 'Large Indel Rate (≥3bp del)'}, 
                    linewidths=0, square=True, xticklabels=False, yticklabels=False)
    
    # Add a rectangle border around the chip
    height, width = large_indel_matrix.shape
    rect = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    ax.set_title('Large Indel Distribution (≥3bp del) on Chip', fontsize=14)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add statistics annotation
    plt.annotate(f'Mean: {large_mean:.6f}\nStd: {large_std:.6f}', 
                xy=(0.05, 0.95), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8),
                verticalalignment='top')
    
    # Save large indel heatmap
    large_heatmap_path = os.path.join(output_dir, f'{file_name}_large_indel_heatmap.png')
    plt.savefig(large_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Large indel heatmap saved to: {large_heatmap_path}")
    
    # 3. Combined indel heatmap (side by side)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(28, 12))
    
    # Small indel subplot
    sns.heatmap(small_indel_matrix, mask=small_mask, cmap=small_cmap, 
                cbar_kws={'label': 'Small Indel Rate (<3bp del)'}, 
                linewidths=0, square=True, xticklabels=False, yticklabels=False, ax=ax1)
    
    ax1.set_title('Small Indel Distribution (<3bp del)', fontsize=14)
    ax1.set_xlabel('X Coordinate', fontsize=12)
    ax1.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add border
    height, width = small_indel_matrix.shape
    rect1 = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax1.add_patch(rect1)
    
    # Large indel subplot
    sns.heatmap(large_indel_matrix, mask=large_mask, cmap=large_cmap, 
                cbar_kws={'label': 'Large Indel Rate (≥3bp del)'}, 
                linewidths=0, square=True, xticklabels=False, yticklabels=False, ax=ax2)
    
    ax2.set_title('Large Indel Distribution (≥3bp del)', fontsize=14)
    ax2.set_xlabel('X Coordinate', fontsize=12)
    ax2.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add border
    rect2 = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax2.add_patch(rect2)
    
    plt.tight_layout()
    
    # Save combined heatmap
    combined_heatmap_path = os.path.join(output_dir, f'{file_name}_combined_indel_heatmap.png')
    plt.savefig(combined_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Combined indel heatmap saved to: {combined_heatmap_path}")
    
    # 4. High indel regions analysis
    # Define high indel thresholds (mean + 2*std)
    small_high_threshold = small_mean + 2 * small_std
    large_high_threshold = large_mean + 2 * large_std
    
    # Create binary masks for high indel regions
    small_high_mask = np.logical_or(small_mask, small_indel_matrix < small_high_threshold)
    large_high_mask = np.logical_or(large_mask, large_indel_matrix < large_high_threshold)
    
    # Generate high indel regions heatmap
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(28, 12))
    
    # High small indel regions
    sns.heatmap(small_indel_matrix, mask=small_high_mask, cmap='Oranges', 
                vmin=small_high_threshold, 
                cbar_kws={'label': f'High Small Indel (>{small_high_threshold:.6f})'}, 
                linewidths=0, square=True, xticklabels=False, yticklabels=False, ax=ax1)
    
    ax1.set_title(f'High Small Indel Regions (>Mean+2σ)', fontsize=14)
    ax1.set_xlabel('X Coordinate', fontsize=12)
    ax1.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add border
    rect1 = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax1.add_patch(rect1)
    
    # High large indel regions
    sns.heatmap(large_indel_matrix, mask=large_high_mask, cmap='Reds', 
                vmin=large_high_threshold,
                cbar_kws={'label': f'High Large Indel (>{large_high_threshold:.6f})'}, 
                linewidths=0, square=True, xticklabels=False, yticklabels=False, ax=ax2)
    
    ax2.set_title(f'High Large Indel Regions (>Mean+2σ)', fontsize=14)
    ax2.set_xlabel('X Coordinate', fontsize=12)
    ax2.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add border
    rect2 = Rectangle((0, 0), width, height, linewidth=2, edgecolor='black', facecolor='none')
    ax2.add_patch(rect2)
    
    plt.tight_layout()
    
    # Save high indel regions heatmap
    high_indel_path = os.path.join(output_dir, f'{file_name}_high_indel_regions.png')
    plt.savefig(high_indel_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"High indel regions heatmap saved to: {high_indel_path}")

def analyze_indel_position_bias(small_indel_matrix, large_indel_matrix, ngs_file, output_dir):
    """
    Analyze position bias for indel data
    
    Parameters:
    small_indel_matrix: Chip matrix for small indels
    large_indel_matrix: Chip matrix for large indels
    ngs_file: NGS file name, used to generate output file name
    output_dir: Output directory
    """
    if small_indel_matrix is None or large_indel_matrix is None:
        print("Error: Indel matrices are empty, unable to analyze position bias")
        return
    
    # Extract file name (without extension) for output file naming
    file_name = os.path.splitext(os.path.basename(ngs_file))[0]
    
    # Create position data for small indels
    small_positions = []
    small_values = []
    
    for y in range(small_indel_matrix.shape[0]):
        for x in range(small_indel_matrix.shape[1]):
            if not np.isnan(small_indel_matrix[y, x]):
                small_positions.append((x, y))
                small_values.append(small_indel_matrix[y, x])
    
    # Create position data for large indels
    large_positions = []
    large_values = []
    
    for y in range(large_indel_matrix.shape[0]):
        for x in range(large_indel_matrix.shape[1]):
            if not np.isnan(large_indel_matrix[y, x]):
                large_positions.append((x, y))
                large_values.append(large_indel_matrix[y, x])
    
    if not small_positions or not large_positions:
        print("Not enough data for position bias analysis")
        return
    
    # Extract X and Y coordinates
    small_x_coords = [p[0] for p in small_positions]
    small_y_coords = [p[1] for p in small_positions]
    large_x_coords = [p[0] for p in large_positions]
    large_y_coords = [p[1] for p in large_positions]
    
    # Create combined X position analysis
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Small indel vs X position
    ax1.scatter(small_x_coords, small_values, alpha=0.5, s=10, color='orange')
    ax1.set_title('Small Indel vs X Position')
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Small Indel Rate')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Add trend line
    try:
        z = np.polyfit(small_x_coords, small_values, 1)
        p = np.poly1d(z)
        x_range = np.linspace(min(small_x_coords), max(small_x_coords), 100)
        ax1.plot(x_range, p(x_range), 'r--')
        corr = np.corrcoef(small_x_coords, small_values)[0, 1]
        ax1.annotate(f'Correlation: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"Unable to calculate small indel X trend line: {e}")
    
    # Large indel vs X position
    ax2.scatter(large_x_coords, large_values, alpha=0.5, s=10, color='red')
    ax2.set_title('Large Indel vs X Position')
    ax2.set_xlabel('X Position')
    ax2.set_ylabel('Large Indel Rate')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # Add trend line
    try:
        z = np.polyfit(large_x_coords, large_values, 1)
        p = np.poly1d(z)
        x_range = np.linspace(min(large_x_coords), max(large_x_coords), 100)
        ax2.plot(x_range, p(x_range), 'r--')
        corr = np.corrcoef(large_x_coords, large_values)[0, 1]
        ax2.annotate(f'Correlation: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"Unable to calculate large indel X trend line: {e}")
    
    # Small indel vs Y position
    ax3.scatter(small_y_coords, small_values, alpha=0.5, s=10, color='orange')
    ax3.set_title('Small Indel vs Y Position')
    ax3.set_xlabel('Y Position')
    ax3.set_ylabel('Small Indel Rate')
    ax3.grid(True, linestyle='--', alpha=0.7)
    
    # Add trend line
    try:
        z = np.polyfit(small_y_coords, small_values, 1)
        p = np.poly1d(z)
        y_range = np.linspace(min(small_y_coords), max(small_y_coords), 100)
        ax3.plot(y_range, p(y_range), 'r--')
        corr = np.corrcoef(small_y_coords, small_values)[0, 1]
        ax3.annotate(f'Correlation: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"Unable to calculate small indel Y trend line: {e}")
    
    # Large indel vs Y position
    ax4.scatter(large_y_coords, large_values, alpha=0.5, s=10, color='red')
    ax4.set_title('Large Indel vs Y Position')
    ax4.set_xlabel('Y Position')
    ax4.set_ylabel('Large Indel Rate')
    ax4.grid(True, linestyle='--', alpha=0.7)
    
    # Add trend line
    try:
        z = np.polyfit(large_y_coords, large_values, 1)
        p = np.poly1d(z)
        y_range = np.linspace(min(large_y_coords), max(large_y_coords), 100)
        ax4.plot(y_range, p(y_range), 'r--')
        corr = np.corrcoef(large_y_coords, large_values)[0, 1]
        ax4.annotate(f'Correlation: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"Unable to calculate large indel Y trend line: {e}")
    
    plt.tight_layout()
    
    # Save position analysis figure
    position_analysis_path = os.path.join(output_dir, f'{file_name}_indel_position_analysis.png')
    plt.savefig(position_analysis_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Indel position analysis saved to: {position_analysis_path}")
    
    # Perform regional analysis
    max_x = max(max(small_x_coords), max(large_x_coords))
    max_y = max(max(small_y_coords), max(large_y_coords))
    mid_x = max_x / 2
    mid_y = max_y / 2
    
    # Define regions for both small and large indels
    regions = {
        "Upper Left": lambda x, y: x < mid_x and y > mid_y,
        "Upper Right": lambda x, y: x >= mid_x and y > mid_y,
        "Lower Left": lambda x, y: x < mid_x and y <= mid_y,
        "Lower Right": lambda x, y: x >= mid_x and y <= mid_y,
        "Center": lambda x, y: mid_x/2 <= x <= mid_x*1.5 and mid_y/2 <= y <= mid_y*1.5
    }
    
    # Calculate statistics for each region
    small_region_stats = {}
    large_region_stats = {}
    
    for region_name, region_func in regions.items():
        # Small indel statistics
        small_region_values = [small_values[i] for i, (x, y) in enumerate(small_positions) if region_func(x, y)]
        if small_region_values:
            small_region_stats[region_name] = {
                "Count": len(small_region_values),
                "Mean": np.mean(small_region_values),
                "Median": np.median(small_region_values),
                "Std Dev": np.std(small_region_values),
                "Min": np.min(small_region_values),
                "Max": np.max(small_region_values)
            }
        
        # Large indel statistics
        large_region_values = [large_values[i] for i, (x, y) in enumerate(large_positions) if region_func(x, y)]
        if large_region_values:
            large_region_stats[region_name] = {
                "Count": len(large_region_values),
                "Mean": np.mean(large_region_values),
                "Median": np.median(large_region_values),
                "Std Dev": np.std(large_region_values),
                "Min": np.min(large_region_values),
                "Max": np.max(large_region_values)
            }
    
    # Save region statistics
    if small_region_stats:
        small_region_df = pd.DataFrame(small_region_stats).T
        small_region_path = os.path.join(output_dir, f'{file_name}_small_indel_region_statistics.csv')
        small_region_df.to_csv(small_region_path)
        print(f"Small indel region statistics saved to: {small_region_path}")
    
    if large_region_stats:
        large_region_df = pd.DataFrame(large_region_stats).T
        large_region_path = os.path.join(output_dir, f'{file_name}_large_indel_region_statistics.csv')
        large_region_df.to_csv(large_region_path)
        print(f"Large indel region statistics saved to: {large_region_path}")

def main():
    parser = argparse.ArgumentParser(description='Map NGS indel data to chip physical location and generate heatmaps')
    parser.add_argument('ref_file', help='Path to the reference file for synthetic chip')
    parser.add_argument('ngs_file', help='Path to the NGS sequencing result file')
    parser.add_argument('--output_dir', '-o', default='chip_indel_analysis_results', 
                        help='Output directory name')
    parser.add_argument('--chip_type', '-t', choices=['small', 'large'], default='small',
                        help='Chip type: small(540x635) or large(1080x636)')
    
    args = parser.parse_args()
    
    # Create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Create chip layout
    sequences, seq_to_pos = create_chip_layout(args.ref_file, args.chip_type)
    
    # Map indel data to chip position
    small_indel_matrix, large_indel_matrix, mapped_count = map_indel_to_chip(args.ngs_file, seq_to_pos, args.output_dir)
    
    if mapped_count > 0:
        # Generate heatmaps
        generate_indel_heatmaps(small_indel_matrix, large_indel_matrix, args.ngs_file, args.output_dir)
        
        # Analyze position bias
        analyze_indel_position_bias(small_indel_matrix, large_indel_matrix, args.ngs_file, args.output_dir)
        
        print(f"\nIndel analysis complete! Results saved to '{args.output_dir}' directory")
    else:
        print("No sequences were successfully mapped, cannot generate heatmaps")

if __name__ == "__main__":
    main() 