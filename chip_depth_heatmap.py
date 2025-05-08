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
from matplotlib.colors import LogNorm

def create_chip_layout(ref_file, chip_type='small'):
    """
    读取参考文件并创建芯片布局矩阵
    
    参数:
    ref_file: 合成芯片的参考文件
    chip_type: 芯片类型 'small'(635x540) 或 'large'(636x1080)
    
    返回:
    sequences: 所有合成序列列表
    seq_to_pos: 序列到芯片位置的映射
    """
    print(f"正在读取参考文件: {ref_file}")
    
    # 设置芯片尺寸
    if chip_type == 'small':
        width, height = 635, 540
        expected_lines = 342583
    else:  # large
        width, height = 636, 1080
        expected_lines = 686880
    
    # 读取参考文件
    with open(ref_file, 'r') as f:
        ref_sequences = [line.strip() for line in f]
    
    print(f"读取了 {len(ref_sequences)} 行序列")
    
    # 确保序列数量符合预期
    if len(ref_sequences) < expected_lines:
        print(f"警告: 文件中序列数量少于预期值 ({len(ref_sequences)} < {expected_lines})")
        # 使用空序列填充到预期行数
        ref_sequences.extend(["0" * 150] * (expected_lines - len(ref_sequences)))
    
    # 提取每个序列的16-120位子串
    sequences = [seq[15:120] if len(seq) >= 120 else seq for seq in ref_sequences]
    
    # 创建序列到芯片位置的映射
    seq_to_pos = {}
    
    index = 0
    for x in range(width):
        for y in range(height):
            if index < len(sequences):
                if sequences[index] != "0" * 105:  # 排除填充序列
                    # 注意：排列顺序是从下到上，从左到右
                    seq_to_pos[sequences[index]] = (x, y)
                index += 1
    
    print(f"创建了 {len(seq_to_pos)} 个序列到位置的映射")
    return sequences, seq_to_pos

def map_depth_to_chip(ngs_file, seq_to_pos, output_dir):
    """
    将NGS深度数据映射到芯片位置并生成热图
    
    参数:
    ngs_file: NGS测序结果文件
    seq_to_pos: 序列到芯片位置的映射
    output_dir: 输出目录
    
    返回:
    chip_matrix: 芯片深度矩阵
    mapped_count: 成功映射的序列数量
    """
    print(f"正在读取NGS文件: {ngs_file}")
    
    # 读取NGS测序结果
    try:
        # 首先尝试以TSV格式读取
        df = pd.read_csv(ngs_file, sep='\t')
        print("成功以TSV格式读取文件")
    except Exception as e1:
        try:
            # 再尝试作为Excel文件读取
            if ngs_file.endswith('.xls') or ngs_file.endswith('.xlsx'):
                df = pd.read_excel(ngs_file)
            else:
                raise Exception("无法识别文件格式")
        except Exception as e2:
            print(f"读取文件出错: {e1}\n第二次尝试错误: {e2}")
            return None, 0
    
    # 检查必需列
    if 'Oligo seq' not in df.columns or 'Depth' not in df.columns:
        print(f"错误: 未找到必需的列。可用列: {', '.join(df.columns)}")
        return None, 0
    
    # 获取芯片尺寸
    max_x = max([pos[0] for pos in seq_to_pos.values()]) + 1
    max_y = max([pos[1] for pos in seq_to_pos.values()]) + 1
    
    # 创建芯片深度矩阵，初始值为NaN（表示没有数据）
    chip_matrix = np.full((max_y, max_x), np.nan)
    
    # 映射深度到芯片位置
    mapped_count = 0
    not_found = 0
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="映射深度数据"):
        oligo_seq = row['Oligo seq']
        depth = row['Depth']
        
        if oligo_seq in seq_to_pos:
            x, y = seq_to_pos[oligo_seq]
            # 注意：矩阵索引是(行,列)，对应(y,x)
            chip_matrix[y, x] = depth
            mapped_count += 1
        else:
            not_found += 1
    
    mapping_rate = mapped_count / len(df) * 100 if len(df) > 0 else 0
    print(f"总序列数: {len(df)}")
    print(f"成功映射: {mapped_count} ({mapping_rate:.2f}%)")
    print(f"未找到映射: {not_found}")
    
    # 保存芯片深度矩阵到CSV
    df_matrix = pd.DataFrame(chip_matrix)
    matrix_file = os.path.join(output_dir, 'chip_depth_matrix.csv')
    df_matrix.to_csv(matrix_file)
    print(f"芯片深度矩阵已保存到: {matrix_file}")
    
    return chip_matrix, mapped_count

def generate_heatmaps(chip_matrix, ngs_file, output_dir):
    """
    生成深度分布热图
    
    参数:
    chip_matrix: 芯片深度矩阵
    ngs_file: NGS文件名，用于生成输出文件名
    output_dir: 输出目录
    """
    if chip_matrix is None:
        print("错误: 芯片深度矩阵为空，无法生成热图")
        return
    
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 提取文件名（不含扩展名）用于输出文件命名
    file_name = os.path.splitext(os.path.basename(ngs_file))[0]
    
    # 为有效数据点创建掩码
    mask = np.isnan(chip_matrix)
    
    # 1. 标准线性热图
    plt.figure(figsize=(12, 10))
    ax = sns.heatmap(chip_matrix, mask=mask, cmap='viridis', 
                    cbar_kws={'label': 'Depth'})
    ax.set_title('芯片上的测序深度分布')
    ax.set_xlabel('X 坐标')
    ax.set_ylabel('Y 坐标')
    
    # 保存标准热图
    linear_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_linear.png')
    plt.savefig(linear_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"线性热图已保存到: {linear_heatmap_path}")
    
    # 2. 对数尺度热图（便于观察极端值）
    plt.figure(figsize=(12, 10))
    # 为对数尺度创建掩码 - 排除负值和零值
    log_mask = np.logical_or(mask, chip_matrix <= 0)
    
    # 对有效值取对数
    with np.errstate(invalid='ignore'):
        log_matrix = np.log10(chip_matrix)
    
    ax = sns.heatmap(log_matrix, mask=log_mask, cmap='viridis',
                    cbar_kws={'label': 'log10(Depth)'})
    ax.set_title('芯片上的测序深度分布 (对数尺度)')
    ax.set_xlabel('X 坐标')
    ax.set_ylabel('Y 坐标')
    
    # 保存对数热图
    log_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_log.png')
    plt.savefig(log_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"对数热图已保存到: {log_heatmap_path}")
    
    # 3. 分位数归一化热图 - 突出相对高低
    plt.figure(figsize=(12, 10))
    
    # 提取有效数据点
    valid_data = chip_matrix[~mask]
    if len(valid_data) > 0:
        # 计算分位数
        q_low = np.nanpercentile(chip_matrix, 5)
        q_high = np.nanpercentile(chip_matrix, 95)
        
        # 创建自定义色标，突出显示异常值
        ax = sns.heatmap(chip_matrix, mask=mask, cmap='coolwarm', 
                        vmin=q_low, vmax=q_high,
                        cbar_kws={'label': 'Depth (5%-95% 范围)'})
        ax.set_title('芯片上的测序深度分布 (分位数归一化)')
        ax.set_xlabel('X 坐标')
        ax.set_ylabel('Y 坐标')
        
        # 保存分位数热图
        quantile_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_quantile.png')
        plt.savefig(quantile_heatmap_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"分位数热图已保存到: {quantile_heatmap_path}")
        
    # 4. 创建Z-score热图，更好地表示偏差
    plt.figure(figsize=(12, 10))
    # 计算Z分数
    mean = np.nanmean(chip_matrix)
    std = np.nanstd(chip_matrix)
    if std > 0:
        z_matrix = (chip_matrix - mean) / std
        
        # 限制Z分数范围，便于可视化
        z_matrix = np.clip(z_matrix, -3, 3)
        
        ax = sns.heatmap(z_matrix, mask=mask, cmap='coolwarm', 
                        vmin=-3, vmax=3,
                        cbar_kws={'label': 'Z-score'})
        ax.set_title('芯片上的测序深度分布 (Z-score标准化)')
        ax.set_xlabel('X 坐标')
        ax.set_ylabel('Y 坐标')
        
        # 保存Z-score热图
        zscore_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_zscore.png')
        plt.savefig(zscore_heatmap_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Z-score热图已保存到: {zscore_heatmap_path}")
    
    # 5. 生成相对深度分布热图
    # 按X轴（列）和Y轴（行）计算均值
    row_means = np.nanmean(chip_matrix, axis=1)
    col_means = np.nanmean(chip_matrix, axis=0)
    
    # 创建相对热图
    relative_matrix = np.zeros_like(chip_matrix)
    for i in range(chip_matrix.shape[0]):
        for j in range(chip_matrix.shape[1]):
            if not np.isnan(chip_matrix[i, j]):
                # 计算相对值: 实际值/位置均值
                row_mean = row_means[i] if not np.isnan(row_means[i]) else mean
                col_mean = col_means[j] if not np.isnan(col_means[j]) else mean
                position_mean = (row_mean + col_mean) / 2
                if position_mean > 0:
                    relative_matrix[i, j] = chip_matrix[i, j] / position_mean
                else:
                    relative_matrix[i, j] = np.nan
            else:
                relative_matrix[i, j] = np.nan
    
    plt.figure(figsize=(12, 10))
    rel_mask = np.isnan(relative_matrix)
    ax = sns.heatmap(relative_matrix, mask=rel_mask, cmap='coolwarm', 
                    vmin=0.5, vmax=1.5, center=1.0,
                    cbar_kws={'label': 'Relative Depth'})
    ax.set_title('芯片上的相对测序深度分布')
    ax.set_xlabel('X 坐标')
    ax.set_ylabel('Y 坐标')
    
    # 保存相对深度热图
    relative_heatmap_path = os.path.join(output_dir, f'{file_name}_depth_heatmap_relative.png')
    plt.savefig(relative_heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"相对深度热图已保存到: {relative_heatmap_path}")

def analyze_position_bias(chip_matrix, ngs_file, output_dir):
    """
    分析位置偏差，查看深度与芯片位置的关系
    
    参数:
    chip_matrix: 芯片深度矩阵
    ngs_file: NGS文件名，用于生成输出文件名
    output_dir: 输出目录
    """
    if chip_matrix is None:
        print("错误: 芯片深度矩阵为空，无法分析位置偏差")
        return
    
    # 提取文件名（不含扩展名）用于输出文件命名
    file_name = os.path.splitext(os.path.basename(ngs_file))[0]
    
    # 创建位置数据
    positions = []
    depths = []
    
    for y in range(chip_matrix.shape[0]):
        for x in range(chip_matrix.shape[1]):
            if not np.isnan(chip_matrix[y, x]):
                positions.append((x, y))
                depths.append(chip_matrix[y, x])
    
    if not positions:
        print("没有足够的数据进行位置偏差分析")
        return
    
    # 提取X和Y坐标
    x_coords = [p[0] for p in positions]
    y_coords = [p[1] for p in positions]
    
    # 创建基于X位置的分析
    plt.figure(figsize=(12, 6))
    plt.scatter(x_coords, depths, alpha=0.5, s=10)
    plt.title('深度与X位置关系')
    plt.xlabel('X位置')
    plt.ylabel('深度')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # 添加趋势线
    try:
        z = np.polyfit(x_coords, depths, 1)
        p = np.poly1d(z)
        x_range = np.linspace(min(x_coords), max(x_coords), 100)
        plt.plot(x_range, p(x_range), 'r--')
        corr = np.corrcoef(x_coords, depths)[0, 1]
        plt.annotate(f'相关系数: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"无法计算X趋势线: {e}")
    
    # 保存X位置分析图
    x_analysis_path = os.path.join(output_dir, f'{file_name}_x_position_analysis.png')
    plt.savefig(x_analysis_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 创建基于Y位置的分析
    plt.figure(figsize=(12, 6))
    plt.scatter(y_coords, depths, alpha=0.5, s=10)
    plt.title('深度与Y位置关系')
    plt.xlabel('Y位置')
    plt.ylabel('深度')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # 添加趋势线
    try:
        z = np.polyfit(y_coords, depths, 1)
        p = np.poly1d(z)
        y_range = np.linspace(min(y_coords), max(y_coords), 100)
        plt.plot(y_range, p(y_range), 'r--')
        corr = np.corrcoef(y_coords, depths)[0, 1]
        plt.annotate(f'相关系数: {corr:.4f}', 
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    except Exception as e:
        print(f"无法计算Y趋势线: {e}")
    
    # 保存Y位置分析图
    y_analysis_path = os.path.join(output_dir, f'{file_name}_y_position_analysis.png')
    plt.savefig(y_analysis_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 进行四分位和中心区域的统计分析
    max_x = max(x_coords)
    max_y = max(y_coords)
    mid_x = max_x / 2
    mid_y = max_y / 2
    
    # 定义区域
    regions = {
        "左上": [(x, y) for x, y in positions if x < mid_x and y > mid_y],
        "右上": [(x, y) for x, y in positions if x >= mid_x and y > mid_y],
        "左下": [(x, y) for x, y in positions if x < mid_x and y <= mid_y],
        "右下": [(x, y) for x, y in positions if x >= mid_x and y <= mid_y],
        "中心": [(x, y) for x, y in positions if mid_x/2 <= x <= mid_x*1.5 and mid_y/2 <= y <= mid_y*1.5]
    }
    
    # 计算每个区域的统计数据
    region_stats = {}
    for region, pos in regions.items():
        if pos:
            region_depths = [depths[positions.index(p)] for p in pos if p in positions]
            if region_depths:
                region_stats[region] = {
                    "数量": len(region_depths),
                    "平均深度": np.mean(region_depths),
                    "中位数深度": np.median(region_depths),
                    "标准差": np.std(region_depths),
                    "最小值": np.min(region_depths),
                    "最大值": np.max(region_depths)
                }
    
    # 保存区域统计数据
    region_stats_df = pd.DataFrame(region_stats).T
    region_stats_path = os.path.join(output_dir, f'{file_name}_region_statistics.csv')
    region_stats_df.to_csv(region_stats_path)
    print(f"区域统计数据已保存到: {region_stats_path}")
    
    # 创建区域统计条形图
    if region_stats:
        regions = list(region_stats.keys())
        means = [region_stats[r]["平均深度"] for r in regions]
        
        plt.figure(figsize=(10, 6))
        bars = plt.bar(regions, means)
        
        # 添加数据标签
        for bar, mean in zip(bars, means):
            plt.text(bar.get_x() + bar.get_width()/2, mean + 5, 
                    f'{mean:.1f}', 
                    ha='center', va='bottom')
        
        plt.title('不同区域的平均测序深度')
        plt.ylabel('平均深度')
        plt.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # 保存区域统计图
        region_plot_path = os.path.join(output_dir, f'{file_name}_region_depth_analysis.png')
        plt.savefig(region_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"区域统计图已保存到: {region_plot_path}")

def main():
    parser = argparse.ArgumentParser(description='将NGS深度数据映射到芯片物理位置并生成热图')
    parser.add_argument('ref_file', help='合成芯片的参考文件路径')
    parser.add_argument('ngs_file', help='NGS测序结果文件路径')
    parser.add_argument('--output_dir', '-o', default='chip_analysis_results', 
                        help='输出目录名称')
    parser.add_argument('--chip_type', '-t', choices=['small', 'large'], default='small',
                        help='芯片类型: small(635x540) 或 large(636x1080)')
    
    args = parser.parse_args()
    
    # 创建输出目录
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # 创建芯片布局
    sequences, seq_to_pos = create_chip_layout(args.ref_file, args.chip_type)
    
    # 映射深度到芯片位置
    chip_matrix, mapped_count = map_depth_to_chip(args.ngs_file, seq_to_pos, args.output_dir)
    
    if mapped_count > 0:
        # 生成热图
        generate_heatmaps(chip_matrix, args.ngs_file, args.output_dir)
        
        # 分析位置偏差
        analyze_position_bias(chip_matrix, args.ngs_file, args.output_dir)
        
        print(f"\n分析完成! 结果已保存到 '{args.output_dir}' 目录")
    else:
        print("没有成功映射任何序列，无法生成热图")

if __name__ == "__main__":
    main() 