#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import re
from datetime import datetime
import base64
from io import BytesIO
from matplotlib.ticker import PercentFormatter

# 设置绘图风格
plt.style.use('ggplot')
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def count_nucleotides(seq):
    """统计序列中碱基的数量"""
    seq = seq.upper()
    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')
    
    length = len(seq)
    gc_content = (g_count + c_count) / length if length > 0 else 0
    ag_content = (a_count + g_count) / length if length > 0 else 0
    
    return {
        'A': a_count / length if length > 0 else 0,
        'T': t_count / length if length > 0 else 0,
        'G': g_count / length if length > 0 else 0,
        'C': c_count / length if length > 0 else 0,
        'GC': gc_content,
        'AG': ag_content
    }

def fig_to_base64(fig):
    """将matplotlib图形转换为base64编码，用于HTML报告"""
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode('utf-8')
    return img_str

def analyze_oligo_ngs(file_path, output_dir=None):
    """
    分析NGS oligopool测序数据
    
    参数:
    file_path: 输入文件路径
    output_dir: 输出结果根目录
    """
    # 提取文件名作为输出子目录名
    file_name = os.path.basename(file_path)
    file_name_no_ext = os.path.splitext(file_name)[0]
    
    # 如果未指定output_dir，则在当前目录下创建results目录
    if output_dir is None:
        output_dir = 'results'
    
    # 为当前分析数据创建子目录
    analysis_dir = os.path.join(output_dir, file_name_no_ext)
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    # 在子目录中创建images和data目录
    img_dir = os.path.join(analysis_dir, 'images')
    data_dir = os.path.join(analysis_dir, 'data')
    
    for dir_path in [img_dir, data_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    
    # 读取数据文件
    print(f"正在读取文件: {file_path}...")
    try:
        # 首先尝试以TSV格式读取文件(无论扩展名)
        df = pd.read_csv(file_path, sep='\t')
        print(f"成功以TSV格式读取文件")
    except Exception as e1:
        print(f"以TSV格式读取失败，尝试其他格式...")
        try:
            # 再尝试作为Excel文件读取
            if file_path.endswith('.xls') or file_path.endswith('.xlsx'):
                df = pd.read_excel(file_path)
            elif file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                print(f"无法确定文件格式，尝试作为通用文本文件读取...")
                with open(file_path, 'r') as f:
                    first_line = f.readline().strip()
                    if '\t' in first_line:
                        df = pd.read_csv(file_path, sep='\t')
                    elif ',' in first_line:
                        df = pd.read_csv(file_path)
                    else:
                        raise Exception("无法自动检测文件分隔符")
        except Exception as e2:
            print(f"读取文件出错: {e1}\n第二次尝试错误: {e2}")
            return
    
    # 检查必需列
    required_columns = ['Oligo seq', 'Depth']
    for col in required_columns:
        if col not in df.columns:
            print(f"错误: 未找到'{col}'列。")
            print(f"可用列: {', '.join(df.columns)}")
            return
    
    # 保存原始数据副本用于分析
    oligo_data = df.copy()
    
    # 分析步骤1: 计算碱基组成
    print("正在分析序列碱基组成...")
    nucleotide_data = []
    
    for _, row in oligo_data.iterrows():
        seq = row['Oligo seq']
        nucleotide_counts = count_nucleotides(seq)
        nucleotide_data.append({
            'Oligo name': row.get('Oligo name', f"Oligo_{_}"),
            'Depth': row['Depth'],
            'A_content': nucleotide_counts['A'],
            'T_content': nucleotide_counts['T'],
            'G_content': nucleotide_counts['G'],
            'C_content': nucleotide_counts['C'],
            'GC_content': nucleotide_counts['GC'],
            'AG_content': nucleotide_counts['AG']
        })
    
    # 创建包含碱基含量的DataFrame
    composition_df = pd.DataFrame(nucleotide_data)
    
    # 将碱基组成数据保存到CSV
    composition_csv = os.path.join(data_dir, 'nucleotide_composition.csv')
    composition_df.to_csv(composition_csv, index=False)
    print(f"已保存碱基组成数据: {composition_csv}")
    
    # 分析步骤2: Depth 分布分析
    print("正在分析深度分布...")
    depth_data = oligo_data['Depth']
    
    # 基本统计信息
    depth_stats = depth_data.describe()
    print("\n深度统计:")
    print(depth_stats)
    
    # 保存统计数据
    stats_csv = os.path.join(data_dir, 'depth_statistics.csv')
    depth_stats.to_csv(stats_csv)
    
    # 图片列表，用于生成报告
    report_images = {}
    
    # 2.1 深度分布直方图与CDF曲线并排
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 直方图
    bins = min(50, int(np.sqrt(len(depth_data))))
    sns.histplot(depth_data, bins=bins, kde=True, ax=ax1)
    ax1.axvline(depth_stats['mean'], color='r', linestyle='--', 
               label=f'平均值: {depth_stats["mean"]:.2f}')
    ax1.axvline(depth_stats['50%'], color='g', linestyle='--', 
               label=f'中位数: {depth_stats["50%"]:.2f}')
    ax1.set_title('测序深度分布直方图')
    ax1.set_xlabel('深度')
    ax1.set_ylabel('频率')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # 计算CDF (累积分布函数)
    sorted_depths = np.sort(depth_data.values)
    cumulative_prob = np.arange(1, len(sorted_depths) + 1) / len(sorted_depths)
    
    # 保存CDF数据
    cdf_df = pd.DataFrame({
        'Depth': sorted_depths,
        'Cumulative_Probability': cumulative_prob
    })
    cdf_csv = os.path.join(data_dir, 'depth_cdf.csv')
    cdf_df.to_csv(cdf_csv, index=False)
    
    # 计算比值 0.95/0.05
    depth_95 = np.percentile(sorted_depths, 95)
    depth_5 = np.percentile(sorted_depths, 5)
    ratio_95_5 = depth_95 / depth_5 if depth_5 > 0 else float('inf')
    
    print(f"深度分析结果:")
    print(f"  95%分位数深度: {depth_95:.2f}")
    print(f"  5%分位数深度: {depth_5:.2f}")
    print(f"  95%/5%深度比值: {ratio_95_5:.2f}")
    
    # CDF曲线
    ax2.plot(sorted_depths, cumulative_prob, '-', linewidth=2)
    ax2.axhline(0.95, color='r', linestyle='--', 
               label=f'95%分位: {depth_95:.2f}')
    ax2.axhline(0.05, color='g', linestyle='--', 
               label=f'5%分位: {depth_5:.2f}')
    ax2.set_title(f'深度累积分布函数 (CDF) - 95%/5%比值: {ratio_95_5:.2f}')
    ax2.set_xlabel('深度')
    ax2.set_ylabel('累积概率')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()
    
    plt.tight_layout()
    
    # 保存图像
    depth_dist_path = os.path.join(img_dir, 'depth_distribution.png')
    plt.savefig(depth_dist_path, dpi=300, bbox_inches='tight')
    report_images['depth_distribution'] = fig_to_base64(plt.gcf())
    plt.close()
    
    # 3. 分析深度与碱基含量关系
    print("正在分析深度与碱基含量的关系...")
    
    # 3.1 四种碱基含量与深度关系图拼在一起
    base_colors = {'A': 'green', 'T': 'red', 'G': 'purple', 'C': 'blue'}
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for i, base in enumerate(['A', 'T', 'G', 'C']):
        content_col = f'{base}_content'
        axes[i].scatter(composition_df[content_col], composition_df['Depth'], 
                      alpha=0.5, s=10, color=base_colors[base])
        
        # 添加拟合线
        z = np.polyfit(composition_df[content_col], composition_df['Depth'], 1)
        p = np.poly1d(z)
        axes[i].plot(np.unique(composition_df[content_col]), 
                   p(np.unique(composition_df[content_col])), 
                   color='black', linestyle='--')
        
        axes[i].set_title(f'深度与{base}含量的关系')
        axes[i].set_xlabel(f'{base}含量')
        axes[i].set_ylabel('深度')
        axes[i].grid(True, linestyle='--', alpha=0.7)
        
        # 计算相关系数
        corr = np.corrcoef(composition_df[content_col], composition_df['Depth'])[0, 1]
        axes[i].annotate(f'相关系数: {corr:.4f}', 
                       xy=(0.05, 0.95), xycoords='axes fraction',
                       bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    
    plt.tight_layout()
    
    # 保存单碱基关系图
    bases_path = os.path.join(img_dir, 'depth_vs_bases.png')
    plt.savefig(bases_path, dpi=300, bbox_inches='tight')
    report_images['depth_vs_bases'] = fig_to_base64(plt.gcf())
    plt.close()
    
    # 3.2 GC和AG含量与深度关系图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # GC含量关系
    ax1.scatter(composition_df['GC_content'], composition_df['Depth'], 
              alpha=0.5, s=10, color='teal')
    
    # 添加拟合线
    z = np.polyfit(composition_df['GC_content'], composition_df['Depth'], 1)
    p = np.poly1d(z)
    ax1.plot(np.unique(composition_df['GC_content']), 
           p(np.unique(composition_df['GC_content'])), 
           color='black', linestyle='--')
    
    ax1.set_title('深度与GC含量的关系')
    ax1.set_xlabel('GC含量')
    ax1.set_ylabel('深度')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # 计算相关系数
    corr = np.corrcoef(composition_df['GC_content'], composition_df['Depth'])[0, 1]
    ax1.annotate(f'相关系数: {corr:.4f}', 
               xy=(0.05, 0.95), xycoords='axes fraction',
               bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    
    # AG含量关系
    ax2.scatter(composition_df['AG_content'], composition_df['Depth'], 
              alpha=0.5, s=10, color='orange')
    
    # 添加拟合线
    z = np.polyfit(composition_df['AG_content'], composition_df['Depth'], 1)
    p = np.poly1d(z)
    ax2.plot(np.unique(composition_df['AG_content']), 
           p(np.unique(composition_df['AG_content'])), 
           color='black', linestyle='--')
    
    ax2.set_title('深度与AG含量的关系')
    ax2.set_xlabel('AG含量')
    ax2.set_ylabel('深度')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # 计算相关系数
    corr = np.corrcoef(composition_df['AG_content'], composition_df['Depth'])[0, 1]
    ax2.annotate(f'相关系数: {corr:.4f}', 
               xy=(0.05, 0.95), xycoords='axes fraction',
               bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))
    
    plt.tight_layout()
    
    # 保存GC/AG关系图
    gc_ag_path = os.path.join(img_dir, 'depth_vs_gc_ag.png')
    plt.savefig(gc_ag_path, dpi=300, bbox_inches='tight')
    report_images['depth_vs_gc_ag'] = fig_to_base64(plt.gcf())
    plt.close()
    
    # 计算相关性矩阵(用于分析，但不再展示)
    correlation_data = pd.DataFrame({
        'Depth': composition_df['Depth'],
        'A_content': composition_df['A_content'],
        'T_content': composition_df['T_content'],
        'G_content': composition_df['G_content'],
        'C_content': composition_df['C_content'],
        'GC_content': composition_df['GC_content'],
        'AG_content': composition_df['AG_content']
    })
    
    corr_matrix = correlation_data.corr()
    
    # 保存到CSV，但不生成图表
    corr_csv = os.path.join(data_dir, 'correlation_matrix.csv')
    corr_matrix.to_csv(corr_csv)
    
    # 分析步骤4: 生成HTML报告
    print("正在生成分析报告...")
    
    # 创建HTML报告
    html_report = generate_html_report(
        file_name=file_name,
        stats=depth_stats,
        depth_95=depth_95,
        depth_5=depth_5,
        ratio_95_5=ratio_95_5,
        corr_matrix=corr_matrix,
        images=report_images
    )
    
    # 保存HTML报告
    report_path = os.path.join(analysis_dir, f'{file_name_no_ext}_analysis_report.html')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_report)
    
    print(f"\n分析完成！结果已保存到 '{analysis_dir}' 目录")
    print(f"分析报告: {report_path}")
    
    return analysis_dir

def generate_html_report(file_name, stats, depth_95, depth_5, ratio_95_5, corr_matrix, images):
    """生成分析报告的HTML"""
    
    # 基本统计摘要表格
    stats_table = '<table class="table table-striped table-bordered">\n'
    stats_table += '<thead><tr><th>统计量</th><th>值</th></tr></thead>\n<tbody>\n'
    
    for idx in stats.index:
        stats_table += f'<tr><td>{idx}</td><td>{stats[idx]:.4f}</td></tr>\n'
    
    # 添加分位数比值
    stats_table += f'<tr><td>95%分位数</td><td>{depth_95:.4f}</td></tr>\n'
    stats_table += f'<tr><td>5%分位数</td><td>{depth_5:.4f}</td></tr>\n'
    stats_table += f'<tr><td>95%/5%比值</td><td>{ratio_95_5:.4f}</td></tr>\n'
    stats_table += '</tbody></table>'
    
    # 获取当前时间
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 根据比值评估质量等级
    if ratio_95_5 < 3.0:
        quality_class = "text-success"
        quality_text = "(优)"
        quality_description = "表明序列之间的覆盖度非常均匀。"
    elif ratio_95_5 < 4.5:
        quality_class = "text-success"
        quality_text = "(良)"
        quality_description = "表明序列之间的覆盖度较为均匀。"
    elif ratio_95_5 < 6.5:
        quality_class = "text-warning"
        quality_text = "(一般)"
        quality_description = "表明不同序列间覆盖度存在一定差异。"
    else:
        quality_class = "text-danger"
        quality_text = "(较差)"
        quality_description = "表明不同序列间覆盖度差异较大，可能存在偏好性扩增或测序偏好性。"
    
    # 创建HTML报告
    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NGS Oligopool分析报告</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 20px; }}
        h1, h2, h3 {{ color: #2c3e50; }}
        .img-fluid {{ max-width: 100%; height: auto; }}
        .chart-container {{ margin-bottom: 30px; }}
        .footer {{ margin-top: 30px; text-align: center; font-size: 0.8em; color: #7f8c8d; }}
    </style>
</head>
<body>
    <div class="container">
        <header class="my-4">
            <h1 class="text-center">NGS Oligopool测序分析报告</h1>
            <p class="text-center text-muted">文件: {file_name} | 生成时间: {now}</p>
            <hr>
        </header>

        <section id="summary">
            <h2>分析摘要</h2>
            <div class="row">
                <div class="col-md-6">
                    <h3>基本统计数据</h3>
                    {stats_table}
                </div>
                <div class="col-md-6">
                    <h3>分析结论</h3>
                    <ul class="list-group">
                        <li class="list-group-item">
                            <strong>覆盖均匀性:</strong> 
                            <span class="{quality_class}">
                                95%/5%深度比值为 {ratio_95_5:.2f} {quality_text}
                            </span>
                        </li>
                        <li class="list-group-item">
                            <strong>平均测序深度:</strong> {stats['mean']:.2f}
                        </li>
                        <li class="list-group-item">
                            <strong>深度分布范围:</strong> {stats['min']:.0f} - {stats['max']:.0f}
                        </li>
                    </ul>
                </div>
            </div>
        </section>

        <section id="depth-distribution" class="my-5">
            <h2>深度分布分析</h2>
            
            <div class="chart-container">
                <h3>深度分布直方图与累积分布函数</h3>
                <div class="text-center">
                    <img src="data:image/png;base64,{images['depth_distribution']}" class="img-fluid" alt="深度分布与CDF">
                </div>
                <div class="my-3">
                    <p class="text-center">
                        <strong>95%/5%深度比值:</strong> {ratio_95_5:.2f} 
                        <span class="{quality_class}">
                            {quality_text}
                        </span>
                    </p>
                </div>
            </div>
        </section>

        <section id="base-composition" class="my-5">
            <h2>碱基组成与深度关系分析</h2>
            
            <div class="chart-container">
                <h3>深度与ATCG各碱基含量关系</h3>
                <div class="text-center">
                    <img src="data:image/png;base64,{images['depth_vs_bases']}" class="img-fluid" alt="深度与ATCG含量关系">
                </div>
            </div>
            
            <div class="chart-container">
                <h3>深度与GC含量和AG含量关系</h3>
                <div class="text-center">
                    <img src="data:image/png;base64,{images['depth_vs_gc_ag']}" class="img-fluid" alt="深度与GC/AG含量关系">
                </div>
            </div>
        </section>

        <section id="conclusions" class="my-5">
            <h2>结论与建议</h2>
            <div class="card">
                <div class="card-body">
                    <h5 class="card-title">数据质量评估</h5>
                    <p class="card-text">
                        基于分析结果，我们可以得出以下结论：
                    </p>
                    <ul>
                        <li>测序深度中位数为 {stats['50%']:.2f}，平均值为 {stats['mean']:.2f}。</li>
                        <li>95%/5%深度比值为 {ratio_95_5:.2f}，{quality_description}</li>
                        <li>
                            {'深度与GC含量呈现明显相关性，建议在下一步分析中考虑GC偏好性的影响。' if abs(corr_matrix.loc['Depth', 'GC_content']) > 0.2 else '深度与GC含量相关性不明显，说明测序过程中GC偏好性较小。'}
                        </li>
                    </ul>
                </div>
            </div>
        </section>
        
        <footer class="footer">
            <p>NGS Oligopool QC 分析报告 | 生成时间: {now}</p>
        </footer>
    </div>
</body>
</html>
"""
    return html

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python analyze_oligo_ngs.py <输入文件> [输出目录]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    analyze_oligo_ngs(input_file, output_dir) 