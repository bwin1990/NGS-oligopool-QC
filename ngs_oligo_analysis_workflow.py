#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import subprocess
import webbrowser
import time
import json
import tkinter as tk
from tkinter import filedialog
from datetime import datetime
from pathlib import Path
import tempfile
import threading
import uuid

# 尝试导入分析模块
try:
    import analyze_oligo_ngs
    import chip_depth_heatmap
except ImportError:
    print("警告: 无法导入分析模块。请确保 analyze_oligo_ngs.py 和 chip_depth_heatmap.py 文件在当前目录中。")
    print("程序将继续运行，但某些功能可能不可用。")

class NGSOligoAnalysisWorkflow:
    """NGS Oligopool分析工作流"""
    
    def __init__(self):
        """初始化分析工作流"""
        self.input_dir = None
        self.output_dir = None
        self.config_data = None
        self.analysis_results = []
    
    def scan_directory(self, directory_path):
        """扫描目录获取可分析文件，生成简单的配置文本文件"""
        if not os.path.exists(directory_path):
            print(f"错误: 目录不存在 - {directory_path}")
            return False
        
        # 设置输入目录
        self.input_dir = directory_path
        self.output_dir = os.path.join(directory_path, "results")
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 查找目录中的所有90.oligostat.xls文件
        ngs_files = []
        for ext in ['*90.oligostat.xls', '*90.oligostat.xlsx']:
            ngs_files.extend(glob.glob(os.path.join(directory_path, ext)))
        
        # 按修改时间排序
        ngs_files.sort(key=os.path.getmtime, reverse=True)
        
        if not ngs_files:
            print(f"在目录 {directory_path} 中未找到任何 *90.oligostat.xls 文件")
            return False
        
        print(f"找到 {len(ngs_files)} 个测序结果文件")
        
        # 生成配置文本文件
        config_txt_path = os.path.join(self.output_dir, "analysis_config.txt")
        
        with open(config_txt_path, 'w', encoding='utf-8') as f:
            f.write("# NGS Oligopool测序分析配置文件\n")
            f.write("# 格式: NGS测序文件的绝对路径 [Tab] 芯片打印序列文件的绝对路径 [Tab] 芯片类型(small或large)\n")
            f.write("# 请在每行的NGS文件路径后添加Tab键和对应的芯片打印序列文件路径，再添加Tab键和芯片类型\n")
            f.write("# 芯片类型: small(540x635)或large(1080x636)\n")
            f.write("# 保存文件后，程序将自动开始分析\n\n")
            
            for ngs_file in ngs_files:
                # 默认值: [NGS文件路径] [Tab] [留空供用户填写] [Tab] small
                f.write(f"{os.path.abspath(ngs_file)}\t\tsmall\n")
        
        print(f"配置文件已创建: {config_txt_path}")
        
        # 尝试打开文本文件供用户编辑
        try:
            if sys.platform == "win32":
                os.startfile(config_txt_path)
            elif sys.platform == "darwin":  # macOS
                subprocess.run(["open", config_txt_path])
            else:  # Linux
                subprocess.run(["xdg-open", config_txt_path])
            
            print("\n请在打开的文本文件中：")
            print("1. 在每行的NGS文件路径后添加Tab键和对应的芯片打印序列文件绝对路径")
            print("2. 根据需要修改芯片类型(small或large)")
            print("3. 保存文件")
            print("\n配置完成后，请运行以下命令开始分析:")
            print(f"python {os.path.basename(__file__)} --config {config_txt_path}")
            
            return True
            
        except Exception as e:
            print(f"无法自动打开配置文件: {e}")
            print(f"请手动打开并编辑配置文件: {config_txt_path}")
            return True
    
    def load_config(self, config_file):
        """从配置文本文件加载分析设置"""
        if not os.path.exists(config_file):
            print(f"错误: 配置文件不存在 - {config_file}")
            return False
        
        # 判断是JSON还是TXT格式
        config_ext = os.path.splitext(config_file)[1].lower()
        
        try:
            if config_ext == ".json":
                # 处理JSON格式配置文件
                with open(config_file, 'r', encoding='utf-8') as f:
                    self.config_data = json.load(f)
                
                # 验证基本配置
                if not self.config_data.get('input_dir'):
                    print("错误: 配置文件中缺少input_dir字段")
                    return False
                
                if not self.config_data.get('files'):
                    print("错误: 配置文件中缺少files字段或文件列表为空")
                    return False
                
                # 设置目录
                self.input_dir = self.config_data.get('input_dir')
                self.output_dir = self.config_data.get('output_dir', os.path.join(self.input_dir, "results"))
                
                # 准备文件配置列表
                analysis_files = self.config_data['files']
                
            else:
                # 处理TXT格式配置文件
                analysis_files = []
                
                with open(config_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        # 跳过注释行和空行
                        if not line or line.startswith("#"):
                            continue
                        
                        # 按Tab分割行
                        parts = line.split('\t')
                        
                        # 确保至少有NGS文件路径
                        if len(parts) >= 1:
                            ngs_file = parts[0].strip()
                            chip_file = parts[1].strip() if len(parts) > 1 else ""
                            chip_type = parts[2].strip() if len(parts) > 2 else "small"
                            
                            analysis_files.append({
                                'ngs_file': ngs_file,
                                'chip_file': chip_file,
                                'chip_type': chip_type
                            })
                
                # 获取输入输出目录
                if analysis_files:
                    # 使用第一个NGS文件的目录作为输入目录
                    self.input_dir = os.path.dirname(analysis_files[0]['ngs_file'])
                    self.output_dir = os.path.join(self.input_dir, "results")
                    os.makedirs(self.output_dir, exist_ok=True)
                
                # 创建配置数据结构
                self.config_data = {
                    'input_dir': self.input_dir,
                    'output_dir': self.output_dir,
                    'files': analysis_files
                }
            
            # 验证每个文件配置
            for i, file_info in enumerate(analysis_files):
                if not file_info.get('ngs_file'):
                    print(f"警告: 文件 #{i+1} 缺少ngs_file字段")
                elif not os.path.exists(file_info['ngs_file']):
                    print(f"警告: NGS文件不存在 - {file_info['ngs_file']}")
                
                # 芯片文件是可选的
                if file_info.get('chip_file') and not os.path.exists(file_info['chip_file']):
                    print(f"警告: 芯片合成文件不存在 - {file_info['chip_file']}")
            
            print(f"成功加载配置，包含 {len(analysis_files)} 个文件的分析任务")
            return True
        
        except Exception as e:
            print(f"加载配置文件时出错: {e}")
            return False
    
    def run_analysis(self):
        """运行分析流程"""
        if not self.config_data or not self.config_data.get('files'):
            print("错误: 未加载配置或配置中没有要分析的文件")
            return False
        
        # 确保输出目录存在
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # 获取文件列表
        analysis_files = self.config_data.get('files', [])
        
        # 分析进度
        total_files = len(analysis_files)
        print(f"\n开始分析 {total_files} 个文件...")
        
        # 逐个分析文件
        for i, file_info in enumerate(analysis_files):
            ngs_file = file_info.get('ngs_file')
            chip_file = file_info.get('chip_file')
            chip_type = file_info.get('chip_type', 'small')
            
            if not ngs_file or not os.path.exists(ngs_file):
                print(f"错误: NGS文件不存在 - {ngs_file}")
                continue
                
            if not chip_file or not os.path.exists(chip_file):
                print(f"警告: 芯片合成文件不存在 - {chip_file}")
                print(f"      跳过芯片位置分析")
                chip_file = None
            
            print(f"\n[{i+1}/{total_files}] 分析文件: {os.path.basename(ngs_file)}")
            
            try:
                # 步骤1: 深度分布分析
                print("  步骤1: 执行深度分布分析...")
                depth_result_dir = analyze_oligo_ngs.analyze_oligo_ngs(ngs_file, self.output_dir)
                
                # 收集结果信息
                result_info = {
                    'ngs_file': ngs_file,
                    'chip_file': chip_file,
                    'chip_type': chip_type,
                    'depth_dir': depth_result_dir,
                    'depth_img_dir': os.path.join(depth_result_dir, 'images') if depth_result_dir else None,
                }
                
                # 查找深度分析报告
                if depth_result_dir:
                    file_name_no_ext = os.path.splitext(os.path.basename(ngs_file))[0]
                    depth_report = os.path.join(depth_result_dir, f'{file_name_no_ext}_analysis_report.html')
                    if os.path.exists(depth_report):
                        result_info['depth_report'] = depth_report
                
                # 步骤2: 芯片位置分析 (如果有芯片文件)
                if chip_file:
                    print("  步骤2: 执行芯片位置分析...")
                    # 创建芯片分析输出目录
                    chip_result_dir = os.path.join(self.output_dir, f"{os.path.splitext(os.path.basename(ngs_file))[0]}_chip_analysis")
                    if not os.path.exists(chip_result_dir):
                        os.makedirs(chip_result_dir)
                    
                    # 调用芯片位置分析
                    cmd = [
                        sys.executable,
                        os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chip_depth_heatmap.py'),
                        chip_file,
                        ngs_file,
                        '--output_dir', chip_result_dir,
                        '--chip_type', chip_type
                    ]
                    print(f"  执行命令: {' '.join(cmd)}")
                    process = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if process.returncode == 0:
                        print("  芯片位置分析完成")
                        result_info['chip_dir'] = chip_result_dir
                        result_info['chip_img_dir'] = chip_result_dir
                        
                        # 查找可能的芯片分析报告
                        for html_file in glob.glob(os.path.join(chip_result_dir, '*.html')):
                            result_info['chip_report'] = html_file
                            break
                    else:
                        print(f"  芯片位置分析失败: {process.stderr}")
                
                # 添加到结果列表
                self.analysis_results.append(result_info)
                print(f"  文件 {os.path.basename(ngs_file)} 分析完成!")
                
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"  分析过程中出错: {e}")
                continue
        
        # 生成汇总报告
        if self.analysis_results:
            self.generate_summary_report()
            return True
        else:
            print("未能成功分析任何文件，无法生成报告")
            return False
    
    def generate_summary_report(self):
        """生成分析结果汇总报告"""
        if not self.analysis_results:
            print("没有可用的分析结果，无法生成汇总报告")
            return
        
        # 报告文件路径
        report_path = os.path.join(self.output_dir, "analysis_summary.html")
        
        # 生成报告内容
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html_content = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NGS Oligopool分析结果</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 20px; }}
        h1, h2, h3 {{ color: #2c3e50; }}
        .result-card {{ border-left: 4px solid #3498db; margin-bottom: 20px; }}
        .footer {{ margin-top: 30px; text-align: center; font-size: 0.8em; color: #7f8c8d; }}
        .img-gallery {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: 10px; }}
        .img-gallery img {{ max-width: 250px; max-height: 250px; object-fit: contain; }}
    </style>
</head>
<body>
    <div class="container">
        <header class="my-4">
            <h1 class="text-center">NGS Oligopool测序分析报告</h1>
            <p class="text-center text-muted">生成时间: {now}</p>
            <hr>
        </header>

        <section id="results-summary">
            <h2>分析结果汇总</h2>
            <p>共分析了 {len(self.analysis_results)} 个测序结果文件。以下是详细结果：</p>
        """
        
        # 添加每个分析结果
        for i, result in enumerate(self.analysis_results):
            ngs_file = os.path.basename(result['ngs_file'])
            chip_type = result['chip_type']
            
            html_content += f"""
            <div class="card result-card">
                <div class="card-header">
                    <h3>文件 {i+1}: {ngs_file}</h3>
                </div>
                <div class="card-body">
                    <p><strong>芯片类型:</strong> {chip_type}</p>
            """
            
            # 添加报告链接
            if 'depth_report' in result:
                rel_path = os.path.relpath(result['depth_report'], self.output_dir)
                html_content += f'<p><a href="{rel_path}" class="btn btn-primary" target="_blank">查看深度分析报告</a></p>\n'
            
            if 'chip_report' in result:
                rel_path = os.path.relpath(result['chip_report'], self.output_dir)
                html_content += f'<p><a href="{rel_path}" class="btn btn-info" target="_blank">查看芯片位置分析报告</a></p>\n'
            
            # 添加图像预览
            html_content += '<h4 class="mt-3">深度分析图像:</h4>\n<div class="img-gallery">\n'
            
            # 添加深度分析图像
            if result.get('depth_img_dir') and os.path.exists(result['depth_img_dir']):
                img_files = glob.glob(os.path.join(result['depth_img_dir'], '*.png'))
                for img_file in img_files[:4]:  # 限制显示前4张图片
                    rel_path = os.path.relpath(img_file, self.output_dir)
                    img_name = os.path.basename(img_file)
                    html_content += f'<a href="{rel_path}" target="_blank"><img src="{rel_path}" alt="{img_name}" title="{img_name}" class="img-thumbnail"></a>\n'
            
            html_content += '</div>\n'
            
            # 添加芯片分析图像
            if result.get('chip_img_dir') and os.path.exists(result['chip_img_dir']):
                html_content += '<h4 class="mt-3">芯片位置分析图像:</h4>\n<div class="img-gallery">\n'
                
                # 添加芯片位置图像
                img_files = glob.glob(os.path.join(result['chip_img_dir'], '*heatmap*.png'))
                for img_file in img_files[:4]:  # 限制显示前4张图片
                    rel_path = os.path.relpath(img_file, self.output_dir)
                    img_name = os.path.basename(img_file)
                    html_content += f'<a href="{rel_path}" target="_blank"><img src="{rel_path}" alt="{img_name}" title="{img_name}" class="img-thumbnail"></a>\n'
                
                html_content += '</div>\n'
            
            html_content += """
                </div>
            </div>
            """
        
        # 完成HTML
        html_content += f"""
        </section>
        
        <footer class="footer">
            <p>NGS Oligopool测序分析报告 | 生成时间: {now}</p>
        </footer>
    </div>
</body>
</html>
"""
        
        # 保存报告
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"\n分析汇总报告已生成: {report_path}")
        
        # 尝试打开报告
        try:
            webbrowser.open(f"file://{report_path}")
            print("已在浏览器中打开报告")
        except Exception as e:
            print(f"无法打开浏览器: {e}")
            print(f"请手动打开报告文件: {report_path}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='NGS Oligopool分析工作流')
    subparsers = parser.add_subparsers(dest='command', help='命令')
    
    # 扫描命令
    scan_parser = subparsers.add_parser('scan', help='扫描目录并创建配置文件')
    scan_parser.add_argument('input_dir', nargs='?', help='输入目录路径，包含测序结果文件')
    scan_parser.add_argument('--output_dir', '-o', help='输出目录路径 (默认为input_dir/results)')
    
    # 分析命令
    analyze_parser = subparsers.add_parser('analyze', help='使用配置文件执行分析')
    analyze_parser.add_argument('--config', '-c', required=True, help='分析配置文件路径')
    
    # 配置参数
    parser.add_argument('--config', '-c', help='分析配置文件路径')
    parser.add_argument('--gui', '-g', action='store_true', help='使用图形界面选择目录')
    parser.add_argument('input_dir', nargs='?', help='输入目录路径 (如果不使用配置文件)')
    
    args = parser.parse_args()
    
    # 创建工作流对象
    workflow = NGSOligoAnalysisWorkflow()
    
    # 解析命令
    if args.command == 'scan':
        # 如果没有指定输入目录，使用弹窗选择
        input_dir = args.input_dir
        if not input_dir:
            print("打开文件夹选择对话框...")
            root = tk.Tk()
            root.withdraw()  # 隐藏主窗口
            input_dir = filedialog.askdirectory(title="选择包含测序结果的目录")
            root.destroy()
            if not input_dir:
                print("未选择目录，退出程序")
                return
        
        output_dir = args.output_dir
        workflow.scan_directory(input_dir)
    elif args.command == 'analyze' or args.config:
        config_file = args.config
        if workflow.load_config(config_file):
            workflow.run_analysis()
    elif args.input_dir:
        # 使用命令行指定的目录
        workflow.scan_directory(args.input_dir)
    else:
        # 默认行为：弹窗选择目录
        print("打开文件夹选择对话框...")
        root = tk.Tk()
        root.withdraw()  # 隐藏主窗口
        input_dir = filedialog.askdirectory(title="选择包含测序结果的目录")
        root.destroy()
        if input_dir:
            print(f"已选择目录: {input_dir}")
            workflow.scan_directory(input_dir)
        else:
            print("未选择目录，退出程序")
            parser.print_help()

if __name__ == "__main__":
    main() 