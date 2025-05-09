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

# 用于创建一个独立的文件选择进程
FILE_SELECTOR_SCRIPT = """
import os
import sys
import tkinter as tk
from tkinter import filedialog

# 创建一个独立的文件选择对话框
def select_file():
    # 创建根窗口但不显示
    root = tk.Tk()
    root.withdraw()
    
    # 打开文件选择对话框
    file_path = filedialog.askopenfilename(
        title="选择芯片合成文件",
        filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
    )
    
    # 销毁根窗口
    root.destroy()
    
    return file_path

if __name__ == "__main__":
    temp_file = sys.argv[1]
    
    # 获取文件路径
    file_path = select_file()
    
    # 将文件路径写入临时文件
    with open(temp_file, 'w', encoding='utf-8') as f:
        f.write(file_path)
"""

# 创建独立的文件选择器脚本
def create_file_selector_script():
    script_path = os.path.join(tempfile.gettempdir(), "file_selector.py")
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(FILE_SELECTOR_SCRIPT)
    return script_path

# 创建分析UI HTML页面
def create_analysis_ui(input_dir, output_dir=None):
    """创建分析用户界面HTML文件"""
    if output_dir is None:
        output_dir = os.path.join(input_dir, "results")
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 查找目录中的所有90.oligostat.xls文件
    ngs_files = []
    for ext in ['*90.oligostat.xls', '*90.oligostat.xlsx']:
        ngs_files.extend(glob.glob(os.path.join(input_dir, ext)))
    
    # 按修改时间排序
    ngs_files.sort(key=os.path.getmtime, reverse=True)
    
    if not ngs_files:
        print(f"在目录 {input_dir} 中未找到任何 *90.oligostat.xls 文件")
        return None
    
    print(f"找到 {len(ngs_files)} 个测序结果文件")
    
    # 创建文件选择器脚本
    file_selector = create_file_selector_script()
    print(f"文件选择器脚本已创建: {file_selector}")
    
    # 生成HTML
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_content = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NGS Oligopool分析工作流</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 20px; }}
        h1, h2, h3 {{ color: #2c3e50; }}
        .file-item {{ border-left: 4px solid #3498db; }}
        .footer {{ margin-top: 30px; text-align: center; font-size: 0.8em; color: #7f8c8d; }}
        #status-area {{ display: none; max-height: 300px; overflow-y: auto; }}
        .btn-outline-secondary:hover {{ background-color: #e9ecef; }}
    </style>
</head>
<body>
    <div class="container">
        <header class="my-4">
            <h1 class="text-center">NGS Oligopool测序分析工作流</h1>
            <p class="text-center text-muted">扫描目录: {input_dir} | 生成时间: {now}</p>
            <hr>
        </header>

        <section id="file-list">
            <h2>检测到的文件 <small class="text-muted">({len(ngs_files)}个文件)</small></h2>
            <p>请为每个测序文件指定对应的芯片合成文件路径，并选择芯片类型。点击"浏览"按钮可以打开文件选择对话框。</p>
            
            <form id="analysis-form" action="javascript:void(0);">
                <input type="hidden" name="input_dir" value="{input_dir}">
                <input type="hidden" name="output_dir" value="{output_dir}">
                
                """
    
    # 添加文件条目
    for i, ngs_file in enumerate(ngs_files):
        file_name = os.path.basename(ngs_file)
        html_content += f"""
            <div class="file-item card mb-3">
                <div class="card-body">
                    <h5 class="card-title"><span class="text-primary">{i+1}.</span> {file_name}</h5>
                    <p class="card-text text-muted">{file_name}</p>
                    <div class="row mb-2">
                        <div class="col-md-6">
                            <label class="form-label">芯片合成文件路径:</label>
                            <div class="input-group">
                                <input type="text" class="form-control chip-file" name="chip_file_{i}" id="chip_file_{i}" placeholder="选择对应的芯片合成文件...">
                                <button class="btn btn-outline-secondary browse-btn" type="button" data-index="{i}">浏览...</button>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <label class="form-label">芯片类型:</label>
                            <select class="form-select chip-type" name="chip_type_{i}">
                                <option value="small" selected>小片(540x635)</option>
                                <option value="large">大片(1080x636)</option>
                            </select>
                        </div>
                        <div class="col-md-3">
                            <label class="form-label">添加到分析:</label>
                            <div class="form-check mt-2">
                                <input class="form-check-input file-select" type="checkbox" name="select_{i}" checked>
                                <label class="form-check-label">选择此文件</label>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            """
    
    # 添加提交按钮和剩余的HTML
    html_content += f"""
                <div class="d-grid gap-2 col-md-6 mx-auto mt-4">
                    <button type="button" id="start-analysis" class="btn btn-primary btn-lg">开始分析</button>
                </div>
            </form>
        </section>
        
        <section id="status-area" class="my-4">
            <h3>分析状态</h3>
            <div class="progress mb-3">
                <div id="progress-bar" class="progress-bar progress-bar-striped progress-bar-animated" role="progressbar" style="width: 0%"></div>
            </div>
            <div class="card">
                <div class="card-body">
                    <pre id="status-log" class="mb-0"></pre>
                </div>
            </div>
        </section>
        
        <section id="results-section" class="my-4" style="display: none;">
            <h3>分析结果</h3>
            <div id="results-list" class="list-group">
                <!-- 结果将在这里动态添加 -->
            </div>
        </section>
        
        <footer class="footer">
            <p>NGS Oligopool测序分析工作流 | 生成时间: {now}</p>
        </footer>
    </div>
    
    <script>
        // 文件选择脚本路径
        const fileSelector = "{file_selector.replace('\\\\', '\\\\\\\\')}";
        
        // 浏览按钮点击事件处理
        document.addEventListener('DOMContentLoaded', function() {{
            // 调试信息
            console.log("页面加载完成，初始化事件处理...");
            
            // 为所有浏览按钮添加点击事件
            document.querySelectorAll('.browse-btn').forEach(button => {{
                button.addEventListener('click', function() {{
                    const index = this.getAttribute('data-index');
                    const inputField = document.getElementById(`chip_file_${{index}}`);
                    
                    // 显示正在选择文件的状态
                    this.disabled = true;
                    this.textContent = "选择中...";
                    inputField.placeholder = "正在打开文件选择对话框...";
                    
                    // 创建一个唯一的标识符，用于临时文件
                    const tempId = Date.now().toString() + Math.random().toString(36).substr(2, 5);
                    
                    // 创建嵌入式iframe来调用系统的文件选择对话框
                    // 这里我们使用一个iframe来提交表单，模拟调用外部程序
                    const formHtml = `
                    <form id="fileForm_${{tempId}}" action="app://select-file/" method="GET" target="_blank">
                        <input type="hidden" name="id" value="${{tempId}}">
                        <input type="hidden" name="script" value="${{fileSelector}}">
                    </form>
                    `;
                    
                    const container = document.createElement('div');
                    container.style.display = 'none';
                    container.innerHTML = formHtml;
                    document.body.appendChild(container);
                    
                    // 提交表单
                    setTimeout(() => {{
                        // 重置按钮状态
                        this.disabled = false;
                        this.textContent = "浏览...";
                        inputField.placeholder = "选择对应的芯片合成文件...";
                        
                        // 提示用户手动输入
                        alert("请手动运行文件选择程序:\\n" + 
                              "1. 按下Ctrl+C终止当前程序\\n" +
                              "2. 重新运行以下命令选择文件:\\n" +
                              `   python "${{fileSelector}}" "temp_file_${{index}}.txt"\\n` +
                              "3. 然后将选择的文件路径复制到此输入框中");
                    }}, 1000);
                }});
            }});
            
            // 开始分析按钮点击事件
            const startButton = document.getElementById('start-analysis');
            console.log("找到开始分析按钮:", startButton ? "是" : "否");
            
            startButton.addEventListener('click', function(event) {{
                console.log("开始分析按钮被点击");
                event.preventDefault(); // 阻止默认行为
                
                try {{
                    // 显示加载状态
                    startButton.disabled = true;
                    startButton.textContent = "处理中...";
                    
                    // 显示状态区域
                    document.getElementById('status-area').style.display = 'block';
                    const statusLog = document.getElementById('status-log');
                    statusLog.textContent = '正在收集分析配置...\\n';
                    
                    // 收集表单数据
                    const form = document.getElementById('analysis-form');
                    const formData = new FormData(form);
                    const analysisData = {{}};
                    
                    analysisData.input_dir = formData.get('input_dir');
                    analysisData.output_dir = formData.get('output_dir');
                    analysisData.files = [];
                    
                    statusLog.textContent += '读取目录信息完成\\n';
                    
                    // 收集文件配置
                    const fileCount = {len(ngs_files)};
                    statusLog.textContent += `处理 ${{fileCount}} 个文件的配置...\\n`;
                    
                    // 文件路径数组
                    const ngsFiles = [
                        {", ".join([f"'{ngs_file.replace('\\', '\\\\')}'" for ngs_file in ngs_files])}
                    ];
                    
                    for (let i = 0; i < fileCount; i++) {{
                        if (formData.get(`select_${{i}}`) === "on") {{
                            const ngsFile = ngsFiles[i];
                            const chipFile = formData.get(`chip_file_${{i}}`);
                            const chipType = formData.get(`chip_type_${{i}}`);
                            
                            analysisData.files.push({{
                                ngs_file: ngsFile,
                                chip_file: chipFile,
                                chip_type: chipType
                            }});
                            
                            statusLog.textContent += `  文件 ${{i+1}}: ${{chipType}} 芯片类型\\n`;
                        }}
                    }}
                    
                    statusLog.textContent += '配置数据收集完成\\n';
                    
                    // 生成JSON配置
                    const configJson = JSON.stringify(analysisData, null, 2);
                    statusLog.textContent += '生成JSON配置完成\\n';
                    
                    // 方法1: 尝试使用下载链接保存文件
                    try {{
                        const configElement = document.createElement('a');
                        configElement.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(configJson));
                        configElement.setAttribute('download', 'analysis_config.json');
                        configElement.style.display = 'none';
                        document.body.appendChild(configElement);
                        
                        statusLog.textContent += '正在尝试通过浏览器下载配置文件...\\n';
                        configElement.click();
                        document.body.removeChild(configElement);
                        
                        statusLog.textContent += '配置文件应已下载\\n';
                    }} catch (downloadError) {{
                        statusLog.textContent += `浏览器下载方法失败: ${{downloadError.message}}\\n`;
                        statusLog.textContent += '尝试备用方法...\\n';
                    }}
                    
                    // 方法2: 在页面上显示配置供复制
                    statusLog.textContent += '将配置显示在页面上以便手动复制保存...\\n';
                    
                    // 创建配置文本区
                    const configArea = document.createElement('div');
                    configArea.className = 'mt-4';
                    configArea.innerHTML = `
                    <div class="card">
                        <div class="card-header bg-warning">
                            <h5 class="mb-0">配置数据 - 请复制并保存为 analysis_config.json</h5>
                        </div>
                        <div class="card-body">
                            <button class="btn btn-sm btn-primary mb-2" id="copy-config">复制配置</button>
                            <pre style="max-height: 300px; overflow: auto; background: #f8f9fa; padding: 10px; border-radius: 5px;">${{configJson}}</pre>
                        </div>
                    </div>
                    `;
                    
                    const statusArea = document.getElementById('status-area');
                    statusArea.appendChild(configArea);
                    
                    // 添加复制按钮功能
                    document.getElementById('copy-config').addEventListener('click', function() {{
                        navigator.clipboard.writeText(configJson).then(function() {{
                            this.textContent = '复制成功!';
                            setTimeout(() => {{ this.textContent = '复制配置'; }}, 2000);
                        }}.bind(this)).catch(function(err) {{
                            alert('复制失败: ' + err);
                        }});
                    }});
                    
                    // 恢复按钮状态
                    setTimeout(() => {{
                        startButton.disabled = false;
                        startButton.textContent = "开始分析";
                    }}, 2000);
                    
                    // 提醒用户如何继续
                    statusLog.textContent += '\\n配置已保存，请运行分析脚本。\\n';
                    statusLog.textContent += '执行命令: python ngs_oligo_analysis_workflow.py --config analysis_config.json\\n';
                    
                    alert('配置已保存 (或显示在页面上供复制)。\\n\\n请将该文件放在与分析脚本相同的目录，然后在命令行中运行:\\npython ngs_oligo_analysis_workflow.py --config analysis_config.json');
                }} catch (error) {{
                    console.error("处理配置时出错:", error);
                    alert(`处理配置时出错: ${{error.message || error}}\\n请检查控制台获取更多信息`);
                    
                    // 恢复按钮状态
                    startButton.disabled = false;
                    startButton.textContent = "开始分析";
                    
                    // 显示错误信息
                    const statusLog = document.getElementById('status-log');
                    statusLog.textContent += `\\n错误: ${{error.message || error}}\\n`;
                    statusLog.textContent += '请尝试手动创建配置文件，或联系开发人员获取帮助。\\n';
                }}
            }});
        }});
    </script>
</body>
</html>
"""
    
    # 保存HTML文件
    html_path = os.path.join(output_dir, "analysis_ui.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return html_path

class NGSOligoAnalysisWorkflow:
    """NGS Oligopool分析工作流"""
    
    def __init__(self):
        """初始化分析工作流"""
        self.input_dir = None
        self.output_dir = None
        self.config_data = None
        self.analysis_results = []
    
    def scan_directory(self, directory_path):
        """扫描目录获取可分析文件"""
        if not os.path.exists(directory_path):
            print(f"错误: 目录不存在 - {directory_path}")
            return False
        
        # 设置输入目录
        self.input_dir = directory_path
        self.output_dir = os.path.join(directory_path, "results")
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 创建Web界面
        html_path = create_analysis_ui(self.input_dir, self.output_dir)
        
        if html_path:
            print(f"分析界面已创建: {html_path}")
            # 打开浏览器显示界面
            try:
                webbrowser.open(f"file://{html_path}")
                print("已在浏览器中打开分析界面")
            except Exception as e:
                print(f"无法打开浏览器: {e}")
                print(f"请手动打开文件: {html_path}")
            
            return True
        else:
            print("未能创建分析界面")
            return False
    
    def load_config(self, config_file):
        """从配置文件加载分析设置"""
        if not os.path.exists(config_file):
            print(f"错误: 配置文件不存在 - {config_file}")
            return False
        
        try:
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
            
            # 验证每个文件配置
            for i, file_info in enumerate(self.config_data['files']):
                if not file_info.get('ngs_file'):
                    print(f"警告: 文件 #{i+1} 缺少ngs_file字段")
                elif not os.path.exists(file_info['ngs_file']):
                    print(f"警告: NGS文件不存在 - {file_info['ngs_file']}")
                
                # 芯片文件是可选的
                if file_info.get('chip_file') and not os.path.exists(file_info['chip_file']):
                    print(f"警告: 芯片合成文件不存在 - {file_info['chip_file']}")
            
            print(f"成功加载配置，包含 {len(self.config_data['files'])} 个文件的分析任务")
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
    scan_parser = subparsers.add_parser('scan', help='扫描目录并创建分析界面')
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