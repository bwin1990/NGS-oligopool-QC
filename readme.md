# NGS Oligopool 质量控制分析工具

## 项目概述

本项目是一套用于分析NGS oligopool测序数据的综合质量控制工具，提供深度分布分析、碱基组成关系分析和芯片位置热力图分析功能。

## 数据格式要求

- **输入数据文件**：以`90.oligostat.xls`或`90.oligostat.xlsx`结尾的文件
- **示例文件名**：`NovaS05933-RD-YF-Nova-20250311-WYJ-2haoji-190.oligostat.xls`
- **必需字段**：
  - `Oligo seq`：DNA序列信息
  - `Depth`：对应的测序深度
  - 其他相关内容

## 核心功能模块

### 1. 深度分布与碱基组成分析 (`analyze_oligo_ngs.py`)

**分析内容：**
- **深度分布分析**：
  - 生成测序深度分布直方图
  - 计算累积分布函数(CDF)曲线
  - 计算95%/5%深度比值，评估覆盖均匀性
  - 提供质量等级评估（优/良/一般/较差）

- **碱基组成关系分析**：
  - 分析深度与各碱基（A、T、G、C）含量的相关性
  - 分析深度与GC含量的关系
  - 分析深度与AG含量的关系
  - 生成散点图和相关系数

- **输出结果**：
  - HTML格式的详细分析报告
  - 碱基组成数据CSV文件
  - 深度统计数据CSV文件
  - 相关性矩阵CSV文件
  - 高质量PNG图像文件

### 2. 芯片位置热力图分析 (`chip_depth_heatmap.py`)

**分析内容：**
- **芯片规格支持**：
  - Small芯片：540×635像素
  - Large芯片：1080×636像素

- **热力图类型**：
  - **标准深度分布热力图**：显示整体深度分布
  - **低覆盖度强调热力图**：突出显示低于平均深度50%的区域
  - **低覆盖度聚焦热力图**：仅显示低覆盖度区域
  - **高覆盖度区域热力图**：显示超过平均深度2倍的区域

- **位置偏好性分析**：
  - X轴位置与深度关系分析
  - Y轴位置与深度关系分析
  - 不同区域（上下左右、中心）统计分析

- **输出结果**：
  - 芯片深度矩阵CSV文件
  - 多类型热力图PNG文件
  - 位置偏好性分析图
  - 区域统计数据CSV文件

### 3. 芯片Indel热力图分析 (`chip_indel_heatmap.py`)

**分析内容：**
- **Indel分类分析**：
  - Small Indel (<3bp del)：小于3bp缺失的平均indel率
  - Large Indel (≥3bp del)：大于等于3bp缺失的平均indel率

- **热力图类型**：
  - **小缺失热力图**：显示<3bp缺失的分布（YlOrRd色彩）
  - **大缺失热力图**：显示≥3bp缺失的分布（Reds色彩）
  - **并排对比热力图**：两种indel类型的直观对比
  - **高indel区域强调图**：使用对比色和轮廓线突出显示>1.5x均值的区域
  - **聚焦高indel区域图**：仅显示高indel区域，类似深度分析中的低覆盖度聚焦图

- **位置偏好性分析**：
  - X轴、Y轴位置与indel率的相关性分析
  - 不同区域（上下左右、中心）的indel统计

- **输出结果**：
  - 小缺失/大缺失矩阵CSV文件
  - 多类型indel热力图PNG文件（包含强调图和聚焦图）
  - 高indel区域轮廓线标注（类似深度分析的低覆盖度强调）
  - 位置偏好性分析图
  - 区域统计数据CSV文件

### 4. 批量分析工作流 (`ngs_oligo_analysis_workflow.py`)

**功能特点：**
- **智能文件扫描**：自动扫描目录中的所有`*90.oligostat.xls`文件
- **配置文件生成**：自动生成配置文件供用户编辑
- **批量处理**：支持多个文件的批量分析
- **三重分析**：集成深度分析、芯片位置分析和芯片indel分析
- **汇总报告**：生成包含所有分析结果的HTML汇总报告

## 使用方法

### 方法1：一站式分析（推荐）

```bash
# 使用图形界面选择目录
python ngs_oligo_analysis_workflow.py

# 或直接指定目录
python ngs_oligo_analysis_workflow.py oneshot /path/to/data/directory
```

### 方法2：分步骤分析

```bash
# 步骤1：扫描目录并生成配置文件
python ngs_oligo_analysis_workflow.py scan /path/to/data/directory

# 步骤2：编辑生成的配置文件（添加芯片合成文件路径）

# 步骤3：运行分析
python ngs_oligo_analysis_workflow.py analyze --config /path/to/config.txt
```

### 方法3：单文件分析

```bash
# 仅深度分析
python analyze_oligo_ngs.py input_file.xls [output_directory]

# 仅芯片位置分析
python chip_depth_heatmap.py chip_synthesis_file.txt ngs_file.xls --chip_type small

# 仅芯片indel分析
python chip_indel_heatmap.py chip_synthesis_file.txt ngs_file.xls --chip_type small
```

## 配置文件格式

配置文件采用Tab分隔的文本格式：
```
# 格式：NGS文件路径 [Tab] 芯片合成文件路径 [Tab] 芯片类型
/path/to/ngs_file1.xls	/path/to/chip_file1.txt	small
/path/to/ngs_file2.xls	/path/to/chip_file2.txt	large
```

## 输出结构

```
results/
├── analysis_summary.html                    # 汇总报告
├── file1_name/                             # 单个文件分析结果
│   ├── file1_name_analysis_report.html     # 深度分析报告
│   ├── images/                             # 图像文件
│   │   ├── depth_distribution.png          # 深度分布图
│   │   ├── depth_vs_bases.png             # 深度与碱基关系图
│   │   └── depth_vs_gc_ag.png             # 深度与GC/AG关系图
│   └── data/                              # 数据文件
│       ├── nucleotide_composition.csv      # 碱基组成数据
│       ├── depth_statistics.csv           # 深度统计数据
│       └── correlation_matrix.csv         # 相关性矩阵
├── file1_name_chip_analysis/              # 芯片位置分析结果
│   ├── chip_depth_matrix.csv              # 芯片深度矩阵
│   ├── file1_name_depth_heatmap_standard.png
│   ├── file1_name_depth_heatmap_low_emph.png
│   └── file1_name_region_statistics.csv   # 区域统计数据
├── file1_name_chip_indel_analysis/        # 芯片indel分析结果
│   ├── chip_small_indel_matrix.csv        # 小缺失矩阵
│   ├── chip_large_indel_matrix.csv        # 大缺失矩阵
│   ├── file1_name_small_indel_heatmap.png # 小缺失热力图
│   ├── file1_name_large_indel_heatmap.png # 大缺失热力图
│   ├── file1_name_combined_indel_heatmap.png # 对比热力图
│   ├── file1_name_high_indel_emphasis.png # 高indel区域强调图
│   ├── file1_name_focused_high_indel_regions.png # 聚焦高indel区域图
│   ├── file1_name_indel_position_analysis.png # 位置分析图
│   └── file1_name_*_indel_region_statistics.csv # 区域统计
└── analysis_config.txt                    # 分析配置文件
```

## 质量评估指标

- **95%/5%深度比值**：
  - < 3.0：优秀（覆盖度非常均匀）
  - 3.0-4.5：良好（覆盖度较为均匀）
  - 4.5-6.5：一般（存在一定差异）
  - > 6.5：较差（差异较大，可能存在偏好性扩增）

## 依赖环境

- Python 3.7+
- pandas
- numpy
- matplotlib
- seaborn
- tqdm
- tkinter（用于图形界面）

## 安装依赖

```bash
pip install pandas numpy matplotlib seaborn tqdm
```

## 注意事项

1. 确保输入文件包含必需列：
   - `Oligo seq`：DNA序列
   - `Depth`：测序深度
   - `Mean indel(<3bp del)`：小缺失indel率（用于indel分析）
   - `Mean indel`：总indel率（用于indel分析）
2. 芯片合成文件为可选，如无此文件将跳过芯片位置和indel分析
3. 支持Excel (.xls, .xlsx)、CSV和TSV格式的输入文件
4. 建议使用一站式分析模式以获得完整的分析结果
5. 生成的HTML报告可直接在浏览器中查看
6. Indel分析会自动计算大缺失率（Mean indel - Mean indel(<3bp del)）

