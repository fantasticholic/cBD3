# 抗菌肽优化项目

## 项目概述
本项目用于抗菌肽的优化研究，包含从数据准备、结构分析、模型训练到结果可视化的完整工作流程。

## 文件夹结构

```
peptide_optimization/
├── 01_data/                 # 数据目录
│   ├── raw/                # 原始未处理数据
│   ├── processed/          # 处理后的数据集
│   └── reference/          # 参考序列和数据库
├── 02_analysis/             # 分析目录
│   ├── structure/          # 结构分析结果
│   ├── heliquest/          # HeliQuest分析（原heliquest_analysis）
│   └── physico_chemical/   # 物理化学特性分析
├── 03_models/               # 模型目录
│   ├── training/           # 模型训练相关文件
│   └── evaluation/         # 模型评估结果
├── 04_results/              # 结果目录
│   ├── summaries/          # 结果摘要和统计
│   └── figures/            # 可视化图表和图形
├── 05_scripts/              # 脚本目录
│   ├── preprocessing/      # 数据预处理脚本
│   ├── analysis/           # 分析脚本
│   └── visualization/      # 可视化脚本
├── 06_documentation/        # 文档目录
│   ├── protocols/          # 实验协议和方法
│   └── reports/            # 研究报告和论文
└── venv/                    # Python虚拟环境
```

## 目录用途说明

### 01_data/
- **raw/**: 存储原始肽序列数据、实验数据等未处理的输入数据
- **processed/**: 存储经过清洗、标准化、特征提取后的数据集
- **reference/**: 存储参考序列库、数据库和标准数据集

### 02_analysis/
- **structure/**: 存储肽结构预测、3D建模和结构分析结果
- **heliquest/**: 存储HeliQuest软件分析结果，包括螺旋轮图、疏水性分析等
- **physico_chemical/**: 存储物理化学特性分析，如电荷、疏水性、等电点等

### 03_models/
- **training/**: 存储模型定义、训练代码、超参数配置和训练过程文件
- **evaluation/**: 存储模型性能评估结果、预测输出和验证报告

### 04_results/
- **summaries/**: 存储结果汇总表格、统计分析和重要发现
- **figures/**: 存储所有可视化图表、数据展示和发表级图像

### 05_scripts/
- **preprocessing/**: 数据清洗、转换和特征工程脚本
- **analysis/**: 数据分析、建模和评估脚本
- **visualization/**: 结果可视化和绘图脚本

### 06_documentation/
- **protocols/**: 实验方法、分析流程和标准操作程序
- **reports/**: 研究报告、演示文稿和论文草稿

## 关于虚拟环境（venv）
项目中包含的venv目录是Python虚拟环境，用于隔离项目依赖。使用时请激活此虚拟环境以确保使用正确的库版本。

```bash
# Windows
venv\Scripts\activate

# Linux/Mac
source venv/bin/activate
```

## 使用指南
1. 将原始数据放入`01_data/raw/`目录
2. 使用`05_scripts/preprocessing/`中的脚本处理数据
3. 将处理后的数据保存到`01_data/processed/`
4. 使用`05_scripts/analysis/`中的脚本进行分析
5. 将分析结果保存到对应目录
6. 使用`05_scripts/visualization/`中的脚本生成可视化
7. 将图表保存到`04_results/figures/`

## 注意事项
- 请按照目录结构组织文件，确保项目的可维护性
- 重要数据和结果请及时备份
- 大型文件或敏感数据请勿提交到版本控制系统

## 联系人
如有问题，请联系项目负责人。