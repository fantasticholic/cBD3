# Sequence Evolution.py 脚本改进文档

## 改进概述

本文档记录了对 `Sequence Evolution.py` 脚本的主要改进，重点增强了交叉验证评估部分，提供更详细的模型性能指标。

## 主要改进

### 1. 增强的交叉验证评估系统

- **替换了简单的 cross_val_score 调用**，改为手动实现的交叉验证循环，以获取更详细的评估指标
- **增加了多维度性能指标**：
  - 精确率 (Precision)
  - 召回率 (Recall)  
  - F1 分数 (F1-Score)
  - 混淆矩阵 (Confusion Matrix) 详细信息
  - 完整分类报告 (Classification Report)

### 2. 详细的折叠级评估

- 对每个交叉验证折叠单独训练模型并计算性能指标
- 显示每个折叠的关键指标（准确率、AUC、F1分数）
- 计算并显示所有折叠的平均指标及其标准差

### 3. 改进的数据分布可视化

- 添加了训练数据类别分布显示（Label 0 和 Label 1 的样本数量）
- 提供特征维度信息的更完整显示

### 4. 增强的输出格式

- 重新格式化的输出，使用分隔线增强可读性
- 为不同阶段的输出添加了清晰的标题和组织
- 改进了训练集整体表现的指标显示，添加了 AUC 和 F1 分数

### 5. 增加的导入项

- 添加了新的 scikit-learn 评估指标导入：
  ```python
  from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix, classification_report
  ```

## 如何观察改进效果

1. **首次运行** - 当脚本训练新模型时，将显示完整的交叉验证详细信息
2. **后续运行** - 加载现有模型时，将显示增强的训练集表现统计
3. **进化输出** - 进化搜索部分的输出格式也经过了优化，添加了分隔线和标题

## 技术细节

### 交叉验证实现

- 使用手动循环替代 cross_val_score，为每个折叠创建独立的 pipeline
- 使用不同的随机种子确保结果可复现同时增加多样性
- 收集所有预测结果以生成整体分类报告和混淆矩阵

### 性能指标计算

- 使用 pandas DataFrame 存储折叠结果，便于统计分析
- 计算平均值和标准差以评估模型稳定性
- 生成针对标签 0（低活性肽）和标签 1（高活性肽）的详细分类报告

## 使用说明

脚本使用方式保持不变，直接运行即可：

```bash
python "Sequence Evolution.py"
```

如需强制重新训练模型以查看交叉验证详情，请先删除 `rf_amp_model.joblib` 文件：

```bash
# Windows PowerShell
Remove-Item rf_amp_model.joblib

# 然后重新运行脚本
python "Sequence Evolution.py"
```