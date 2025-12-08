[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Science](https://img.shields.io/badge/Science-Bioinformatics-brightgreen.svg)]()

## 项目简介
cBD3-DeepEvo-Pept是一个基于机器学习与进化算法的抗菌肽设计与优化平台。该项目整合了多种特征提取方法、随机森林分类器和先进的遗传算法，用于设计和优化具有高抗菌活性的多肽序列。平台不仅支持常规序列进化，还提供了带约束条件的进化算法，确保生成的肽序列具有更好的溶解性和安全性。

## 主要特性

- **多维度特征提取**：包含氨基酸组成、二肽频率、Kyte-Doolittle疏水性、疏水矩、螺旋倾向性和Biopython理化性质等多种特征
- **强大的分类模型**：基于随机森林的高效分类器，支持交叉验证和性能评估
- **高级进化算法**：包含种群进化、锦标赛选择、交叉、多点突变、精英保留和种群多样性维护
- **约束优化**：支持硬约束（如禁用特定氨基酸）、疏水性惩罚和两亲性交替引导
- **可视化分析**：提供详细的模型评估指标和进化过程可视化
- **模型持久化**：支持模型保存和加载，方便后续使用

## 项目结构

```
cBD3-DeepEvo-Pept/
├── cBD3/                      # 核心代码与数据目录
│   ├── Sequence Evolution.py        # 基础进化算法实现
│   ├── Sequence Evolution-v2.py     # 带约束条件的进化算法
│   ├── process_amp_data.py          # 抗菌肽数据处理
│   ├── final_train_set.csv          # 训练数据集
│   ├── rf_amp_model.joblib          # 训练好的模型
│   └── ...                           # 其他辅助文件
├── peptide_optimization/       # 结构化项目目录
│   ├── 01_data/                 # 数据目录
│   │   ├── raw/                # 原始数据
│   │   ├── processed/          # 处理后数据
│   │   └── reference/          # 参考数据
│   ├── 02_analysis/             # 分析目录
│   │   ├── structure/          # 结构分析
│   │   ├── heliquest/          # HeliQuest分析
│   │   └── physico_chemical/   # 物理化学特性
│   ├── 03_models/               # 模型目录
│   │   ├── training/           # 模型训练
│   │   └── evaluation/         # 模型评估
│   ├── 04_results/              # 结果目录
│   ├── 05_scripts/              # 脚本目录
│   └── 06_documentation/        # 文档目录
├── candidates.csv              # 进化算法生成的候选序列
├── candidates_v2.csv            # 约束进化算法生成的候选序列
├── excel_to_csv.py             # Excel转CSV工具
├── validate_amp_dataset.py     # 数据集验证工具
└── README.md                   # 本文件
```

## 安装说明

### 系统要求

- Python 3.8 或更高版本
- 至少 4GB RAM（推荐 8GB 或更多）
- 2GB 可用磁盘空间

### 依赖安装

1. 克隆仓库：
```bash
git clone https://github.com/yourusername/cBD3-DeepEvo-Pept.git
cd cBD3-DeepEvo-Pept
```

2. 创建虚拟环境（推荐）：
```bash
# 使用venv创建虚拟环境
python -m venv venv

# 激活虚拟环境
# Windows
venv\Scripts\activate
# Linux/Mac
source venv/bin/activate
```

3. 安装依赖：
```bash
pip install numpy pandas scikit-learn biopython joblib matplotlib seaborn
```

或者使用requirements.txt（如果提供）：
```bash
pip install -r requirements.txt
```

## 使用指南

### 1. 数据准备

确保您的训练数据格式为CSV，包含以下列：
- `Sequence`: 肽序列（大写字母，不含空格）
- `Label`: 标签（0表示低活性，1表示高活性）

### 2. 基础进化算法

运行基础进化算法：
```bash
cd cBD3
python "Sequence Evolution.py"
```

该脚本将：
- 训练随机森林模型（如果不存在已训练的模型）
- 进行交叉验证评估
- 运行进化算法生成优化的肽序列
- 保存前50个候选序列到`candidates.csv`

### 3. 约束进化算法

运行带约束条件的进化算法（推荐）：
```bash
cd cBD3
python "Sequence Evolution-v2.py"
```

该算法具有以下特点：
- 禁止使用半胱氨酸(C)残基
- 惩罚连续3个或更多强疏水残基的序列
- 控制整体疏水性比例
- 引导两亲性交替模式

### 4. 数据验证

验证AMP数据集质量：
```bash
python validate_amp_dataset.py
```

### 5. Excel数据转换

将Excel格式的肽数据转换为CSV：
```bash
python excel_to_csv.py
```

## 配置选项

### 进化算法参数

可以在脚本中调整以下参数：

```python
# 种群与进化参数
POPULATION_SIZE = 120      # 种群大小
GENERATIONS = 80           # 进化代数
ELITE_SIZE = 10            # 精英个体数量
TOURNAMENT_SIZE = 5        # 锦标赛选择大小
MUTATION_RATE = 0.12       # 突变概率
CROSSOVER_RATE = 0.6       # 交叉概率

# 约束参数（v2版本）
HYDROPHOBIC_STRONG = set(['L', 'V', 'I', 'F', 'W'])  # 强疏水残基
ALLOWED_AA = ['A','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']  # 允许的氨基酸
```

### 模型参数

```python
# 随机森林参数
RandomForestClassifier(
    n_estimators=200,                    # 树的数量
    class_weight='balanced_subsample',   # 类别权重
    random_state=RANDOM_SEED,           # 随机种子
    n_jobs=-1                           # 并行核心数
)
```

## 使用示例

### 自定义数据训练和进化

```python
from Sequence_Evolution import train_model, run_evolution
import pandas as pd

# 加载自定义数据
custom_data = pd.read_csv('my_peptide_data.csv')

# 训练模型
model = train_model('my_peptide_data.csv', 'my_model.joblib')

# 运行进化
custom_seed = "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK"
candidates = run_evolution(model, seed_seq=custom_seed, save_candidates='my_candidates.csv')
```

### 评估已训练模型

```python
from Sequence_Evolution import load_model, get_features
import pandas as pd

# 加载模型
model = load_model('rf_amp_model.joblib')

# 测试新序列
test_seq = "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK"
features = get_features(test_seq)
probability = model.predict_proba([features])[0][1]
print(f"抗菌活性概率: {probability:.4f}")
```

## 输出解读

### 模型评估指标

- **准确率 (Accuracy)**: 正确预测的比例
- **ROC-AUC**: 受试者工作特征曲线下面积
- **精确率 (Precision)**: 预测为正例中实际为正例的比例
- **召回率 (Recall)**: 实际正例中被正确预测为正例的比例
- **F1分数**: 精确率和召回率的调和平均

### 进化输出

候选序列文件包含以下列：
- `Rank`: 排名（1-50）
- `Score`: 模型预测得分
- `GenFound`: 发现该序列的代数
- `Sequence`: 肽序列

## 实验设计指南

### 数据集准备建议

1. **数据质量**: 确保序列准确无误，活性测量方法一致
2. **平衡性**: 尽量保持高/低活性样本数量平衡
3. **多样性**: 包含不同长度、来源和作用机制的肽
4. **验证集**: 保留独立验证集用于最终模型评估

### 参数调整策略

1. **保守开始**: 使用默认参数作为基准
2. **逐步调整**: 一次只调整少数参数，评估效果
3. **多轮运行**: 使用不同随机种子运行多次，评估稳定性
4. **领域知识**: 结合肽的生物学知识调整约束条件

### 结果验证

1. **湿实验验证**: 对排名靠前的序列进行实验验证
2. **交叉验证**: 使用不同训练/测试划分验证模型稳健性
3. **外部数据集**: 在独立数据集上测试模型泛化能力
4. **可解释性分析**: 分析特征重要性，理解模型决策依据

## 常见问题解答

### Q: 模型训练需要多长时间？
A: 训练时间取决于数据集大小和硬件配置。对于1000条样本，通常需要1-5分钟。

### Q: 如何调整算法以优化特定性质？
A: 可以修改约束条件和适应度函数。例如，若要增强对某种细菌的活性，可调整训练数据集标签定义。

### Q: 生成的序列长度如何确定？
A: 生成的序列长度与种子序列相同。如需不同长度，可修改进化算法中的序列生成逻辑。

### Q: 模型可以预测哪些性质？
A: 当前模型主要预测抗菌活性（二元分类）。可以通过修改标签数据扩展到其他性质预测。

## 贡献指南

我们欢迎社区贡献！以下是一些参与方式：

1. **问题报告**: 在GitHub上提交问题报告，详细描述问题和复现步骤
2. **功能建议**: 提出新功能或改进建议
3. **代码贡献**: 
   - Fork项目
   - 创建功能分支
   - 提交更改
   - 发起Pull Request
4. **文档改进**: 修正文档错误或改进说明

### 代码贡献规范

1. 遵循PEP 8代码风格
2. 添加适当的注释和文档字符串
3. 确保所有测试通过
4. 更新相关文档

## 版本历史

- **v2.0**: 添加约束进化算法，支持溶解性优化和两亲性引导
- **v1.0**: 基础进化算法实现，支持随机森林分类和特征工程

## 引用

如果您在研究中使用了本项目，请引用：

```
cBD3-DeepEvo-Pept: A Deep Evolution Platform for Antimicrobial Peptide Optimization
[Your Name], [Year]
GitHub Repository: https://github.com/yourusername/cBD3-DeepEvo-Pept
```

## 许可证

本项目采用MIT许可证 - 详见[LICENSE](LICENSE)文件

## 联系方式

- 项目维护者: []
- 邮箱: []
- 项目主页: [https://github.com/yourusername/cBD3-DeepEvo-Pept]

## 致谢

感谢以下开源项目和资源：
- [scikit-learn](https://scikit-learn.org/): 机器学习库
- [Biopython](https://biopython.org/): 生物信息学工具
- [pandas](https://pandas.pydata.org/): 数据处理库
- [numpy](https://numpy.org/): 数值计算库

---

**免责声明**: 本工具仅用于研究目的。生成的肽序列需要经过严格的实验验证才能用于实际应用。
