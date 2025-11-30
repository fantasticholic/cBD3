"""
sequence_evolution_v2.py

专业版序列优化框架（Python）
功能：
- 丰富的特征工程（AAC、dipeptide、Kyte-Doolittle疏水性、hydrophobic moment、helix propensity、Biopython理化量）
- 可复现的训练流程（StratifiedKFold、cross-validation、StandardScaler pipeline）
- 更稳健的分类器（RandomForest + class_weight），支持模型保存/加载
- 更完善的进化策略：种群、锦标赛选择、交叉（crossover）、多点突变（multi-point mutation）、精英保留、种群多样性维护
- 输出候选列表与每代日志

依赖：
- Python 3.8+
- numpy, pandas, scikit-learn, biopython, joblib

用法：直接运行脚本，将 TRAIN_FILE 与 SEED_SEQ 配置为你的数据和种子序列

"""

import math
import random
import itertools
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix, classification_report
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import joblib
import os

# ============== 配置区 ==============
TRAIN_FILE = 'final_train_set.csv'  # CSV should have columns: Sequence, Label (0/1)
MODEL_PATH = 'rf_amp_model.joblib'
SEED_SEQ = "KAWNLRGSAREKAIKNEKLYIFATSGKLAALKPK"
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

# 进化参数
POPULATION_SIZE = 120
GENERATIONS = 80
OFFSPRING_PER_PARENT = 8
ELITE_SIZE = 10
TOURNAMENT_SIZE = 5
MUTATION_RATE = 0.12  # 每个残基发生突变的概率
CROSSOVER_RATE = 0.6
MAX_MUTATION_POINTS = 3

# 氨基酸表（标准 20 种）
AA_LIST = list('ACDEFGHIKLMNPQRSTVWY')
AA_INDEX = {aa:i for i,aa in enumerate(AA_LIST)}

# Kyte-Doolittle 疏水性标量
KD_HYDRO = {
    'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,
    'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5
}

# 近似的 Chou-Fasman alpha helix propensities (常用近似值)
HELIX_PROP = {
    'A':1.45,'R':0.79,'N':0.73,'D':0.98,'C':0.77,'Q':1.17,'E':1.53,'G':0.53,'H':1.00,'I':1.00,
    'L':1.34,'K':1.07,'M':1.20,'F':1.12,'P':0.59,'S':0.79,'T':0.82,'W':1.14,'Y':0.61,'V':1.14
}

# 计算 hydrophobic moment 所用的每残基角度（alpha-helix: 100 deg）
HELICAL_ANGLE_RAD = math.radians(100)

# ============== 特征工程 ==============

def compute_aac(seq):
    """氨基酸组成频率（20 维）"""
    n = len(seq)
    cnt = Counter(seq)
    return np.array([cnt.get(aa,0)/n for aa in AA_LIST], dtype=float)


def compute_dipeptide_freq(seq):
    """二肽频率 (400 维) — 归一化计数"""
    n = len(seq)
    pairs = [seq[i:i+2] for i in range(n-1)]
    cnt = Counter(pairs)
    vec = np.zeros((20*20,), dtype=float)
    for i, a in enumerate(AA_LIST):
        for j, b in enumerate(AA_LIST):
            key = a+b
            vec[i*20 + j] = cnt.get(key, 0) / max(1, len(pairs))
    return vec


def compute_kd_stats(seq):
    """Kyte-Doolittle 相关统计量：均值、方差、亲水/疏水比"""
    vals = [KD_HYDRO.get(aa, 0.0) for aa in seq]
    arr = np.array(vals, dtype=float)
    mean = arr.mean()
    std = arr.std()
    frac_hydrophobic = np.sum(arr > 0) / len(arr)
    return np.array([mean, std, frac_hydrophobic], dtype=float)


def compute_hydrophobic_moment(seq):
    """计算全序列的 hydrophobic moment（基于 alpha 螺旋每残基 100°）
    muH = sqrt((sum(h_i cos(i*theta)))^2 + (sum(h_i sin(i*theta)))^2) / n
    这是一个整体指标，也可按滑窗进一步计算。"""
    n = len(seq)
    hs = [KD_HYDRO.get(aa, 0.0) for aa in seq]
    sum_cos = 0.0
    sum_sin = 0.0
    for i, h in enumerate(hs):
        angle = i * HELICAL_ANGLE_RAD
        sum_cos += h * math.cos(angle)
        sum_sin += h * math.sin(angle)
    mu = math.sqrt(sum_cos**2 + sum_sin**2) / max(1, n)
    return np.array([mu], dtype=float)


def compute_helix_propensity(seq):
    """返回平均 helix propensity 与最大连续 helix propensity 窗口"""
    vals = [HELIX_PROP.get(aa, 0.0) for aa in seq]
    arr = np.array(vals, dtype=float)
    mean = arr.mean()
    # 最大连续窗口的平均（滑窗长度 9，为经验值）
    win = 9
    if len(arr) < win:
        max_avg = arr.mean()
    else:
        max_avg = max(arr[i:i+win].mean() for i in range(len(arr)-win+1))
    return np.array([mean, max_avg], dtype=float)


def compute_basic_biopy(seq):
    """使用 Biopython 提取理化量（分子量、芳香性、instability、pI、净电荷）"""
    pa = ProteinAnalysis(seq)
    mw = pa.molecular_weight()
    arom = pa.aromaticity()
    instab = pa.instability_index()
    pI = pa.isoelectric_point()
    # 净电荷估计 (K+R) - (D+E)
    pos = seq.count('K') + seq.count('R')
    neg = seq.count('D') + seq.count('E')
    net_charge = pos - neg
    return np.array([mw, arom, instab, pI, net_charge], dtype=float)


def get_features(seq):
    """整合特征，返回 1-d numpy array
    特征组成示例： [AAC(20), dipeptide(400), KD_stats(3), muH(1), helix(2), biopy(5), pos_frac(1)]
    总维度约为 432 个特征（根据实现略有差别）。"""
    seq = seq.upper()
    aac = compute_aac(seq)
    dipep = compute_dipeptide_freq(seq)
    kd = compute_kd_stats(seq)
    mu = compute_hydrophobic_moment(seq)
    helix = compute_helix_propensity(seq)
    bio = compute_basic_biopy(seq)
    pos_frac = np.array([(seq.count('K') + seq.count('R')) / len(seq)], dtype=float)
    features = np.concatenate([aac, dipep, kd, mu, helix, bio, pos_frac])
    return features

# ============== 模型训练 ==============

def train_model(train_file=TRAIN_FILE, save_path=MODEL_PATH):
    print('> 载入数据...')
    df = pd.read_csv(train_file)
    assert 'Sequence' in df.columns and 'Label' in df.columns, 'CSV 需包含 Sequence, Label 列'

    X = np.vstack([get_features(s) for s in df['Sequence'].astype(str)])
    y = df['Label'].astype(int).values

    print(f'> 特征维度: {X.shape}')
    print(f'> 类别分布: Label 0: {np.sum(y==0)}, Label 1: {np.sum(y==1)}')

    # pipeline: StandardScaler + RandomForest
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('clf', RandomForestClassifier(n_estimators=200, class_weight='balanced_subsample', random_state=RANDOM_SEED, n_jobs=-1))
    ])

    # 交叉验证评估 - 获取更多详细指标
    print('> 采用 StratifiedKFold 交叉验证评估模型 (5-fold) ...')
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    
    # 收集各折的评估指标
    fold_results = []
    all_y_true = []
    all_y_pred = []
    all_y_pred_proba = []
    
    for fold, (train_idx, val_idx) in enumerate(cv.split(X, y), 1):
        print(f'  Fold {fold}/{cv.get_n_splits()}...')
        
        # 分割数据
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        # 训练模型
        fold_pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', RandomForestClassifier(n_estimators=200, class_weight='balanced_subsample', 
                                         random_state=RANDOM_SEED + fold, n_jobs=-1))
        ])
        fold_pipeline.fit(X_train, y_train)
        
        # 预测
        y_pred = fold_pipeline.predict(X_val)
        y_pred_proba = fold_pipeline.predict_proba(X_val)[:, 1]
        
        # 计算指标
        accuracy = accuracy_score(y_val, y_pred)
        roc_auc = roc_auc_score(y_val, y_pred_proba)
        precision = precision_score(y_val, y_pred)
        recall = recall_score(y_val, y_pred)
        f1 = f1_score(y_val, y_pred)
        
        # 混淆矩阵
        tn, fp, fn, tp = confusion_matrix(y_val, y_pred).ravel()
        
        fold_result = {
            'fold': fold,
            'accuracy': accuracy,
            'roc_auc': roc_auc,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn
        }
        fold_results.append(fold_result)
        
        # 收集所有预测结果
        all_y_true.extend(y_val)
        all_y_pred.extend(y_pred)
        all_y_pred_proba.extend(y_pred_proba)
        
        print(f'    准确率: {accuracy:.3f}, AUC: {roc_auc:.3f}, F1: {f1:.3f}')
    
    # 计算平均指标
    df_results = pd.DataFrame(fold_results)
    avg_metrics = df_results[['accuracy', 'roc_auc', 'precision', 'recall', 'f1_score']].mean()
    std_metrics = df_results[['accuracy', 'roc_auc', 'precision', 'recall', 'f1_score']].std()
    
    # 打印详细的交叉验证结果
    print('\n> 交叉验证综合结果:')
    print(f'  平均准确率 (Accuracy): {avg_metrics["accuracy"]:.3f} ± {std_metrics["accuracy"]:.3f}')
    print(f'  平均AUC (ROC-AUC): {avg_metrics["roc_auc"]:.3f} ± {std_metrics["roc_auc"]:.3f}')
    print(f'  平均精确率 (Precision): {avg_metrics["precision"]:.3f} ± {std_metrics["precision"]:.3f}')
    print(f'  平均召回率 (Recall): {avg_metrics["recall"]:.3f} ± {std_metrics["recall"]:.3f}')
    print(f'  平均F1分数 (F1-Score): {avg_metrics["f1_score"]:.3f} ± {std_metrics["f1_score"]:.3f}')
    
    # 整体混淆矩阵
    overall_cm = confusion_matrix(all_y_true, all_y_pred)
    print('\n> 整体混淆矩阵:')
    print(f'    TP: {overall_cm[1,1]}, FP: {overall_cm[0,1]}')
    print(f'    FN: {overall_cm[1,0]}, TN: {overall_cm[0,0]}')
    
    # 分类报告
    print('\n> 分类报告:')
    report = classification_report(all_y_true, all_y_pred, target_names=['低活性肽 (Label 0)', '高活性肽 (Label 1)'])
    print(report)
    
    # 在全量数据上训练最终模型
    pipeline.fit(X, y)
    joblib.dump(pipeline, save_path)
    print(f'> 模型已保存到: {save_path}')
    return pipeline


def load_model(path=MODEL_PATH):
    if os.path.exists(path):
        print('> 从文件加载模型...')
        return joblib.load(path)
    else:
        return None

# ============== 进化算法辅助函数 ==============

def random_sequence_like(seq):
    return ''.join(random.choice(AA_LIST) for _ in range(len(seq)))


def hamming_distance(a, b):
    assert len(a) == len(b)
    return sum(x!=y for x,y in zip(a,b))


def tournament_selection(population, scores, k=TOURNAMENT_SIZE):
    """简单锦标赛选择，返回一个被选中的序列（index）"""
    idxs = random.sample(range(len(population)), k)
    best = max(idxs, key=lambda i: scores[i])
    return population[best]


def crossover(parent1, parent2):
    """单点或两点交叉（视序列长度）"""
    L = len(parent1)
    if random.random() > CROSSOVER_RATE or L < 4:
        return parent1  # 不交叉，复制父本1
    # 两点交叉
    p1 = random.randint(1, L-3)
    p2 = random.randint(p1+1, L-1)
    child = parent1[:p1] + parent2[p1:p2] + parent1[p2:]
    return child


def mutate(seq, mutation_rate=MUTATION_RATE):
    seq_list = list(seq)
    L = len(seq_list)
    # 决定 mutation points 数量
    for i in range(L):
        if random.random() < mutation_rate:
            # 随机替换，不允许替换成原本残基
            choices = [aa for aa in AA_LIST if aa != seq_list[i]]
            seq_list[i] = random.choice(choices)
    # 额外的多点小概率突变
    if random.random() < 0.02:
        points = random.randint(1, MAX_MUTATION_POINTS)
        for _ in range(points):
            pos = random.randrange(L)
            choices = [aa for aa in AA_LIST if aa != seq_list[pos]]
            seq_list[pos] = random.choice(choices)
    return ''.join(seq_list)

# ============== 进化主循环 ==============

def run_evolution(model, seed_seq=SEED_SEQ, save_candidates='candidates.csv'):
    print('\n> 开始进化搜索...')
    # 初始种群：以 seed + 若干随机序列 + 轻微扰动
    population = [seed_seq]
    # 轻微扰动版本
    for _ in range(int(POPULATION_SIZE*0.3)):
        population.append(mutate(seed_seq, mutation_rate=0.08))
    # 其余随机初始化
    while len(population) < POPULATION_SIZE:
        population.append(random_sequence_like(seed_seq))

    # 记录最优
    best_history = []  # tuples of (generation, score, seq)

    # 计算特征批量化函数
    def batch_score(seqs):
        feats = np.vstack([get_features(s) for s in seqs])
        probs = model.predict_proba(feats)[:,1]
        return probs

    # 评估初始种群
    pop_scores = batch_score(population)

    for gen in range(1, GENERATIONS+1):
        # 保留精英
        elite_idx = np.argsort(pop_scores)[-ELITE_SIZE:]
        elites = [population[i] for i in elite_idx]
        elites_scores = pop_scores[elite_idx]

        # 生成新一代
        new_population = elites.copy()
        # 生出 offspring
        while len(new_population) < POPULATION_SIZE:
            # 选择两个父本
            p1 = tournament_selection(population, pop_scores)
            p2 = tournament_selection(population, pop_scores)
            child = crossover(p1, p2)
            child = mutate(child)
            new_population.append(child)

        # 引入少量随机移民以维持多样性
        n_immigrants = max(1, int(0.03 * POPULATION_SIZE))
        for i in range(n_immigrants):
            new_population[-1-i] = random_sequence_like(seed_seq)

        # 评估
        new_scores = batch_score(new_population)

        # 多样性检查（基于平均 Hamming 距离）
        # 如果种群过于相似，强制注入若干随机序列
        avg_hd = np.mean([hamming_distance(a,b) for a,b in zip(new_population[:-1], new_population[1:])])
        if avg_hd < max(1, len(seed_seq)*0.15):
            for i in range(3):
                idx = random.randrange(ELITE_SIZE, POPULATION_SIZE)
                new_population[idx] = random_sequence_like(seed_seq)
                new_scores[idx] = 0.5  # 中性分数

        # 更新
        population = new_population
        pop_scores = new_scores

        # 记录并打印当代最好
        best_idx = int(np.argmax(pop_scores))
        best_score = float(pop_scores[best_idx])
        best_seq = population[best_idx]
        best_history.append((gen, best_score, best_seq))

        if gen % 5 == 0 or gen == 1 or gen == GENERATIONS:
            print(f'  Gen {gen:03d} | BestScore: {best_score:.4f} | BestSeq: {best_seq} | AvgHD: {avg_hd:.2f}')

    # 汇总 top candidates
    unique = {}
    for gen, score, seq in best_history:
        if seq not in unique or score > unique[seq][0]:
            unique[seq] = (score, gen)
    sorted_cands = sorted([(v[0], v[1], s) for s,v in unique.items()], key=lambda x: x[0], reverse=True)

    # 保存结果
    out_rows = []
    for rank, (score, gen, seq) in enumerate(sorted_cands[:50], start=1):
        out_rows.append({'Rank': rank, 'Score': score, 'GenFound': gen, 'Sequence': seq})
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(save_candidates, index=False)
    print(f'> Top {len(out_rows)} 候选已保存到: {save_candidates}')
    print('> 结束进化搜索。')
    return out_df

# ============== 主流程 ==============

def main():
    model = load_model()
    if model is None:
        model = train_model()
    # 快速做一次训练集上的表现检查（可选）
    df = pd.read_csv(TRAIN_FILE)
    X = np.vstack([get_features(s) for s in df['Sequence'].astype(str)])
    y = df['Label'].astype(int).values
    preds = model.predict(X)
    probs = model.predict_proba(X)[:, 1]
    
    print(f"\n> 训练集整体表现:")
    print(f"  准确率: {accuracy_score(y, preds):.3f}")
    print(f"  AUC: {roc_auc_score(y, probs):.3f}")
    print(f"  F1分数: {f1_score(y, preds):.3f}")

    # 运行进化
    print("\n" + "="*50)
    out = run_evolution(model)
    print("\n" + "="*50)
    print("\n> Top 10 候选序列:")
    print(out.head(10).to_string(index=False))

if __name__ == '__main__':
    main()
