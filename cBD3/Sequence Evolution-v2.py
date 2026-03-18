# professional_evolution.py
"""
更严格的多肽进化设计脚本（带硬约束、疏水惩罚与两亲性交替引导）
- 绝对禁止生成 C (Cysteine)
- 连续 >=3 强疏水残基 (L, V, I, F, W) 的序列会被直接惩罚为分数 0
- 总疏水比例 > 50% 会被按比例扣分
- 突变函数优先引导交替排列（亲水-疏水-亲水...）
依赖:
  numpy, pandas, scikit-learn, biopython, joblib
"""

import math
import random
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import joblib
import os

# -------------------- 配置 --------------------
TRAIN_FILE = 'final_train_set.csv'   # 需包含 Sequence, Label
MODEL_PATH = 'rf_amp_model_constrained.joblib'
SEED_SEQ = "KAWNLRGSAREKAIKNEKLYIFATSGKLAALKPK"  # 默认种子（不含 C）
RANDOM_SEED = 42

np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

POPULATION_SIZE = 120
GENERATIONS = 60
ELITE_SIZE = 10
TOURNAMENT_SIZE = 4
CROSSOVER_RATE = 0.6
MUTATION_RATE = 0.12

# 强疏水残基集合 (会触发连续惩罚)
HYDROPHOBIC_STRONG = set(['L', 'V', 'I', 'F', 'W'])
# 允许用于生成/突变的氨基酸（**移除 'C'**）
ALLOWED_AA = ['A','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
AA_LIST = ALLOWED_AA.copy()

# Kyte-Doolittle 简单疏水性标量（用于总体疏水比例计算）
KD_HYDRO = {
    'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,
    'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5
}
# 亲水/带电残基集合（引导交替）
POLAR_CHARGED = set(['K','R','D','E','Q','N','S','T','H','Y'])

# -------------------- 特征/工具函数 --------------------
def compute_aac(seq):
    n = len(seq)
    cnt = Counter(seq)
    return np.array([cnt.get(aa,0)/n for aa in AA_LIST], dtype=float)

def compute_basic_biopy(seq):
    pa = ProteinAnalysis(seq)
    mw = pa.molecular_weight()
    arom = pa.aromaticity()
    instab = pa.instability_index()
    pI = pa.isoelectric_point()
    pos = seq.count('K') + seq.count('R')
    neg = seq.count('D') + seq.count('E')
    net_charge = pos - neg
    return np.array([mw, arom, instab, pI, net_charge], dtype=float)

def get_features(seq):
    # 简化版本：AAC + Biopython 基本理化量（训练评分目的）
    seq = seq.upper()
    aac = compute_aac(seq)
    bio = compute_basic_biopy(seq)
    return np.concatenate([aac, bio])

# -------------------- 模型训练/加载 --------------------
def train_model(train_file=TRAIN_FILE, save_path=MODEL_PATH):
    if not os.path.exists(train_file):
        raise FileNotFoundError(f"训练文件未找到: {train_file}")
    df = pd.read_csv(train_file)
    assert 'Sequence' in df.columns and 'Label' in df.columns, "CSV 需包含 Sequence 和 Label 列"

    X = np.vstack([get_features(s) for s in df['Sequence'].astype(str)])
    y = df['Label'].astype(int).values

    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('clf', RandomForestClassifier(n_estimators=150, class_weight='balanced_subsample', random_state=RANDOM_SEED, n_jobs=-1))
    ])

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    aucs = cross_val_score(pipeline, X, y, cv=cv, scoring='roc_auc', n_jobs=-1)
    print(f"[train] CV ROC-AUC: {aucs.mean():.3f} ± {aucs.std():.3f}")

    pipeline.fit(X, y)
    joblib.dump(pipeline, save_path)
    print(f"[train] 模型保存到 {save_path}")
    return pipeline

def load_model(path=MODEL_PATH):
    if os.path.exists(path):
        return joblib.load(path)
    return None

# -------------------- Solubility / Safety 检查 --------------------
def has_consecutive_strong_hydro(seq, threshold=3):
    """若存在连续 threshold 个或以上的强疏水残基 -> True"""
    run = 0
    for aa in seq:
        if aa in HYDROPHOBIC_STRONG:
            run += 1
            if run >= threshold:
                return True
        else:
            run = 0
    return False

def hydrophobic_fraction(seq):
    """返回强疏水残基占比（L,V,I,F,W）"""
    if len(seq) == 0:
        return 0.0
    count = sum(1 for aa in seq if aa in HYDROPHOBIC_STRONG)
    return count / len(seq)

def check_solubility_penalty(seq):
    """
    返回一个 penalty_multiplier (0.0 - 1.0)
    - 若存在连续 >=3 强疏水残基：直接返回 0.0（Rule A: 严厉惩罚）
    - 否则若总疏水比例 > 0.5：按比例扣分，线性从 1 -> 0.2 (在 0.5 -> 1.0 区间)
    - 否则返回 1.0（无惩罚）
    """
    seq = seq.upper()
    if has_consecutive_strong_hydro(seq, threshold=3):
        return 0.0
    frac = hydrophobic_fraction(seq)
    if frac > 0.5:
        # frac==0.5 -> multiplier 1.0; frac==1.0 -> multiplier 0.2
        mul = max(0.2, 1.0 - (frac - 0.5) * 1.6)  # 线性压缩
        return mul
    return 1.0

# -------------------- 生成与突变（硬约束） --------------------
def random_sequence_like(seq_length):
    """生成不含 C 的随机序列"""
    return ''.join(random.choice(AA_LIST) for _ in range(seq_length))

def guided_mutation_single_residue(seq, pos):
    """
    针对位置 pos 给出一个更“带两亲性倾向”的替换：
    - 若当前位置或邻居导致连续疏水风险，则优先选择亲水/带电残基 (K,R,Q,S,T,D,E,N,H)
    - 若附近为亲水，允许或鼓励选择疏水以维持交替（但不形成连续 >=3）
    - 严格禁止选择 'C'
    """
    L = len(seq)
    left = seq[pos-1] if pos-1 >= 0 else None
    right = seq[pos+1] if pos+1 < L else None

    # 统计邻近是否为强疏水
    neighbor_hydro = sum(1 for x in (left, right) if x in HYDROPHOBIC_STRONG)
    # 构建候选集（禁止 C）
    # 如果周边有 >=1 强疏水 -> 倾向亲水/带电
    polar_candidates = [aa for aa in ALLOWED_AA if aa in POLAR_CHARGED]
    hydro_candidates = [aa for aa in ALLOWED_AA if aa in HYDROPHOBIC_STRONG]
    neutral_candidates = [aa for aa in ALLOWED_AA if aa not in POLAR_CHARGED and aa not in HYDROPHOBIC_STRONG]

    # 优先策略
    if neighbor_hydro >= 1:
        # 强制倾向亲水/带电（但仍保留少量概率选中中性/疏水）
        choices = polar_candidates * 6 + neutral_candidates * 2 + hydro_candidates * 1
    else:
        # 周围不多疏水，允许交替：优先选择疏水如果左右为亲水
        neighbor_polar = sum(1 for x in (left, right) if x in POLAR_CHARGED)
        if neighbor_polar >= 1:
            choices = hydro_candidates * 6 + polar_candidates * 2 + neutral_candidates * 1
        else:
            # 周围中性，保持多样性
            choices = ALLOWED_AA.copy()
    # 绝对确保不返回原残基（若可能）
    orig = seq[pos]
    filtered = [c for c in choices if c != orig]
    if not filtered:
        return orig
    return random.choice(filtered)

def mutate(seq, mutation_rate=MUTATION_RATE):
    """逐位突变，但受硬约束（无 C）且带引导策略以鼓励交替"""
    seq_list = list(seq)
    L = len(seq_list)

    for i in range(L):
        if random.random() < mutation_rate:
            # 使用 guided_mutation_single_residue 以鼓励两亲性交替
            new_aa = guided_mutation_single_residue(seq_list, i)
            seq_list[i] = new_aa

    # 小概率进行额外随机替换（仍禁止 C）
    if random.random() < 0.03:
        pos = random.randrange(L)
        choices = [aa for aa in ALLOWED_AA if aa != seq_list[pos]]
        seq_list[pos] = random.choice(choices)

    # 最后保证不会产生连续 >=3 强疏水（如发生则在末尾替换一个为亲水）
    # 这里仅作最后保险（理论上 guided_mutation 使其概率极低）
    if has_consecutive_strong_hydro(''.join(seq_list), threshold=3):
        # 找到第一个触发位置并修正其中一个位为亲水/带电
        s = ''.join(seq_list)
        run = 0
        for idx, aa in enumerate(s):
            if aa in HYDROPHOBIC_STRONG:
                run += 1
                if run >= 3:
                    # 将当前或前一个位置替换为亲水/带电残基
                    replace_pos = idx if random.random() < 0.6 else idx-1
                    seq_list[replace_pos] = random.choice([a for a in ALLOWED_AA if a in POLAR_CHARGED])
                    break
            else:
                run = 0

    return ''.join(seq_list)

# -------------------- 选择/交叉 --------------------
def tournament_selection(population, scores, k=TOURNAMENT_SIZE):
    idxs = random.sample(range(len(population)), k)
    best_idx = max(idxs, key=lambda i: scores[i])
    return population[best_idx]

def crossover(parent1, parent2):
    L = len(parent1)
    if random.random() > CROSSOVER_RATE or L < 4:
        return parent1
    p1 = random.randint(1, L-3)
    p2 = random.randint(p1+1, L-1)
    child = parent1[:p1] + parent2[p1:p2] + parent1[p2:]
    # 交叉可能引入 C（若父本含 C），但我们确保父本不会含 C（初始化和突变保证）
    return child

# -------------------- 进化主循环（带疏水惩罚） --------------------
def run_evolution(model, seed_seq=SEED_SEQ, save_candidates='candidates_constrained.csv'):
    print("[evolution] Start evolution with hard constraints (no C) and solubility filter")
    seq_len = len(seed_seq)

    # 初始化种群：以 seed + 轻微扰动 + 随机（均不含 C）
    population = [seed_seq]
    for _ in range(int(POPULATION_SIZE * 0.25)):
        population.append(mutate(seed_seq))
    while len(population) < POPULATION_SIZE:
        population.append(random_sequence_like(seq_len))

    def batch_score(seqs):
        feats = np.vstack([get_features(s) for s in seqs])
        probs = model.predict_proba(feats)[:,1]
        # 应用 solubility penalty：连续疏水直接 0；总疏水比例高则乘以系数
        adjusted = []
        for seq, p in zip(seqs, probs):
            mul = check_solubility_penalty(seq)
            adjusted.append(p * mul)
        return np.array(adjusted, dtype=float)

    # 初次评估
    pop_scores = batch_score(population)

    best_history = []
    for gen in range(1, GENERATIONS + 1):
        # 精英保留
        elite_idx = np.argsort(pop_scores)[-ELITE_SIZE:]
        elites = [population[i] for i in elite_idx]
        elites_scores = pop_scores[elite_idx]

        # 新一代先放入 elites
        new_pop = elites.copy()

        # 产生 offspring
        while len(new_pop) < POPULATION_SIZE:
            p1 = tournament_selection(population, pop_scores)
            p2 = tournament_selection(population, pop_scores)
            child = crossover(p1, p2)
            child = mutate(child)
            # 确保 child 不含 C（保险）
            if 'C' in child:
                child = child.replace('C', random.choice([aa for aa in ALLOWED_AA if aa != 'C']))
            new_pop.append(child)

        # 注入少量随机移民以维持多样性（也不含 C）
        n_imm = max(1, int(0.03 * POPULATION_SIZE))
        for i in range(n_imm):
            idx = random.randrange(ELITE_SIZE, POPULATION_SIZE)
            new_pop[idx] = random_sequence_like(seq_len)

        # 评估新一代
        new_scores = batch_score(new_pop)

        # 若平均 Hamming 距离过小 -> 用更多移民
        avg_hd = np.mean([sum(a!=b for a,b in zip(new_pop[i], new_pop[(i+1)%len(new_pop)])) for i in range(len(new_pop))])
        if avg_hd < max(1, seq_len * 0.12):
            # 插入更多随机个体
            for i in range(3):
                idx = random.randrange(ELITE_SIZE, POPULATION_SIZE)
                new_pop[idx] = random_sequence_like(seq_len)
                new_scores[idx] = 0.5

        # 更新
        population = new_pop
        pop_scores = new_scores

        # 记录当代最优
        best_idx = int(np.argmax(pop_scores))
        best_score = float(pop_scores[best_idx])
        best_seq = population[best_idx]
        best_history.append((gen, best_score, best_seq))

        if gen % 5 == 0 or gen == 1 or gen == GENERATIONS:
            print(f"[gen {gen:03d}] BestScore={best_score:.4f} | BestSeq={best_seq} | AvgHD={avg_hd:.2f} | HydFrac={hydrophobic_fraction(best_seq):.3f}")

    # 汇总 top candidates（去重并按 score 排序）
    unique = {}
    for gen, score, seq in best_history:
        if seq not in unique or score > unique[seq][0]:
            unique[seq] = (score, gen)
    sorted_cands = sorted([(v[0], v[1], s) for s,v in unique.items()], key=lambda x: x[0], reverse=True)

    # 保存 top50
    out_rows = []
    for rank, (score, gen, seq) in enumerate(sorted_cands[:50], start=1):
        out_rows.append({'Rank': rank, 'Score': score, 'GenFound': gen, 'Sequence': seq})
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(save_candidates, index=False)
    print(f"[evolution] Saved top candidates to {save_candidates}")
    return out_df

# -------------------- 主流程 --------------------
def main():
    # 加载或训练模型
    model = load_model()
    if model is None:
        print("[main] 模型不存在，开始训练（可能耗时）。")
        model = train_model()
    else:
        print("[main] 已加载模型。")

    # （可选）在训练集上检查表现
    if os.path.exists(TRAIN_FILE):
        df = pd.read_csv(TRAIN_FILE)
        X = np.vstack([get_features(s) for s in df['Sequence'].astype(str)])
        y = df['Label'].astype(int).values
        preds = model.predict(X)
        acc = (preds == y).mean()
        print(f"[main] Training set accuracy (full): {acc:.3f}")

    # 运行受限的进化搜索
    out = run_evolution(model, seed_seq=SEED_SEQ, save_candidates='candidates_constrained.csv')
    print(out.head(10).to_string(index=False))

if __name__ == '__main__':
    main()
