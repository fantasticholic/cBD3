import pandas as pd
import re
import os

def clean_and_label_amp_data(input_file, output_file):
    """
    清洗抗菌肽数据并进行二分类标签化
    
    Args:
        input_file: 输入CSV文件路径
        output_file: 输出CSV文件路径
    """
    # 读取CSV文件
    print(f"正在读取CSV文件: {input_file}")
    try:
        df = pd.read_csv(input_file)
        print(f"初始数据行数: {len(df)}")
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return
    
    # 数据预处理函数
    def extract_numeric_value(value):
        """从字符串中提取数值，如果无法提取则返回NaN"""
        if pd.isna(value):
            return float('nan')
        
        # 转换为字符串并去除空白
        value_str = str(value).strip()
        
        # 处理包含'>'或'<'等符号的情况
        if any(char in value_str for char in ['>', '<', '=', '≥', '≤']):
            # 提取数字部分
            numbers = re.findall(r'\d+(?:\.\d+)?', value_str)
            if numbers:
                return float(numbers[0])
            else:
                return float('nan')
        
        # 尝试直接转换为浮点数
        try:
            return float(value_str)
        except:
            return float('nan')
    
    # 预处理MIC列（注意列名可能有空格）
    mic_column = None
    for col in df.columns:
        if 'MIC' in col:
            mic_column = col
            break
    
    hemolysis_column = None
    for col in df.columns:
        if 'Hemolysis' in col and '%' in col:
            hemolysis_column = col
            break
    
    # 如果找不到特定列，使用用户提到的默认列名
    if not mic_column:
        mic_column = 'MIC (µg/M)'
    if not hemolysis_column:
        hemolysis_column = 'Hemolysis (%)'
    
    print(f"使用的MIC列: {mic_column}")
    print(f"使用的Hemolysis列: {hemolysis_column}")
    
    # 将MIC和Hemolysis列转换为数值类型
    print("正在预处理数据...")
    if mic_column in df.columns:
        df[mic_column] = df[mic_column].apply(extract_numeric_value)
    if hemolysis_column in df.columns:
        df[hemolysis_column] = df[hemolysis_column].apply(extract_numeric_value)
    
    # 删除Sequence或MIC为空的行
    initial_rows = len(df)
    df = df.dropna(subset=['Sequence', mic_column])
    removed_rows = initial_rows - len(df)
    print(f"删除空值行: {removed_rows} 行")
    print(f"处理后数据行数: {len(df)}")
    
    # 执行标签化逻辑
    print("正在进行标签化...")
    # 创建空的Label列
    df['Label'] = None
    
    # 好肽（Label = 1）：MIC ≤ 10 且 Hemolysis < 20
    # 检查hemolysis_column是否存在，不存在则只使用MIC条件
    if hemolysis_column in df.columns:
        df.loc[(df[mic_column] <= 10) & (df[hemolysis_column] < 20), 'Label'] = 1
    else:
        df.loc[df[mic_column] <= 10, 'Label'] = 1
    
    # 差肽（Label = 0）：MIC ≥ 32
    df.loc[df[mic_column] >= 32, 'Label'] = 0
    
    # 删除中间地带（Label为空的行）
    before_filter = len(df)
    df = df.dropna(subset=['Label'])
    df['Label'] = df['Label'].astype(int)  # 确保Label是整数类型
    removed_middle = before_filter - len(df)
    print(f"删除中间地带数据: {removed_middle} 行")
    
    # 添加锚点数据（先导肽）
    anchor_sequence = "KAWNLRGSAREKAIKNEKLYIFATSGKLAALKPK"
    anchor_mic = 64
    anchor_hemolysis = 5
    anchor_label = 0
    
    print("\n添加锚点数据...")
    # 检查数据集中是否已存在该序列
    existing_anchor = df[df['Sequence'] == anchor_sequence]
    
    if not existing_anchor.empty:
        # 如果存在，确保其Label为0
        print("锚点序列已存在，更新其标签和属性")
        df.loc[df['Sequence'] == anchor_sequence, 'Label'] = anchor_label
        if mic_column in df.columns:
            df.loc[df['Sequence'] == anchor_sequence, mic_column] = anchor_mic
        if hemolysis_column in df.columns:
            df.loc[df['Sequence'] == anchor_sequence, hemolysis_column] = anchor_hemolysis
    else:
        # 如果不存在，添加新行
        print("锚点序列不存在，添加新行")
        anchor_data = {'Sequence': anchor_sequence, 'Label': anchor_label}
        if mic_column in df.columns:
            anchor_data[mic_column] = anchor_mic
        if hemolysis_column in df.columns:
            anchor_data[hemolysis_column] = anchor_hemolysis
        
        # 确保新行包含所有必要的列
        for col in df.columns:
            if col not in anchor_data:
                anchor_data[col] = None
        
        # 添加到DataFrame
        df = pd.concat([df, pd.DataFrame([anchor_data])], ignore_index=True)
    
    # 打印统计信息
    print("\n数据统计信息:")
    print(f"清洗后总行数: {len(df)}")
    print(f"正样本数量 (Label=1): {(df['Label'] == 1).sum()}")
    print(f"负样本数量 (Label=0): {(df['Label'] == 0).sum()}")
    
    # 保存结果
    print(f"\n正在保存结果到: {output_file}")
    try:
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"保存成功！文件大小: {os.path.getsize(output_file) / 1024:.2f} KB")
    except Exception as e:
        print(f"保存文件时出错: {e}")
    
    # 返回处理后的数据（用于调试）
    return df

if __name__ == "__main__":
    # 定义文件路径
    input_csv = "d:/trea/trae_projects/cBD3-ABU/my_500_peptides.csv"
    output_csv = "d:/trea/trae_projects/cBD3-ABU/final_train_set.csv"
    
    # 运行处理函数
    clean_and_label_amp_data(input_csv, output_csv)