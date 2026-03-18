import pandas as pd
import os

def validate_amp_dataset(dataset_file):
    """
    验证抗菌肽训练数据集的正确性
    
    Args:
        dataset_file: 数据集文件路径
    """
    # 检查文件是否存在
    if not os.path.exists(dataset_file):
        print(f"错误：找不到数据集文件 {dataset_file}")
        return False
    
    print(f"正在验证数据集：{dataset_file}")
    
    try:
        # 读取数据集
        df = pd.read_csv(dataset_file)
        
        print(f"\n基本信息：")
        print(f"数据集行数: {len(df)}")
        print(f"数据集列数: {len(df.columns)}")
        print(f"\n列名列表：")
        for i, col in enumerate(df.columns, 1):
            print(f"{i}. {col}")
        
        # 检查必要的列是否存在
        required_columns = ['Sequence', 'Label']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"\n错误：缺少必要的列: {missing_columns}")
            return False
        
        # 检查标签分布
        label_counts = df['Label'].value_counts()
        print(f"\n标签分布：")
        for label, count in sorted(label_counts.items()):
            percentage = (count / len(df)) * 100
            print(f"  Label {label}: {count} 样本 ({percentage:.1f}%)")
        
        # 检查是否有NaN值
        nan_summary = df.isna().sum()
        columns_with_nan = nan_summary[nan_summary > 0]
        if not columns_with_nan.empty:
            print(f"\n存在NaN值的列：")
            for col, count in columns_with_nan.items():
                print(f"  {col}: {count} 个NaN值")
        else:
            print("\n没有发现NaN值")
        
        # 检查锚点数据是否存在
        anchor_sequence = "KAWNLRGSAREKAIKNEKLYIFATSGKLAALKPK"
        anchor_data = df[df['Sequence'] == anchor_sequence]
        
        if anchor_data.empty:
            print(f"\n错误：未找到锚点数据")
            return False
        else:
            print(f"\n找到锚点数据：")
            print(f"  序列: {anchor_sequence}")
            print(f"  标签: {anchor_data['Label'].values[0]}")
            if 'MIC (µg/M)' in df.columns:
                print(f"  MIC: {anchor_data['MIC (µg/M)'].values[0]}")
            if 'Hemolysis (%)' in df.columns:
                print(f"  Hemolysis: {anchor_data['Hemolysis (%)'].values[0]}%")
        
        # 验证标签逻辑是否正确
        print(f"\n验证标签逻辑：")
        if 'MIC (µg/M)' in df.columns and 'Hemolysis (%)' in df.columns:
            # 检查正样本是否符合条件：MIC ≤ 10 且 Hemolysis < 20
            positive_samples = df[df['Label'] == 1]
            invalid_positives = positive_samples[~((positive_samples['MIC (µg/M)'] <= 10) & 
                                                  (positive_samples['Hemolysis (%)'] < 20))]
            
            if len(invalid_positives) > 0:
                print(f"  警告：发现 {len(invalid_positives)} 个不符合条件的正样本")
            else:
                print(f"  所有正样本都符合条件 (MIC ≤ 10 且 Hemolysis < 20)")
            
            # 检查负样本是否符合条件：MIC ≥ 32
            negative_samples = df[df['Label'] == 0]
            invalid_negatives = negative_samples[negative_samples['MIC (µg/M)'] < 32]
            
            # 排除锚点数据的检查，因为它是强制添加的
            invalid_negatives = invalid_negatives[invalid_negatives['Sequence'] != anchor_sequence]
            
            if len(invalid_negatives) > 0:
                print(f"  警告：发现 {len(invalid_negatives)} 个不符合条件的负样本")
            else:
                print(f"  所有负样本都符合条件 (MIC ≥ 32)，排除锚点数据")
        
        # 显示前5行数据预览
        print(f"\n前5行数据预览：")
        preview_columns = ['Sequence', 'Label']
        if 'MIC (µg/M)' in df.columns:
            preview_columns.append('MIC (µg/M)')
        if 'Hemolysis (%)' in df.columns:
            preview_columns.append('Hemolysis (%)')
        
        print(df[preview_columns].head())
        
        print(f"\n数据集验证完成！")
        return True
        
    except Exception as e:
        print(f"验证过程中出现错误：{str(e)}")
        return False

if __name__ == "__main__":
    # 定义数据集文件路径
    dataset_file = "d:/trea/trae_projects/cBD3-ABU/final_train_set.csv"
    
    # 运行验证
    validate_amp_dataset(dataset_file)