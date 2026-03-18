import pandas as pd
import re
import os
import argparse
from collections import Counter

class AMPDataProcessor:
    """
    抗菌肽数据处理器，提供数据清洗、标记、验证和分析功能
    """
    
    def __init__(self):
        # 定义关键列名和标签阈值
        self.mic_threshold_high = 10  # MIC ≤ 10 为好肽
        self.mic_threshold_low = 32   # MIC ≥ 32 为差肽
        self.hemolysis_threshold = 20  # Hemolysis < 20 为好肽
        self.anchor_sequence = "KAWNLRGSAREKAIKNEKLYIFATSGKLAALKPK"
        self.anchor_mic = 64
        self.anchor_hemolysis = 5
        self.anchor_label = 0
    
    def extract_numeric_value(self, value):
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
    
    def find_column(self, df, keywords):
        """根据关键词列表查找列"""
        for col in df.columns:
            for keyword in keywords:
                if keyword in col:
                    return col
        return None
    
    def clean_and_label_data(self, input_file, output_file=None):
        """
        清洗抗菌肽数据并进行二分类标签化
        
        Args:
            input_file: 输入CSV文件路径
            output_file: 输出CSV文件路径，如不提供则返回处理后的数据
        
        Returns:
            处理后的DataFrame或None（如果出错）
        """
        # 读取CSV文件
        print(f"正在读取CSV文件: {input_file}")
        try:
            df = pd.read_csv(input_file)
            print(f"初始数据行数: {len(df)}")
        except Exception as e:
            print(f"读取文件时出错: {e}")
            return None
        
        # 查找关键列
        mic_column = self.find_column(df, ['MIC', 'mic'])
        hemolysis_column = self.find_column(df, ['Hemolysis', 'hemolysis'])
        sequence_column = self.find_column(df, ['Sequence', 'sequence'])
        
        # 如果找不到特定列，使用默认列名
        if not sequence_column:
            sequence_column = 'Sequence'
        if not mic_column:
            mic_column = 'MIC (µg/M)'
        if not hemolysis_column:
            hemolysis_column = 'Hemolysis (%)'
        
        print(f"使用的序列列: {sequence_column}")
        print(f"使用的MIC列: {mic_column}")
        print(f"使用的Hemolysis列: {hemolysis_column}")
        
        # 将MIC和Hemolysis列转换为数值类型
        print("正在预处理数据...")
        if mic_column in df.columns:
            df[mic_column] = df[mic_column].apply(self.extract_numeric_value)
        if hemolysis_column in df.columns:
            df[hemolysis_column] = df[hemolysis_column].apply(self.extract_numeric_value)
        
        # 删除序列或MIC为空的行
        initial_rows = len(df)
        required_columns = [sequence_column]
        if mic_column in df.columns:
            required_columns.append(mic_column)
            
        df = df.dropna(subset=required_columns)
        removed_rows = initial_rows - len(df)
        print(f"删除空值行: {removed_rows} 行")
        print(f"处理后数据行数: {len(df)}")
        
        # 执行标签化逻辑
        print("正在进行标签化...")
        # 创建空的Label列
        df['Label'] = None
        
        # 好肽（Label = 1）：MIC ≤ 10 且 Hemolysis < 20
        # 差肽（Label = 0）：MIC ≥ 32
        if mic_column in df.columns:
            if hemolysis_column in df.columns:
                df.loc[(df[mic_column] <= self.mic_threshold_high) & 
                      (df[hemolysis_column] < self.hemolysis_threshold), 'Label'] = 1
            else:
                df.loc[df[mic_column] <= self.mic_threshold_high, 'Label'] = 1
            
            df.loc[df[mic_column] >= self.mic_threshold_low, 'Label'] = 0
        
        # 删除中间地带（Label为空的行）
        before_filter = len(df)
        df = df.dropna(subset=['Label'])
        df['Label'] = df['Label'].astype(int)  # 确保Label是整数类型
        removed_middle = before_filter - len(df)
        print(f"删除中间地带数据: {removed_middle} 行")
        
        # 添加锚点数据（先导肽）
        print("\n添加锚点数据...")
        # 检查数据集中是否已存在该序列
        existing_anchor = df[df[sequence_column] == self.anchor_sequence]
        
        if not existing_anchor.empty:
            # 如果存在，确保其Label为0
            print("锚点序列已存在，更新其标签和属性")
            df.loc[df[sequence_column] == self.anchor_sequence, 'Label'] = self.anchor_label
            if mic_column in df.columns:
                df.loc[df[sequence_column] == self.anchor_sequence, mic_column] = self.anchor_mic
            if hemolysis_column in df.columns:
                df.loc[df[sequence_column] == self.anchor_sequence, hemolysis_column] = self.anchor_hemolysis
        else:
            # 如果不存在，添加新行
            print("锚点序列不存在，添加新行")
            anchor_data = {sequence_column: self.anchor_sequence, 'Label': self.anchor_label}
            if mic_column in df.columns:
                anchor_data[mic_column] = self.anchor_mic
            if hemolysis_column in df.columns:
                anchor_data[hemolysis_column] = self.anchor_hemolysis
            
            # 确保新行包含所有必要的列
            for col in df.columns:
                if col not in anchor_data:
                    anchor_data[col] = None
            
            # 添加到DataFrame
            df = pd.concat([df, pd.DataFrame([anchor_data])], ignore_index=True)
        
        # 数据质量分析
        print("\n序列质量分析:")
        # 检查序列长度分布
        df['Sequence_Length'] = df[sequence_column].str.len()
        print(f"序列长度范围: {df['Sequence_Length'].min()} - {df['Sequence_Length'].max()} 氨基酸")
        print(f"平均序列长度: {df['Sequence_Length'].mean():.1f} 氨基酸")
        
        # 统计氨基酸组成（前5行示例）
        print("\n氨基酸组成分析 (前50个序列):")
        sample_sequences = df[sequence_column].head(50)
        amino_acid_counts = Counter(''.join(sample_sequences))
        total_amino_acids = sum(amino_acid_counts.values())
        for aa, count in sorted(amino_acid_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
            percentage = (count / total_amino_acids) * 100
            print(f"  {aa}: {count} ({percentage:.1f}%)")
        
        # 打印统计信息
        print("\n数据统计信息:")
        print(f"清洗后总行数: {len(df)}")
        print(f"正样本数量 (Label=1): {(df['Label'] == 1).sum()}")
        print(f"负样本数量 (Label=0): {(df['Label'] == 0).sum()}")
        
        # 保存结果
        if output_file:
            print(f"\n正在保存结果到: {output_file}")
            try:
                df.to_csv(output_file, index=False, encoding='utf-8-sig')
                print(f"保存成功！文件大小: {os.path.getsize(output_file) / 1024:.2f} KB")
            except Exception as e:
                print(f"保存文件时出错: {e}")
        
        # 返回处理后的数据
        return df
    
    def validate_dataset(self, dataset_file):
        """
        验证抗菌肽训练数据集的正确性
        
        Args:
            dataset_file: 数据集文件路径
            
        Returns:
            验证是否通过 (True/False)
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
            
            # 查找关键列
            sequence_column = self.find_column(df, ['Sequence', 'sequence'])
            mic_column = self.find_column(df, ['MIC', 'mic'])
            hemolysis_column = self.find_column(df, ['Hemolysis', 'hemolysis'])
            
            # 检查必要的列是否存在
            required_columns = ['Label']
            if sequence_column:
                required_columns.append(sequence_column)
            
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
            if sequence_column:
                anchor_data = df[df[sequence_column] == self.anchor_sequence]
                
                if anchor_data.empty:
                    print(f"\n错误：未找到锚点数据")
                    return False
                else:
                    print(f"\n找到锚点数据：")
                    print(f"  序列: {self.anchor_sequence}")
                    print(f"  标签: {anchor_data['Label'].values[0]}")
                    if mic_column and mic_column in df.columns:
                        print(f"  MIC: {anchor_data[mic_column].values[0]}")
                    if hemolysis_column and hemolysis_column in df.columns:
                        print(f"  Hemolysis: {anchor_data[hemolysis_column].values[0]}%")
            
            # 验证标签逻辑是否正确
            print(f"\n验证标签逻辑：")
            if mic_column and mic_column in df.columns and hemolysis_column and hemolysis_column in df.columns:
                # 检查正样本是否符合条件：MIC ≤ 10 且 Hemolysis < 20
                positive_samples = df[df['Label'] == 1]
                invalid_positives = positive_samples[~((positive_samples[mic_column] <= self.mic_threshold_high) & 
                                                      (positive_samples[hemolysis_column] < self.hemolysis_threshold))]
                
                if len(invalid_positives) > 0:
                    print(f"  警告：发现 {len(invalid_positives)} 个不符合条件的正样本")
                else:
                    print(f"  所有正样本都符合条件 (MIC ≤ {self.mic_threshold_high} 且 Hemolysis < {self.hemolysis_threshold})")
                
                # 检查负样本是否符合条件：MIC ≥ 32
                negative_samples = df[df['Label'] == 0]
                invalid_negatives = negative_samples[negative_samples[mic_column] < self.mic_threshold_low]
                
                # 排除锚点数据的检查，因为它是强制添加的
                if sequence_column:
                    invalid_negatives = invalid_negatives[invalid_negatives[sequence_column] != self.anchor_sequence]
                
                if len(invalid_negatives) > 0:
                    print(f"  警告：发现 {len(invalid_negatives)} 个不符合条件的负样本")
                else:
                    print(f"  所有负样本都符合条件 (MIC ≥ {self.mic_threshold_low})，排除锚点数据")
            
            # 显示前5行数据预览
            print(f"\n前5行数据预览：")
            preview_columns = ['Label']
            if sequence_column:
                preview_columns.append(sequence_column)
            if mic_column and mic_column in df.columns:
                preview_columns.append(mic_column)
            if hemolysis_column and hemolysis_column in df.columns:
                preview_columns.append(hemolysis_column)
            
            print(df[preview_columns].head())
            
            print(f"\n数据集验证完成！")
            return True
            
        except Exception as e:
            print(f"验证过程中出现错误：{str(e)}")
            return False
    
    def batch_process(self, input_folder, output_folder):
        """
        批量处理文件夹中的所有CSV文件
        
        Args:
            input_folder: 输入文件夹路径
            output_folder: 输出文件夹路径
        """
        # 确保输出文件夹存在
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            print(f"创建输出文件夹: {output_folder}")
        
        # 获取输入文件夹中的所有CSV文件
        csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]
        print(f"找到 {len(csv_files)} 个CSV文件待处理")
        
        # 逐个处理文件
        for csv_file in csv_files:
            input_path = os.path.join(input_folder, csv_file)
            output_path = os.path.join(output_folder, f"processed_{csv_file}")
            
            print(f"\n{'-'*50}")
            print(f"处理文件: {csv_file}")
            print(f"{'-'*50}")
            
            # 处理文件
            df = self.clean_and_label_data(input_path, output_path)
            
            if df is not None:
                print("处理完成！")
            else:
                print("处理失败！")

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='抗菌肽数据处理工具')
    
    # 定义子命令
    subparsers = parser.add_subparsers(dest='command', help='可用命令')
    
    # 处理单个文件命令
    process_parser = subparsers.add_parser('process', help='处理单个CSV文件')
    process_parser.add_argument('-i', '--input', required=True, help='输入CSV文件路径')
    process_parser.add_argument('-o', '--output', required=True, help='输出CSV文件路径')
    
    # 验证文件命令
    validate_parser = subparsers.add_parser('validate', help='验证数据集')
    validate_parser.add_argument('-d', '--dataset', required=True, help='数据集文件路径')
    
    # 批量处理命令
    batch_parser = subparsers.add_parser('batch', help='批量处理文件夹中的CSV文件')
    batch_parser.add_argument('-i', '--input_folder', required=True, help='输入文件夹路径')
    batch_parser.add_argument('-o', '--output_folder', required=True, help='输出文件夹路径')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 创建处理器实例
    processor = AMPDataProcessor()
    
    # 根据命令执行相应操作
    if args.command == 'process':
        processor.clean_and_label_data(args.input, args.output)
    elif args.command == 'validate':
        processor.validate_dataset(args.dataset)
    elif args.command == 'batch':
        processor.batch_process(args.input_folder, args.output_folder)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
