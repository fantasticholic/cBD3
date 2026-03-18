import os
import sys

# 添加当前目录到Python路径，确保可以导入amp_data_processor
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from amp_data_processor import AMPDataProcessor

def create_example_amp_file(output_path):
    """
    创建一个示例的AMP数据CSV文件，模拟AMP1023.csv的格式
    """
    header = "ID,COMPLEXITY,Peptide_Name,N TERMINUS,Sequence,C TERMINUS,SYNTHESIS TYPE,TARGET GROUP,TARGET OBJECT,MIC (µg/M),Target_Organism,Hemolysis (%),Hemolysis_Concentration (µM),CC50 (µg/mL),Cell_Line,Source_Paper (DOI)"
    
    # 示例数据行
    example_data = [
        "1,Monomer,Example_Peptide_1,,KAWNLRGSAREKAIKNEKLYIFATSGKLAALKPK,,Synthetic,Gram+,Bacteria,8,Staphylococcus aureus,5,10,,,",
        "2,Monomer,Example_Peptide_2,,GLGKFLHSAKRFGKAFVGEAMNS,AMD,Synthetic,Gram+,Bacteria,>32,Staphylococcus aureus,15,20,,,",
        "3,Monomer,Example_Peptide_3,,GIGKFIHAAKKFGKLFIGEIMNS,AMD,Synthetic,Gram+,Bacteria,4,Staphylococcus aureus,2,5,,,",
        "4,Monomer,Example_Peptide_4,,FKRIVQRIKDFLR,,Synthetic,Gram+,Bacteria,64,Staphylococcus aureus,10,15,,,",
        "5,Monomer,Example_Peptide_5,,ALWKTLLKKVLKAAAKAALKAVLVGANA,,Synthetic,Gram+,Bacteria,2,Staphylococcus aureus,0.5,0.5,,,",
        "6,Monomer,Example_Peptide_6,,GIGKFLHSAKKFGKAWVGEIMNS,AMD,Synthetic,Gram+,Bacteria,40,Staphylococcus aureus,25,30,,,",
        "7,Monomer,Example_Peptide_7,,GLGKFIHAAKRFGKLFVGEAMNS,AMD,Synthetic,Gram+,Bacteria,7,Staphylococcus aureus,3,8,,,",
        "8,Monomer,Example_Peptide_8,,KAWNLRGSAREKAIKNEKLYIFATSGKLAALK,AMD,Synthetic,Gram+,Bacteria,12,Staphylococcus aureus,1,2,,,",
        "9,Monomer,Example_Peptide_9,,FKRIVQRIKDFLRKAWN,,Synthetic,Gram+,Bacteria,36,Staphylococcus aureus,18,22,,,",
        "10,Monomer,Example_Peptide_10,,GIGKFLHSAKKFGKAWVGEIMNSKAWN,,Synthetic,Gram+,Bacteria,5,Staphylococcus aureus,4,6,,,",
    ]
    
    # 写入文件
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(header + '\n')
        for line in example_data:
            f.write(line + '\n')
    
    print(f"示例文件已创建: {output_path}")
    return output_path

def main():
    """
    示例：如何使用AMPDataProcessor处理和验证AMP数据
    """
    print("="*60)
    print("抗菌肽数据处理工具使用示例")
    print("="*60)
    
    # 创建示例数据文件
    example_file = create_example_amp_file("example_amp_data.csv")
    processed_file = "processed_amp_data.csv"
    
    # 创建处理器实例
    processor = AMPDataProcessor()
    
    print("\n" + "="*60)
    print("1. 处理单个AMP数据文件")
    print("="*60)
    
    # 处理数据
    df = processor.clean_and_label_data(example_file, processed_file)
    
    if df is not None:
        print("\n处理结果摘要：")
        print(f"- 输入文件: {example_file}")
        print(f"- 输出文件: {processed_file}")
        print(f"- 总样本数: {len(df)}")
        print(f"- 正样本数 (好肽): {(df['Label'] == 1).sum()}")
        print(f"- 负样本数 (差肽): {(df['Label'] == 0).sum()}")
    
    print("\n" + "="*60)
    print("2. 验证处理后的数据")
    print("="*60)
    
    # 验证数据
    is_valid = processor.validate_dataset(processed_file)
    print(f"\n数据验证{'通过' if is_valid else '失败'}")
    
    print("\n" + "="*60)
    print("3. 命令行使用方法")
    print("="*60)
    
    print("处理单个文件:")
    print(f"python amp_data_processor.py process -i {example_file} -o output.csv")
    
    print("\n验证数据集:")
    print(f"python amp_data_processor.py validate -d {processed_file}")
    
    print("\n批量处理文件夹:")
    print("python amp_data_processor.py batch -i input_folder -o output_folder")
    
    print("\n" + "="*60)
    print("使用提示:")
    print("="*60)
    print("1. 确保输入文件包含必要的列：Sequence和MIC相关列")
    print("2. 可以自定义阈值参数：mic_threshold_high, mic_threshold_low, hemolysis_threshold")
    print("3. 工具会自动识别包含'MIC'、'Hemolysis'、'Sequence'的列名")
    print("4. 批量处理功能适用于处理多个数据文件")
    print("\n处理完成！")

if __name__ == "__main__":
    main()
