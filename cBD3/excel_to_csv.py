import pandas as pd
import os

# 定义文件路径
excel_file = 'd:/trea/trae_projects/cBD3-ABU/my_500_peptides.xlsx'
csv_file = 'd:/trea/trae_projects/cBD3-ABU/my_500_peptides.csv'

# 检查Excel文件是否存在
if not os.path.exists(excel_file):
    print(f"错误：找不到文件 {excel_file}")
    exit(1)

try:
    # 读取Excel文件
    print(f"正在读取Excel文件：{excel_file}")
    df = pd.read_excel(excel_file)
    
    # 显示文件的基本信息
    print(f"文件包含 {len(df)} 行和 {len(df.columns)} 列")
    print("\n列名：")
    for col in df.columns:
        print(f"- {col}")
    
    print("\n前5行数据预览：")
    print(df.head())
    
    # 转换为CSV
    print(f"\n正在转换为CSV文件：{csv_file}")
    df.to_csv(csv_file, index=False, encoding='utf-8-sig')
    
    print(f"\n转换完成！CSV文件已保存至：{csv_file}")
    print(f"CSV文件大小：{os.path.getsize(csv_file) / 1024:.2f} KB")
    
except Exception as e:
    print(f"转换过程中出现错误：{str(e)}")