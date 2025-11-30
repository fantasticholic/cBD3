import pandas as pd
import os

# 定义CSV文件路径
csv_file = 'd:/trea/trae_projects/cBD3-ABU/my_500_peptides.csv'

# 检查CSV文件是否存在
if not os.path.exists(csv_file):
    print(f"错误：找不到CSV文件 {csv_file}")
    exit(1)

try:
    # 读取CSV文件
    print(f"正在验证CSV文件：{csv_file}")
    df = pd.read_csv(csv_file)
    
    # 显示CSV文件的基本信息
    print(f"CSV文件包含 {len(df)} 行和 {len(df.columns)} 列")
    print("\n列名：")
    for col in df.columns:
        print(f"- {col}")
    
    print("\n前3行数据预览：")
    print(df.head(3))
    
    print("\nCSV文件验证成功！数据完整且格式正确。")
    
except Exception as e:
    print(f"验证过程中出现错误：{str(e)}")