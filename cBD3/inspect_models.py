import joblib
import numpy as np

def inspect_model(model_path):
    """
    检查模型文件的结构和属性
    """
    try:
        print(f"\n正在检查模型: {model_path}")
        model = joblib.load(model_path)
        
        # 检查模型类型
        print(f"模型类型: {type(model).__name__}")
        
        # 如果是pipeline，检查组件
        if hasattr(model, 'named_steps'):
            print("\nPipeline组件:")
            for name, step in model.named_steps.items():
                print(f"  - {name}: {type(step).__name__}")
            
            # 检查分类器
            clf = model.named_steps.get('clf')
            if clf:
                print("\n分类器信息:")
                print(f"  类型: {type(clf).__name__}")
                if hasattr(clf, 'n_estimators'):
                    print(f"  估计器数量: {clf.n_estimators}")
                if hasattr(clf, 'feature_importances_'):
                    print(f"  特征数量: {clf.feature_importances_.shape[0]}")
                    # 显示前5个最重要的特征
                    importances = clf.feature_importances_
                    indices = np.argsort(importances)[::-1]
                    print("  前5个重要特征索引:", indices[:5])
                    print("  前5个重要特征值:", importances[indices[:5]])
        
        print("\n模型检查完成！")
        return model
        
    except Exception as e:
        print(f"检查模型时出错: {str(e)}")
        return None

if __name__ == "__main__":
    # 检查两个模型文件
    inspect_model("rf_amp_model.joblib")
    inspect_model("rf_amp_model_constrained.joblib")
