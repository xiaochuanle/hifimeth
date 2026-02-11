import numpy as np
import sys
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score, confusion_matrix

def calculate_binary_metrics(data):
    """
    专用于二分类模型的评估指标计算函数，包含Specificity计算
    
    参数:
    data -- NumPy二维数组，形状为 (n_samples, 至少2列)。
            默认情况下:
            - 第0列应为真实标签 (ground truth)
            - 第1列应为预测标签 (predicted labels)
            可选地，可以包含第2列作为预测为正类的概率（用于计算AUC和AP）
    
    返回:
    dict -- 包含各项指标值的字典，包括specificity
    """
    # 检查输入数据有效性
    if not isinstance(data, np.ndarray):
        raise ValueError("输入数据必须是NumPy数组")
    if data.ndim != 2:
        raise ValueError("输入数据必须是二维数组")
    if data.shape[0] == 0:
        raise ValueError("输入数据不能为空")
    if data.shape[1] < 2:
        raise ValueError("输入数据至少需要2列（真实标签和预测标签）")
    
    # 从数组中提取真实标签和预测标签
    y_true = data[:, 0]  # 第一列是真实标签
    y_pred = data[:, 1]  # 第二列是预测标签
    
    # 确保标签是整数类型（0和1）
    y_true = y_true.astype(int)
    y_pred = y_pred.astype(int)
    
    # 计算基础分类指标
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average='binary', zero_division=0)
    recall = recall_score(y_true, y_pred, average='binary', zero_division=0)
    f1 = f1_score(y_true, y_pred, average='binary', zero_division=0)
    
    # 计算Specificity（需要从混淆矩阵中获取TN和FP）
    cm = confusion_matrix(y_true, y_pred)
    # 确保是二分类混淆矩阵
    if cm.shape == (2, 2):
        TN = cm[0, 0]  # 真负例
        FP = cm[0, 1]  # 假正例
        specificity = TN / (TN + FP) if (TN + FP) > 0 else 0.0
    else:
        # 处理非标准二分类情况
        specificity = 0.0
    
    # 初始化 AUC 和 AP
    auc_value = np.nan
    ap_value = np.nan
    
    # 如果有第三列（预测概率），则计算AUC和AP
    if data.shape[1] >= 3:
        y_prob = data[:, 2]  # 第三列是预测概率
        try:
            auc_value = roc_auc_score(y_true, y_prob)
            ap_value = average_precision_score(y_true, y_prob)
        except Exception as e:
            print(f"计算AUC或AP时出错: {e}")
            print("请确保提供的是预测概率而不是标签，或检查数据格式")
    
    return {
        'accuracy': round(accuracy, 4),
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1_score': round(f1, 4),
        'specificity': round(specificity, 4),  # 新增的specificity
        'auc': round(auc_value, 4) if not np.isnan(auc_value) else auc_value,
        'average_precision': round(ap_value, 4) if not np.isnan(ap_value) else ap_value,
        'n_samples': len(y_true)
    }

# 使用示例
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"USAGE:")
        print(f"{sys.argv[0]} input-prefix num-evals")
        sys.exit(1)
    
    input_prefix = sys.argv[1]
    num_evals = int(sys.argv[2])

    acc = np.zeros(num_evals, dtype=np.float32)
    prec = np.zeros(num_evals, dtype=np.float32)
    recall = np.zeros(num_evals, dtype=np.float32)
    specificity = np.zeros(num_evals, dtype=np.float32)
    f1 = np.zeros(num_evals, dtype=np.float32)
    auc = np.zeros(num_evals, dtype=np.float32)
    ap = np.zeros(num_evals, dtype=np.float32)

    for i in range(num_evals):
        path = f"{input_prefix}.{i}"
        example_data = np.loadtxt(path, dtype=np.float32, delimiter='\t')
    
        # 计算指标
        results = calculate_binary_metrics(example_data)

        for metric, value in results.items():
            if metric == 'accuracy':
                acc[i] = value
            if metric == 'precision':
                prec[i] = value
            if metric == 'recall':
                recall[i] = value
            if metric == 'specificity':
                specificity[i] = value
            if metric == 'f1_score':
                f1[i] = value
            if metric == 'auc':
                auc[i] = value
            if metric == 'average_precision':
                ap[i] = value
    
    print("Accuracy:")
    print(acc)
    mean = np.mean(acc)
    var = np.var(acc)
    print(f"mean: {mean}, var: {var}")

    print("Precision:")
    print(prec)
    mean = np.mean(prec)
    var = np.var(prec)
    print(f"mean: {mean}, var: {var}")

    print('Recall')
    print(recall)
    mean = np.mean(recall)
    var = np.var(recall)
    print(f"mean: {mean}, var: {var}")

    print('specificity')
    print(specificity)
    mean = np.mean(specificity)
    var = np.var(specificity)
    print(f"mean: {mean}, var: {var}")

    print("F1-score:")
    print(f1)
    mean = np.mean(f1)
    var = np.var(f1)
    print(f"mean: {mean}, var: {var}")

    print("auc:")
    print(auc)
    mean = np.mean(auc)
    var = np.var(auc)
    print(f"mean: {mean}, var: {var}")

    print("average_precision:")
    print(ap)
    mean = np.mean(ap)
    var = np.var(ap)
    print(f"mean: {mean}, var: {var}")