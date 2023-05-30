#!/usr/bin/env python3

import pandas as pd
import sys
from multiprocessing import Pool

INPUT = sys.argv[1]
CORES = int(sys.argv[2])
#pool = Pool(CORES)

# 读取tsv文件，指定分隔符为\t，列名为pred和label
DF = pd.read_csv(INPUT, sep="\t", names=["pred", "label"])
# 定义一个函数，根据给定的阈值，将预测评分转换为二分类标签
def binarize(pred, threshold):
    if pred >= threshold:
        return 1
    else:
        return 0
# 定义一个空列表，用于存储不同阈值下的结果

# 定义一个阈值列表，从0到1，步长为0.01
thresholds = [i/100 for i in range(101)]
# 遍历阈值列表，对每个阈值，计算TPR，FPR和F1 score，并添加到结果列表中
def ROC(t, df=DF.copy()):
    # 应用binarize函数，得到二分类标签列
    df["bin"] = df["pred"].apply(binarize, args=(t,))
    # 计算真阳性（TP），真阴性（TN），假阳性（FP）和假阴性（FN）的数量
    TP = ((df["bin"] == 1) & (df["label"] == 1)).sum()
    TN = ((df["bin"] == 0) & (df["label"] == 0)).sum()
    FP = ((df["bin"] == 1) & (df["label"] == 0)).sum()
    FN = ((df["bin"] == 0) & (df["label"] == 1)).sum()
    # 计算真阳性率（TPR），假阳性率（FPR）和F1 score
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    F1 = 2 * TP / (2 * TP + FP + FN)
    Precision = TP / (TP+FP)
    Accuracy = (TP + TN) / (TP + TN + FP + FN)
    Recall = TP / (TP + FN)
    # 将阈值和对应的指标添加到结果列表中
    return [t, TP, TN, FP, FN, TPR, FPR, F1, Precision, Accuracy, Recall]
if __name__ == '__main__':
    with Pool(CORES) as p:
        results = p.map(ROC, thresholds)
# 将结果列表转换为pandas数据框，并指定列名
results_df = pd.DataFrame(results, columns=["threshold", 'TP', 'TN', 'FP', 'FN', "TPR", "FPR", "F1", "Precision", "Accuracy", "Recall"])
# 将结果数据框保存到一个tsv文件中
results_df.to_csv(f"{INPUT.rstrip('tsv')}ROC_F1.tsv", sep="\t", index=False)
# 导入sklearn库中的roc_auc_score函数
from sklearn.metrics import roc_auc_score
# 计算ROC曲线的AUC值，并打印出来
auc = roc_auc_score(DF["label"], DF["pred"])
print(f"{INPUT} ROC曲线的AUC值为：{auc:.4f}")
