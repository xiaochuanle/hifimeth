from sklearn import metrics
import torch
from tqdm import tqdm
import numpy as np

def eval_model(model, criterion, data_loader, max_batch, device):
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        loss_batch_list, label_batch_list, pred_batch_list = [], [], []
        bi = 0
        for data in tqdm(data_loader):
            data.to(device)
            out = model(data.x, data.edge_index, data.batch)
            loss = criterion(out, data.y)
            out = torch.softmax(out, 1)
            pred = out.argmax(dim=1)
            correct += int((pred == data.y).sum())
            total += len(data.y)
            bi += 1
            
            loss_batch_list.append(loss.item())
            y = data.y.cpu()
            label_batch_list += y.tolist()
            pred = pred.cpu()
            pred_batch_list += pred.tolist()

            if bi == max_batch:
                break
    acc = correct / total
    print('correct = %d, total = %d, acc = %f' % (correct, total, acc))

    acc = metrics.accuracy_score(label_batch_list, pred_batch_list)
    prec = metrics.precision_score(label_batch_list, pred_batch_list)
    recall = metrics.recall_score(label_batch_list, pred_batch_list)
    mean_loss = np.mean(loss_batch_list)
    f1 = metrics.f1_score(label_batch_list, pred_batch_list)
    print('Loss: {:.4f}; Acc: {:.4f}; Prec: {:.4f}, Recall: {:.4f}, F1-score: {:.4f}'.format(mean_loss, acc, prec, recall, f1))
