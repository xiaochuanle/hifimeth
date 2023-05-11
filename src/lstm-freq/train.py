# -*- coding: utf-8 -*-

import numpy as np
import random
import torch
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader

import kmer_features_dataset
from kmer_features_dataset import KmerFeaturesDataset

import model_lstm
from model_lstm import BiLSTM_Attention

batch_size = 256 
num_positions = 100
motif_len = 2
num_classes = 2

train_data_path = ''
test_data_path = ''
if 0:
    train_data_path = '/data2/GpC/chenying/small_features/train.bin.npy'
    test_data_path = '/data2/GpC/chenying/small_features/test.bin.npy'
else:
    test_data_path = '/data2/GpC/chenying/features/test.bin.npy'
    train_data_path = '/data2/GpC/chenying/features/train.bin.npy'

dl_workers = 16
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#device = "cpu"

learning_rate = 0.001
weight_decay = 0.00001
hidden_channels = 128
num_epochs = 5

def fix_seeds():
    seed = 42
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.cuda.manual_seed_all(seed)

def eval_model(model, data_loader, max_batch):
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        bi = 0
        for fx, rx, y in tqdm(data_loader):
            fx = fx.to(device)
            rx = rx.to(device)
            y = y.to(device)
            out = model(fx, rx)
            out = torch.softmax(out, 1)
            pred = out.argmax(dim=1)
            correct += int((pred == y).sum())
            total += len(y)
            bi += 1
            if bi == max_batch:
                break
    acc = correct / total
    print('correct = %d, total = %d, acc = %f' % (correct, total, acc))

def train():
    fix_seeds()
    
    train_data = KmerFeaturesDataset(train_data_path, num_positions, motif_len)
    train_loader = DataLoader(train_data, batch_size = batch_size, shuffle = True, num_workers = dl_workers)

    test_data = KmerFeaturesDataset(test_data_path, num_positions, motif_len)
    test_loader = DataLoader(test_data, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

    model = BiLSTM_Attention(num_positions, train_data.node_features, device = device)
    print(model)
    model.to(device)

    #optimizer = torch.optim.Adam(model.parameters(), lr = learning_rate, weight_decay = weight_decay)
    optimizer = torch.optim.SGD(model.parameters(), lr = 1, weight_decay = 0.0001, momentum = 0.9, nesterov=True)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size = 10, gamma = 0.1)
    #scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=1, T_mult=2)

    criterion = torch.nn.CrossEntropyLoss()

    for epoch in range(num_epochs):
        print('Epoch %d' % (epoch + 1))
        model.train()
        bi = 0
        total_bi = len(train_loader)
        for fx, rx, y in tqdm(train_loader):
            optimizer.zero_grad()
            fx = fx.to(device)
            rx = rx.to(device)
            y = y.to(device)
            out = model(fx, rx)
            loss = criterion(out, y)

            loss.backward()
            #torch.nn.utils.clip_grad_norm_(parameters = model.parameters(), max_norm = 0.1, norm_type=2)
            optimizer.step()          
            bi += 1
            if bi % 1000 == 0 or bi == total_bi:
                print('eval at step %d' % (bi))
                eval_model(model, test_loader, -1)
                model.train()
        eval_model(model, train_loader, 200)

        scheduler.step()

if __name__ == "__main__":
    if device == 'cpu':
        print('CUDA is not available')
    else:
        print('CUDA is available')

    train()
