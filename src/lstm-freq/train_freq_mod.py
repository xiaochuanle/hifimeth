# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import random
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
import torch.nn.functional as F

import freq_dataset
from freq_dataset import FreqFeaturesDataset
from freq_dataset import mod_site_features

import model_lstm
from model_lstm import AggrAttRNN

from sklearn import metrics

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

def fix_seeds():
    seed = 42
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.cuda.manual_seed_all(seed)

def eval_model(model, criterion, data_loader, max_batch, device):
    model.eval()

    eval_loss = 0.0
    bi = 0
    with torch.no_grad():
        for x, y in tqdm(data_loader):
            x = x.to(device)
            y = y.to(device)
            out = model(x)
            #out = torch.squeeze(out, dim=1)
            #torch.clip_(out, 0.0, 1.0)
            #print(out.shape, data.y.shape)
            loss = criterion(out, y)
            eval_loss = eval_loss + loss.item()
            if bi == 0:
                for i in range(20):
                    print('%d\t%f\t%f' % (i, out[i], y[i]))
            bi += 1
            if bi == max_batch:
                break
    avg_loss = eval_loss / bi 
    print('Loss: {:.4f}'.format(avg_loss))

def train(args):
    fix_seeds()

    model_dir = args.model_dir 
    train_data_path = args.train_file 
    test_data_path = args.test_file
    test_data_path2 = args.test_file2
    kmer_size = args.kmer_size 
    batch_size = args.batch_size
    dl_workers = args.dl_workers
    num_epochs = args.epoch

    print('model dir: %s' % (model_dir))
    print('train data: %s' % (train_data_path))
    print('test data: %s' % (test_data_path))
    if test_data_path2 is not None:
        print('test data2: %s' % (test_data_path2))
    print('kmer size: %d' % (kmer_size))
    print('batch size: %d' % (batch_size))
    print('dl workers: %d' % (dl_workers))
    print('epoches: %d' % (num_epochs))

    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    node_features = mod_site_features
    model = AggrAttRNN(kmer_size, node_features, device = device)
    print(model)
    model.to(device)

    #optimizer = torch.optim.SGD(model.parameters(), lr = 0.001, weight_decay = 0.001, momentum = 0.9, nesterov=True)
    #optimizer = torch.optim.Adam(model.parameters(), lr = 0.1, weight_decay = 0.00001)
    optimizer = torch.optim.Adam(model.parameters(), lr = 0.01, weight_decay = 0.00001)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size = 10, gamma = 0.5)

    #criterion = torch.nn.CrossEntropyLoss()
    criterion = torch.nn.MSELoss()
    #criterion = torch.nn.L1Loss()
    
    test_data = FreqFeaturesDataset(test_data_path, kmer_size, 1)
    test_loader = DataLoader(test_data, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

    test_data2 = None 
    test_loader2 = None 
    if test_data_path2 is not None:
        print('Load test data2 %s' % (test_data_path2))
        test_data2 = FreqFeaturesDataset(test_data_path2, kmer_size, 1)
        test_loader2 = DataLoader(test_data2, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

    print('Eval at pretrained model')
    eval_model(model, criterion, test_loader, 200, device)
    if test_loader2 is not None:
        eval_model(model, criterion, test_loader2, 200, device)

    print('Load train data %s' % (train_data_path))
    train_data = FreqFeaturesDataset(train_data_path, kmer_size, 1)
    train_loader = DataLoader(train_data, batch_size = batch_size, shuffle = True, num_workers = dl_workers)

    for epoch in range(num_epochs):
        print('Epoch %d' % (epoch + 1))
        model.train()
        bi = 0
        total_bi = len(train_loader)
        for fx, y in tqdm(train_loader):
            optimizer.zero_grad()
            fx = fx.to(device)
            y = y.to(device)
            out = model(fx)
            loss = criterion(out, y)
            loss.backward()
            #torch.nn.utils.clip_grad_norm_(parameters = model.parameters(), max_norm = 0.1, norm_type=2)
            optimizer.step()          
            bi += 1
            '''
            if bi % 1000 == 0 or bi == total_bi:
                print('eval at step %d' % (bi))
                eval_model(model, test_loader, -1)
                model.train()
            '''

        torch.save(model.state_dict(), '{}/kmer_{}_epoch_{}.ckpt'.format(model_dir, kmer_size, epoch+1))
        print('eval at test data')
        eval_model(model, criterion, test_loader, -1, device)
        if test_loader2 is not None:
            print('eval at test data2')
            eval_model(model, criterion, test_loader2, -1, device)
        #print('eval at train data')
        #eval_model(model, criterion, train_loader, 200, device)

        scheduler.step()

def main():
    parser = argparse.ArgumentParser("train 5mc-gnn model")
    parser.add_argument('--train_file', type=str, required=True)
    parser.add_argument('--test_file', type=str, required=True)
    parser.add_argument('--test_file2', type=str, required=False)
    parser.add_argument('--model_dir', type=str, required=True)
    parser.add_argument('--kmer_size', type=int, default=11, required=False,
                          help="len of kmer. default: 11")
    parser.add_argument('--batch_size', type=int, default=4096, required=False,
                        help="number of examples in one training batch. default: 1024")
    parser.add_argument('--dl_workers', type=int, default=8, required=False,
                        help="number of CPU threads for generating Graph, default: 8")
    parser.add_argument('--epoch', type=int, default=50, required=False,
                        help="number of training epoch, default: 2");

    args = parser.parse_args()

    train(args)

if __name__ == "__main__":
    if device == torch.device("cpu"):
        print('CUDA is not available')
    else:
        print('CUDA is available')

    main()
