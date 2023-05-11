# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import random
import torch
from tqdm import tqdm
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import DataLoader
import torch.nn.functional as F

import kmer_features_dataset
from kmer_features_dataset import KmerFeaturesDataset
from kmer_features_dataset import kmer_base_features

import kmer_model
from kmer_model import KmerModel

import kmer_util
from kmer_util import eval_model

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

def fix_seeds():
    seed = 42
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.cuda.manual_seed_all(seed)

def train(args):
    fix_seeds()

    pretrained_model_path = args.pretrained_model
    model_dir = args.model_dir 
    train_data_path = args.train_file 
    test_data_path = args.test_file
    test_data_path2 = args.test_file2
    test_data_path3 = args.test_file3
    kmer_size = args.kmer_size 
    batch_size = args.batch_size
    dl_workers = args.dl_workers
    num_epochs = args.epoch
    num_classes = 2
    hidden_channels = 128

    if pretrained_model_path is not None:
        print('pretrained model: %s' % (pretrained_model_path))
    print('model dir: %s' % (model_dir))
    print('train data: %s' % (train_data_path))
    print('test data: %s' % (test_data_path))
    if test_data_path2 is not None:
        print('test data2: %s' % (test_data_path2))
    if test_data_path3 is not None:
        print('test data3: %s' % (test_data_path3))
    print('kmer size: %d' % (kmer_size))
    print('batch size: %d' % (batch_size))
    print('dl workers: %d' % (dl_workers))
    print('epoches: %d' % (num_epochs))

    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    node_features = kmer_base_features
    model = KmerModel(num_node_features = node_features, hidden_channels = hidden_channels, num_classes = num_classes)
    if pretrained_model_path is not None:
        print('Load pretrained model %s' % (pretrained_model_path))
        model.load_state_dict(torch.load(pretrained_model_path, map_location=device))
    print(model)
    model.to(device)

    optimizer = torch.optim.SGD(model.parameters(), lr = 0.1, weight_decay = 0.00001, momentum = 0.9, nesterov=True)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size = 1, gamma = 0.1)

    criterion = torch.nn.CrossEntropyLoss()
    
    test_data = KmerFeaturesDataset(test_data_path, kmer_size, 1)
    test_loader = DataLoader(test_data, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

    test_data2 = None 
    test_loader2 = None 
    if test_data_path2 is not None:
        print('Load test data2 %s' % (test_data_path2))
        test_data2 = KmerFeaturesDataset(test_data_path2, kmer_size, 1)
        test_loader2 = DataLoader(test_data2, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

    test_data3 = None 
    test_loader3 = None 
    if test_data_path3 is not None:
        print('Load test data3 %s' % (test_data_path3))
        test_data3 = KmerFeaturesDataset(test_data_path3, kmer_size, 1)
        test_loader3 = DataLoader(test_data3, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

    print('Eval at pretrained model')
    eval_model(model, criterion, test_loader, 200, device)
    if test_loader2 is not None:
        print('eval at test data 2')
        eval_model(model, criterion, test_loader2, 200, device)
    if test_loader3 is not None:
        print('eval at test data 3')
        eval_model(model, criterion, test_loader3, 200, device)

    print('Load train data %s' % (train_data_path))
    train_data = KmerFeaturesDataset(train_data_path, kmer_size, 1)
    train_loader = DataLoader(train_data, batch_size = batch_size, shuffle = True, num_workers = dl_workers)

    for epoch in range(num_epochs):
        print('Epoch %d' % (epoch + 1))
        model.train()
        bi = 0
        total_bi = len(train_loader)
        for data in tqdm(train_loader):
            optimizer.zero_grad()
            data.to(device)
            out = model(data.x, data.edge_index, data.batch)
            loss = criterion(out, data.y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(parameters = model.parameters(), max_norm = 0.1, norm_type=2)
            optimizer.step()          
            bi += 1
            if bi % 5000 == 0 or bi == total_bi:
                print('eval at step %d' % (bi))
                eval_model(model, criterion, test_loader, 200, device)
                if test_loader2 is not None:
                    print('eval at test data 2')
                    eval_model(model, criterion, test_loader2, 200, device)
                if test_loader3 is not None:
                    print('eval at test data 3')
                    eval_model(model, criterion, test_loader3, 200, device)
                model.train()
        torch.save(model.state_dict(), '{}/kmer_{}_epoch_{}.ckpt'.format(model_dir, kmer_size, epoch+1))
        eval_model(model, criterion, train_loader, 200, device)

        scheduler.step()

def main():
    parser = argparse.ArgumentParser("train 5mc-gnn model")
    parser.add_argument('--pretrained_model', type=str, required=False,
                        help="Training the model with a pretrained one")
    parser.add_argument('--train_file', type=str, required=True)
    parser.add_argument('--test_file', type=str, required=True)
    parser.add_argument('--test_file2', type=str, required=False)
    parser.add_argument('--test_file3', type=str, required=False)
    parser.add_argument('--model_dir', type=str, required=True)
    parser.add_argument('--kmer_size', type=int, default=400, required=False,
                          help="len of kmer. default: 400")
    parser.add_argument('--batch_size', type=int, default=256, required=False,
                        help="number of examples in one training batch. default: 256")
    parser.add_argument('--dl_workers', type=int, default=16, required=False,
                        help="number of CPU threads for generating Graph, default: 16")
    parser.add_argument('--epoch', type=int, default=2, required=False,
                        help="number of training epoch, default: 2");

    args = parser.parse_args()

    train(args)

if __name__ == "__main__":
    if device == torch.device("cpu"):
        print('CUDA is not available')
    else:
        print('CUDA is available')

    main()
