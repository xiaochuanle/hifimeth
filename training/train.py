# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import random
import torch
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F

import sample_dataset
from sample_dataset import KmerFeaturesDataset
from sample_dataset import load_features
from sample_dataset import KMER_BASE_FEATURES

import model_cnn
from model_cnn import DNAModNet

from sklearn import metrics

dl_workers = 48
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#device = "cpu"

torch.set_float32_matmul_precision('high')

BATCH_SIZE = 512

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
    features_path = args.feature
    train_samples_path = args.train_sample
    eval_samples_path = args.eval_sample
    kmer_size = args.kmer_size 
    batch_size = args.batch_size
    dl_workers = args.dl_workers
    num_epochs = args.epoch

    model = DNAModNet(kmer_size)
    print(model)
    model.to(device)

    optimizer = torch.optim.SGD(model.parameters(), lr = 0.1, weight_decay = 0.00001, momentum = 0.9, nesterov=True)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size = 1, gamma = 0.2)
    criterion = torch.nn.CrossEntropyLoss()

    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    features_path = '//data3/cy/cpg/features/features'
    offsets_path = '/data3/cy/cpg/train-samples/cpg.offsets'
    train_sample_path = '/data3/cy/cpg/train-samples/cpg.samples'

    print(f'Load features from {features_path}')
    features = load_features(features_path)
    print(f'Load train samples from {train_sample_path}')
    train_dataset = KmerFeaturesDataset(features, train_sample_path, offsets_path, kmer_size)
    train_loader = DataLoader(train_dataset, batch_size = batch_size, shuffle = True, num_workers = dl_workers)
    print('Done')    

    for epoch in range(num_epochs):
        print('Epoch %d' % (epoch + 1))
        if epoch >= 0:
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
                optimizer.step()
                bi += 1
                if bi % 10000 == 0:
                    model.eval()
                    model_path = f"{model_dir}/{epoch}_{bi}.ckpt"
                    torch.save(model.state_dict(), model_path)
                    model.train()
                    
        if epoch >= 0:
            torch.save(model.state_dict(), '{}/kmer_{}_epoch_{}.ckpt'.format(model_dir, kmer_size, epoch+1))
        scheduler.step()

        current_lr = optimizer.param_groups[0]['lr']
        print(f'Epoch {epoch+1:2d}: Learning Rate = {current_lr:.4f}')

def main():
    parser = argparse.ArgumentParser("train 5mc-gnn model")
    parser.add_argument('--pretrained_model', type=str, required=False,
                        help="Training the model with a pretrained one")
    parser.add_argument('--feature', type=str, required=False)
    parser.add_argument('--train_sample', type=str, required=False)
    parser.add_argument('--eval_sample', type=str, required=False)
    parser.add_argument('--model_dir', type=str, required=True)
    parser.add_argument('--kmer_size', type=int, required=True,
                          help="len of kmer.")
    parser.add_argument('--batch_size', type=int, default=BATCH_SIZE, required=False,
                        help="number of examples in one training batch. default: %d" % (BATCH_SIZE))
    parser.add_argument('--dl_workers', type=int, default=dl_workers, required=False,
                        help="number of CPU threads for generating Graph, default: 16")
    parser.add_argument('--epoch', type=int, default=3, required=False,
                        help="number of training epoch, default: 3")

    args = parser.parse_args()

    train(args)

if __name__ == "__main__":
    if device == torch.device("cpu"):
        print('CUDA is not available')
    else:
        print('CUDA is available')

    main()