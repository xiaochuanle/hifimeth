from threading import Thread 

import argparse
import numpy as np
import os 
import sys
import time 
import torch
import torch.nn as nn
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader

import model_lstm
from model_lstm import AggrAttRNN

import freq_dataset
from freq_dataset import mod_site_features
from freq_dataset import mod_site_input_features

#device = torch.device("cpu")
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

torch.multiprocessing.set_sharing_strategy('file_system')

sample_batch_size = 10000000

def get_time_string():
    localtime = time.asctime(time.localtime(time.time()))
    return localtime

def assemble_one_graph_data_call(line, num_sites):

    line = line.strip()
    cols = line.split()
    features = [float(x) for x in cols[4:]]

    num_nodes = num_sites
    nodes = np.zeros([num_nodes, mod_site_features], dtype=np.float32)
    idx = 0
    for i in range(num_sites):
        nodes[i][0:mod_site_input_features] = features[idx:idx+mod_site_input_features]
        nodes[i][mod_site_input_features] = 0.0
        nodes[i][mod_site_input_features + 1] = 1.0
        idx += mod_site_input_features
    target_pos = num_sites // 2
    nodes[target_pos][mod_site_input_features] = 1.0
    nodes[target_pos][mod_site_input_features+1] = 0.0
    t_nodes = torch.from_numpy(nodes)

    return t_nodes

class FreqFeaturesDatasetCall(Dataset):
    def __init__(self, features_path, feature_stream, num_sites):
        self.data_root = os.path.dirname(features_path)
        self.features_path = features_path
        self.num_sites = num_sites
        super(FreqFeaturesDatasetCall, self).__init__()

        self.line_list = []
        cnt = 0
        while True:
            line = feature_stream.readline()
            if line == '':
                break
            
            self.line_list.append(line)
            cnt += 1
            if cnt == sample_batch_size:
                break 

        self.n_samples = cnt 

    def __len__(self):
        return self.n_samples

    def __getitem__(self, idx):
        assert idx < len(self.line_list)
        return assemble_one_graph_data_call(self.line_list[idx], self.num_sites)

def dump_freq_results(line_list, freq_list, out):
    assert len(line_list) == len(freq_list)
    n_samples = len(line_list)
    for i in range(n_samples):
        line = line_list[i]
        line = line.strip()
        cols = line.split('\t')
        soff = int(cols[1])
        send = soff + 1
        prob = float(freq_list[i][0])
        line = '\t'.join((cols[0], str(soff), str(send), str(prob), cols[2], cols[3]))
        line = line + '\n'
        out.write(line)


def mod_freq_call(args):

    model_path = args.model_path
    kmer_size = args.kmer_size 
    features_path = args.input
    output = args.out
    dl_workers = args.dl_workers
    batch_size = args.batch_size

    print('model: %s' % (model_path))
    print('kmer size: %d' % (kmer_size))
    print('features path: %s' % (features_path))
    print('output: %s' % (output))
    print('dl_workers: %d' % (dl_workers))
    print('batch size: %d' % (batch_size))

    node_features = mod_site_features
    #model = BiLSTM_Attention(kmer_size, node_features, device = device)
    model = AggrAttRNN(kmer_size, node_features, device = device)
    print('Load model %s' % (model_path))
    model.load_state_dict(torch.load(model_path, map_location=device))
    #print(model)
    model.to(device)
    model.eval()

    sample_file = open(features_path, 'r')
    out_file = open(output, 'w')
    samples_processed = 0
    while True:
        dataset = FreqFeaturesDatasetCall(features_path, sample_file, kmer_size)
        num_samples = len(dataset)
        if num_samples == 0:
            break
        now = get_time_string()
        print('[%s] Load %d samples' % (now, num_samples))

        data_loader = DataLoader(dataset, batch_size = batch_size, shuffle = False, num_workers = dl_workers)
        freq_list = [0.0] * len(dataset)
        freq_i = 0
        with torch.no_grad():
            for data in tqdm(data_loader):
                data = data.to(device)
                out = model(data)
                out = out.to('cpu')
                y = out.numpy()
                y = np.clip(y, 0, 1)
                y = y.tolist()
                n = len(y)
                freq_list[freq_i:freq_i+n] = y[0:n]
                freq_i = freq_i + n

        now = get_time_string()
        print('[%s] dump freq results' % (now))
        dump_freq_results(dataset.line_list, freq_list, out_file)
        samples_processed += num_samples
        now = get_time_string()
        print('[%s] %d samples processed' % (now, samples_processed))

    sample_file.close()
    out_file.close()

def main():
    parser = argparse.ArgumentParser("5mc_freq.py")
    parser.add_argument('--model_path', type=str, required=True)
    parser.add_argument('--kmer_size', type=int, default=11, required=False,
                          help="len of kmer. default: 11")
    parser.add_argument('--input', type=str, required = True)
    parser.add_argument('--out', type=str, required=True)
    parser.add_argument('--dl_workers', type = int, default = 16, required = False)
    parser.add_argument('--batch_size', type=int, default=512, required=False)

    args = parser.parse_args()

    mod_freq_call(args)

if __name__ == "__main__":
    main()
