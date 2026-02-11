# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import random
import sys
import torch
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader

from sample_dataset import KmerFeaturesDataset
from sample_dataset import load_features
from sample_dataset import KMER_BASE_FEATURES

from model_cnn import DNAModNet

def export_libtorch_model(model_path, scripted_model_path, kmer):
    print('export to libtorch model')
    
    cpu_device = torch.device("cpu")
    model = DNAModNet(kmer)
    model.load_state_dict(torch.load(model_path, map_location=cpu_device), strict=False)
    print(model)
    model.to(cpu_device) 
    model.eval()

    scripted_model = torch.jit.script(model)
    scripted_model = torch.jit.optimize_for_inference(scripted_model.eval())
    scripted_model.save(scripted_model_path)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('USAGE:')
        print(f"  {sys.argv[0]} model output-model-dir kmer")
        sys.exit (1)
    model_path = sys.argv[1]
    output_model = sys.argv[2]
    kmer = int(sys.argv[3])
    export_libtorch_model(model_path, output_model, kmer)