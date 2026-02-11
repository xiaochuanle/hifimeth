# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import random
import sys
import torch
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader
import onnx
import onnxruntime as ort

from sample_dataset import KmerFeaturesDataset
from sample_dataset import load_features
from sample_dataset import KMER_BASE_FEATURES

from model_cnn import DNAModNet

def export_onnx_model(model_path, onnx_model_path, kmer):
    print('export to onnx model')
    cpu_device = torch.device("cpu")
    model = DNAModNet(kmer)
    model.load_state_dict(torch.load(model_path, map_location=cpu_device, weights_only=True))
    print(model)
    model.eval()  
    model.to(cpu_device)  

    dummy_input = torch.randn(1, kmer, KMER_BASE_FEATURES)
    dummy_input = dummy_input.to(cpu_device)

    torch.onnx.export(
        model,
        dummy_input,
        onnx_model_path,
        opset_version = 11,
        input_names=['input'],
        output_names=['output'],
        dynamic_axes={
            'input': {0: 'batch_size'},
            'output': {0: 'batch_size'}
        }
    )
    print('Done')

    print('validate onnx model')

    onnx_model = onnx.load(onnx_model_path)
    onnx.checker.check_model(onnx_model)

    ort_session = ort.InferenceSession(onnx_model_path)
    ort_inputs = {ort_session.get_inputs()[0].name: dummy_input.numpy()}
    ort_outs = ort_session.run(None, ort_inputs)

    torch_out = model(dummy_input)

    np.testing.assert_allclose(torch_out.detach().numpy(), ort_outs[0], rtol = 1e-3, atol=1e-5)

    print('validate pass')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('USAGE:')
        print(f"  {sys.argv[0]} model output-model-dir kmer")
        sys.exit (1)
    model_path = sys.argv[1]
    output_model = sys.argv[2]
    kmer = int(sys.argv[3])
    export_onnx_model(model_path, output_model, kmer)
