# -*- coding: utf-8 -*-
import os
import random
import re
import sys
import time
import math
import numpy as np
import torch
from torch_geometric.data import Data, Dataset 
from torch_geometric.loader import DataLoader
from tqdm import tqdm

import kmer_model
from kmer_model import KmerModel

import read_loader
from read_loader import ModReadLoader

def dump_usage(program):
    print('')
    print('USAGE:')
    print('%s read.bam' % (program))

def main(argv):
    n_param = len(argv)
    if n_param < 2:
        dump_usage(argv[0])
        sys.exit(1)
    read_sam_path = argv[1]
    reads = ModReadLoader(read_sam_path)

    batch_i = 0
    while True:
        read_list = reads.fetch_batch(1000)
        if len(read_list) == 0:
            break 
        batch_i += 1
        if batch_i == 10:
            break


if __name__ == "__main__":
    main(sys.argv)