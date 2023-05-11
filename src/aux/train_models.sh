#!/bin/bash

./make_train_data_k20.sh

python gcn5mc/src/py-mol/train.py --model_dir models --train_file features_k20/p2-train.bin.npy --test_file features_k20/p2-test.bin.npy --kmer_size 20

# hg002-train.bam contains the first 500,000 reads in m64008_201124_002822.bam
# hg002-test.bam contains the next 10,000 reads in m64008_201124_002822.bam
pbmm2 align --preset CCS -N 1 /data1/chenying/cpg/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa hg002-train.bam hg002-train-grch38.bam
pbmm2 align --preset CCS -N 1 /data1/chenying/cpg/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa hg002-test.bam hg002-test-grch38.bam

./gcn5mc/src/scripts/gcn5mc.sh --src /data1/chenying/bs3/gcn5mc \
    --num_threads 48 \
    --model_path models/kmer_20_epoch_2.ckpt \
    --reference /data1/chenying/cpg/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa \
    --out 5mc-call-hg002-train \
    --input hg002-train-grch38.bam \
    --kmer_size 20

./gcn5mc/src/scripts/gcn5mc.sh --src /data1/chenying/bs3/gcn5mc \
    --num_threads 48 \
    --model_path models/kmer_20_epoch_2.ckpt \
    --reference /data1/chenying/cpg/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa \
    --out 5mc-call-hg002-test \
    --input hg002-test-grch38.bam \
    --kmer_size 20

./extract-features-k400.sh

python gcn5mc/src/py-mol/train.py --pretrained_model models/kmer_20_epoch_2.ckpt \
    --train_file features-k400/train.bin.npy \
    --test_file 5mc-call-hg002-test/features.bin.npy \
    --test_file2 features-k400/pacbio-treated-test.bin.npy \
    --model_dir models
