#!/bin/bash

KMER_SIZE=20

features_dir=features_k${KMER_SIZE}

mkdir -p ${features_dir}

gcn5mc=/data1/chenying/bs3/gcn5mc
export PATH=$PATH:${gcn5mc}/src/cpp_tools

pos_bam=/data1/chenying/cpg/bams/pos.hifi.bam
neg_bam=/data1/chenying/cpg/bams/neg.hifi.bam

p2_test_reads_pos=5000
p2_train_reads_pos=400000
p2_test_reads_neg=5000
p2_train_reads_neg=400000


cmd="samtools view -@8 ${pos_bam} | extract_train_features -b 1 -k ${KMER_SIZE} -t 16 -n ${p2_test_reads_pos} - ${features_dir}/p2-test.pos.bin"
echo "${cmd}"
samtools view -@8 ${pos_bam} | extract_train_features -b 1 -k ${KMER_SIZE} -t 16 -n ${p2_test_reads_pos} - ${features_dir}/p2-test.pos.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="samtools view -@8 ${pos_bam} | extract_train_features -b 1 -k ${KMER_SIZE} -t 16 -s ${p2_test_reads_pos} -n ${p2_train_reads_pos} - ${features_dir}/p2-train.pos.bin"
echo "${cmd}"
samtools view -@8 ${pos_bam} | extract_train_features -b 1 -k ${KMER_SIZE} -t 16 -s ${p2_test_reads_pos} -n ${p2_train_reads_pos} - ${features_dir}/p2-train.pos.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="samtools view -@8 ${neg_bam} | extract_train_features -b 0 -k ${KMER_SIZE} -t 16 -n ${p2_test_reads_neg} - ${features_dir}/p2-test.neg.bin"
echo "${cmd}"
samtools view -@8 ${neg_bam} | extract_train_features -b 0 -k ${KMER_SIZE} -t 16 -n ${p2_test_reads_neg} - ${features_dir}/p2-test.neg.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="samtools view -@8 ${neg_bam} | extract_train_features -b 0 -k ${KMER_SIZE} -t 16 -s ${p2_test_reads_neg} -n ${p2_train_reads_neg} - ${features_dir}/p2-train.neg.bin"
echo "${cmd}"
samtools view -@8 ${neg_bam} | extract_train_features -b 0 -k ${KMER_SIZE} -t 16 -s ${p2_test_reads_neg} -n ${p2_train_reads_neg} - ${features_dir}/p2-train.neg.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="merge_features -k ${KMER_SIZE} ${features_dir}/p2-train.pos.bin ${features_dir}/p2-train.neg.bin ${features_dir}/p2-train.bin"
echo "${cmd}"
merge_features -k ${KMER_SIZE} ${features_dir}/p2-train.pos.bin ${features_dir}/p2-train.neg.bin ${features_dir}/p2-train.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="merge_features -k ${KMER_SIZE} ${features_dir}/p2-test.pos.bin ${features_dir}/p2-test.neg.bin ${features_dir}/p2-test.bin"
echo "${cmd}"
merge_features -k ${KMER_SIZE} ${features_dir}/p2-test.pos.bin ${features_dir}/p2-test.neg.bin ${features_dir}/p2-test.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER_SIZE} 1 ${features_dir}/p2-train.bin"
echo "${cmd}"
python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER_SIZE} 1 ${features_dir}/p2-train.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

cmd="python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER_SIZE} 1 ${features_dir}/p2-test.bin"
echo "${cmd}"
python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER_SIZE} 1 ${features_dir}/p2-test.bin
if [ $? -ne 0 ]; then
	echo fail
	exit 1
fi

rm -f ${features_dir}/*.bin
