#!/bin/bash

ccs_train_bam=/data1/chenying/bs3/hg002-train-grch38.bam
ccs_train_5mc=/data1/chenying/bs3/5mc-call-hg002-train

ccs_test_bam=/data1/chenying/bs3/hg002-test-grch38.bam
ccs_test_5mc=/data1/chenying/bs3/5mc-call-hg002-test

pos_bam=/data1/chenying/cpg/bams/pos.hifi.bam
neg_bam=/data1/chenying/cpg/bams/neg.hifi.bam

bs_bed=/data1/chenying/cpg/CpG.gz.bismark.zero.FOR_REV_MERGED.bed

gcn5mc=/data1/chenying/bs3/gcn5mc

pos_train_reads=60000
neg_train_reads=60000

pos_test_reads=5000
neg_test_reads=5000

KMER=400

pos_prob=0.9
neg_prob=0.1
bs_cov=50
bs_pos_prob=0.95
bs_neg_prob=0.05

output=features-k${KMER}

mkdir -p ${output}

function split_5mc_call()
{
	wrk=$1
	pos_call=${wrk}/k20-pos-5mc.txt
	neg_call=${wrk}/k20-neg-5mc.txt
	amb_call=${wrk}/k20-amb-5mc.txt
	all_call=${wrk}/5mc-call.txt

	rm -f ${pos_call}
	rm -f ${neg_call}
	rm -f ${amb_call}

	awk -v pp=${pos_prob} -v np=${neg_prob} -v pc=${pos_call} -v nc=${neg_call} -v ac=${amb_call} \
		'NF==5 {if ($5<=np) {print $0 >> nc}
		        else if ($5>=pp) {print $0 >> pc}}
		NF>5 {if ($8<=np) {print $0 >> nc}
		      else if ($8>=pp) {print $0 >> pc}
		      else {print >> ac}}' ${all_call}

  	amb_call_bed=${wrk}/amb_call.bed
	awk '{print $5"\t"$6"\t"$6+1"\t"$0}' ${amb_call} > ${amb_call_bed}
	rm -f ${amb_call}

	bs_pos_bed=${wrk}/bs_pos.bed
	awk -v mc=${bs_cov} -v pp=${bs_pos_prob} '$5+$6>=mc && $4>=pp' ${bs_bed} > ${bs_pos_bed}
	amb_bs_pos_bed=${wrk}/amb_bs_pos.bed
	bedtools intersect -wo -a ${amb_call_bed} -b ${bs_pos_bed} > ${amb_bs_pos_bed}
	bs_pos_call=${wrk}/bs-pos-5mc.txt
	awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' ${amb_bs_pos_bed} > ${bs_pos_call}
	rm -f ${bs_pos_bed}
	rm -f ${amb_bs_pos_bed}

	bs_neg_bed=${wrk}/bs_neg.bed
	awk -v mc=${bs_cov} -v np=${bs_neg_prob} '$5+$6>=mc && $4<=np' ${bs_bed} > ${bs_neg_bed}
	amb_bs_neg_bed=${wrk}/amb_bs_neg.bed
	bedtools intersect -wo -a ${amb_call_bed} -b ${bs_neg_bed} > ${amb_bs_neg_bed}
	bs_neg_call=${wrk}/bs-neg-5mc.txt
	awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' ${amb_bs_neg_bed} > ${bs_neg_call}
	rm -f ${bs_neg_bed}
	rm -f ${amb_bs_neg_bed}

	rm -f ${amb_call_bed}

	k20_neg_call=${wrk}/k20-neg-5mc.txt
	k20_pos_call=${wrk}/k20-pos-5mc.txt
	NL=$(wc -l ${k20_neg_call} | awk '{print $1}')
	echo ${NL}
	head  -n ${NL} ${k20_pos_call} > ${wrk}/k20-pos-5mc.txt.tmp
	mv ${wrk}/k20-pos-5mc.txt ${wrk}/all-k20-pos-5mc.txt
	mv ${wrk}/k20-pos-5mc.txt.tmp ${wrk}/k20-pos-5mc.txt

	NL=$(wc -l ${bs_neg_call} | awk '{print $1}')
	echo ${NL}
	head -n ${NL} ${bs_pos_call} > ${wrk}/bs-pos-5mc.txt.tmp
	mv ${wrk}/bs-pos-5mc.txt ${wrk}/all-bs-pos-5mc.txt
	mv ${wrk}/bs-pos-5mc.txt.tmp ${wrk}/bs-pos-5mc.txt
}

function extract_ccs_features()
{
	wrk=$1
	bam=$2

	cmd="samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER} -b 1 -t 8 -e ${wrk}/k20-pos-5mc.txt - ${wrk}/k20-pos.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features  -k ${KMER} -b 1 -t 8 -e ${wrk}/k20-pos-5mc.txt - ${wrk}/k20-pos.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER} -b 0 -t 8 -e ${wrk}/k20-neg-5mc.txt - ${wrk}/k20-neg.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER}  -b 0 -t 8 -e ${wrk}/k20-neg-5mc.txt - ${wrk}/k20-neg.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER} -b 1 -t 8 -e ${wrk}/bs-pos-5mc.txt - ${wrk}/bs-pos.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER} -b 1 -t 8 -e ${wrk}/bs-pos-5mc.txt - ${wrk}/bs-pos.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER} -b 0 -t 8 -e ${wrk}/bs-neg-5mc.txt - ${wrk}/bs-neg.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -k ${KMER} -b 0 -t 8 -e ${wrk}/bs-neg-5mc.txt - ${wrk}/bs-neg.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${wrk}/k20-pos.bin ${wrk}/k20-neg.bin ${wrk}/k20.bin"
	echo running command
	echo "  ${cmd}"
	${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${wrk}/k20-pos.bin ${wrk}/k20-neg.bin ${wrk}/k20.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${wrk}/bs-pos.bin ${wrk}/bs-neg.bin ${wrk}/bs.bin"
	echo running command
	echo "  ${cmd}"
	${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${wrk}/bs-pos.bin ${wrk}/bs-neg.bin ${wrk}/bs.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${wrk}/k20.bin ${wrk}/bs.bin ${wrk}/features.bin"
	echo running command
	echo "  ${cmd}"
	${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${wrk}/k20.bin ${wrk}/bs.bin ${wrk}/features.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi
}

function extract_pacbio_treated_features()
{
	cmd="samtools view -@8 ${pos_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 1 -k ${KMER} -t 8 -n ${pos_test_reads} - ${output}/pacbio-treated-test.pos.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${pos_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 1 -k ${KMER} -t 8 -n ${pos_test_reads} - ${output}/pacbio-treated-test.pos.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="samtools view -@8 ${pos_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 1 -k ${KMER} -t 8 -s ${pos_test_reads} -n ${pos_train_reads} - ${output}/pacbio-treated-train.pos.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${pos_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 1 -k ${KMER} -t 8 -s ${pos_test_reads} -n ${pos_train_reads} - ${output}/pacbio-treated-train.pos.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="samtools view -@8 ${neg_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 0 -k ${KMER} -t 8 -n ${neg_test_reads} - ${output}/pacbio-treated-test.neg.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${neg_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 0 -k ${KMER} -t 8 -n ${neg_test_reads} - ${output}/pacbio-treated-test.neg.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="samtools view -@8 ${neg_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 0 -k ${KMER} -t 8 -s ${neg_test_reads} -n ${neg_train_reads} - ${output}/pacbio-treated-train.neg.bin"
	echo running command
	echo "  ${cmd}"
	samtools view -@8 ${neg_bam} | ${gcn5mc}/src/cpp_tools/extract_train_features -b 0 -k ${KMER} -t 8 -s ${neg_test_reads} -n ${neg_train_reads} - ${output}/pacbio-treated-train.neg.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${output}/pacbio-treated-train.pos.bin ${output}/pacbio-treated-train.neg.bin ${output}/pacbio-treated-train.bin"
	echo running command
	echo "  ${cmd}"
	${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${output}/pacbio-treated-train.pos.bin ${output}/pacbio-treated-train.neg.bin ${output}/pacbio-treated-train.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi

	cmd="${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${output}/pacbio-treated-test.pos.bin ${output}/pacbio-treated-test.neg.bin ${output}/pacbio-treated-test.bin"
	echo running command
	echo "  ${cmd}"
	${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${output}/pacbio-treated-test.pos.bin ${output}/pacbio-treated-test.neg.bin ${output}/pacbio-treated-test.bin
	if [ $? -ne 0 ]; then
		exit 1
	fi
}

echo step 1: split ${ccs_train_5mc}/5mc-call.txt
split_5mc_call ${ccs_train_5mc}

echo step 2: extract ccs features
extract_ccs_features ${ccs_train_5mc} ${ccs_train_bam}

echo step 3: extract pacbio treated features
extract_pacbio_treated_features

echo step 4: split ${ccs_test_5mc}/5mc-call.txt
split_5mc_call ${ccs_test_5mc}

echo step 5: extract ccs test features
extract_ccs_features ${ccs_test_5mc} ${ccs_test_bam}

echo step 6: merge ${ccs_train_5mc}/features.bin ${output}/pacbio-treated-train.bin ${output}/train.bin
${gcn5mc}/src/cpp_tools/merge_features -k ${KMER} ${ccs_train_5mc}/features.bin ${output}/pacbio-treated-train.bin ${output}/train.bin
if [ $? -ne 0 ]; then
	exit 1
fi

python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER} 1 ${output}/train.bin
if [ $? -ne 0 ]; then
	exit 1
fi

python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER} 1 ${ccs_test_5mc}/features.bin
if [ $? -ne 0 ]; then
	exit 1
fi

python ${gcn5mc}/src/py-mol/make_np_features.py ${KMER} 1 ${output}/pacbio-treated-test.bin
if [ $? -ne 0 ]; then
	exit 1
fi
