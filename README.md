# Introduction

GCN5mc is a deep neural network model for calling 5mc for PacBio HiFi reads.

1. At single-molecule resolution, GCN5mc achieves F1-score over 94% and AUC over 99%
2. Calling modification frequencies genomoe, GCN5mc achieves correlation with BS-seq results over 0.9 with counting mode.
3. GCN5mc is built upon the [graph neural network operator](https://arxiv.org/abs/1810.02244) that implemented in the [PyG](https://pytorch-geometric.readthedocs.io/en/latest/index.html) library.
4. A GPU accelerator is required to run GNC5mc

## Installation

### Step 1: Install the Pytorch and PyG environment

Install Pytorch with the default conda channel:
```shell
conda create -n gcn5mc python=3.10
conda activate gcn5mc
conda install cudatoolkit=11.7 -c nvidia
pip install --force-reinstall torch=1.13.1 torchvision=0.14.1 torchaudio=0.13.1 wget https://data.pyg.org/whl/torch-1.13.0+cu117.html
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv torch_geometric -f torch-1.13.0+cu117.html
```
Or with the `tuna.tsinghua.edu.cn` channel:
```shell
conda create -n gcn5mc python=3.10 
conda activate gcn5mc
conda install cudatoolkit=11.7 -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
pip install --force-reinstall torch=1.13.1 torchvision=0.14.1 torchaudio=0.13.1 wget https://data.pyg.org/whl/torch-1.13.0+cu117.html
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv torch_geometric -f torch-1.13.0+cu117.html
```

### Step 2: Install `GCN5mc`

``` shell
git clone
cd gcn5mc
./install_cpp_tools.sh
cd ..
```

## Usage

### Prepare data

The `gcn5mc` is downloaded in directory
```shell
/data1/chenying/bs3/gcn5mc/
```

The reads are found in the source of `gcn5mc`:
```shell
/data1/chenying/bs3/gcn5mc/examples/chr1-reads.bam
```

Download the `GRCh38` reference (`pbmm2` does not recognise the `.fna` format):
```shell
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
$ gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
$ mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
$ pbmm2 index --preset CCS GCA_000001405.15_GRCh38_no_alt_analysis_set.fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.mmi
$ samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
```

Download a phased VCF:
```shell
$ wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/SupplementaryFiles/HG002_NA24385_GRCh38_1_22_v4.2_benchmark_phased_StrandSeq_whatshap.vcf.gz
$ tabix -p vcf HG002_NA24385_GRCh38_1_22_v4.2_benchmark_phased_StrandSeq_whatshap.vcf.gz
```

There are three using cases.

### Case 1: Calling 5mc at single-molecule resolution only.

```shell
$ /data1/chenying/bs3/gcn5mc/src/scripts/gcn5mc.sh \
	--src /data1/chenying/bs3/gcn5mc \
    --num_threads 48 \
    --model_path /data1/chenying/bs3/gcn5mc/models/gcn5mc-mol-k400.ckpt \
    --input /data1/chenying/bs3/gcn5mc/examples/chr1-reads.bam \
    --out 5mc-call
```

#### Results

```shell
5mc-call/5mc-call.bam
5mc-call/5mc-call.txt
```

`5mc-call/5mc-call.bam` adds calling results to the input bam `/data1/chenying/bs3/gcn5mc/examples/chr1-reads.bam`.  5mC base modification values are read from the MM and ML auxiliary tags which encode base modifications and confidence values. These tags are further described in the [SAM tag specification document](https://samtools.github.io/hts-specs/SAMtags.pdf).

`5mc-call/5mc-call.txt` lists calling results in human readable format:
```shell
m64008_201124_002822/104531456/ccs	8	2	15052	0.9861
m64008_201124_002822/104531456/ccs	8	12	15042	1.0000
m64008_201124_002822/104531456/ccs	8	80	14974	0.0026
m64008_201124_002822/104531456/ccs	8	170	14884	0.9995
m64008_201124_002822/104531456/ccs	8	181	14873	0.9998
m64008_201124_002822/104531456/ccs	8	213	14841	0.9860
m64008_201124_002822/104531456/ccs	8	273	14781	0.9969
m64008_201124_002822/104531456/ccs	8	281	14773	0.0308
m64008_201124_002822/104531456/ccs	8	288	14766	0.0100
m64008_201124_002822/104531456/ccs	8	311	14743	0.9862
```
Each calling result ocupies one line and contains five colummns:
1. Read name
2. Numeric id of this read
3. The position of this CpG in the forward strand of the read
4. The position of this CpG in the reverse strand of the read
5. Methylation probability.

### Case 2: Calling 5mc at single-molecule resolution and calling modification frequencies both

Step 1: map the bam to a reference:
``` shell
pbmm2 align --preset CCS --sort GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
	/data1/chenying/bs3/gcn5mc/examples/chr1-reads.bam chr1-reads-grch38.bam
```

Step 2: call 5mc with the mapped bam:
```shell
/data1/chenying/bs3/gcn5mc/src/scripts/gcn5mc.sh \
	--src /data1/chenying/bs3/gcn5mc \
    --model_path /data1/chenying/bs3/gcn5mc/models/gcn5mc-mol-k400.ckpt \
    --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --input chr1-reads-grch38.bam \
    --out 5mc-call-grch38
```
#### Results

```shell
5mc-call-grch38/5mc-call.bam
5mc-call-grch38/5mc-call.txt
```
Details of `5mc-call-grch38/5mc-call.bam` are found in Case 1.

`5mc-call-grch38/5mc-call.txt` lists calling results in human readable format:
```shell
m64008_201124_002822/102303964/ccs      0       5069    12204   0.2486
m64008_201124_002822/102303964/ccs      0       5079    12194   0.7334
m64008_201124_002822/102303964/ccs      0       5097    12176   0.9999
m64008_201124_002822/102303964/ccs      0       5129    12144   chr1    1955103 +       1.0000
m64008_201124_002822/102303964/ccs      0       5140    12133   chr1    1955114 +       1.0000
m64008_201124_002822/102303964/ccs      0       5142    12131   chr1    1955116 +       1.0000
m64008_201124_002822/102303964/ccs      0       5157    12116   chr1    1955131 +       1.0000
m64008_201124_002822/102303964/ccs      0       5161    12112   chr1    1955135 +       1.0000
```
Each calling result ocupies one line and contains five or eight colummns.
Details of a calling result containing five colummns are given in Case 1. In this situation, it means that this CpG site fail to be mapped to a CpG site on the reference.

A calling result containing eight colummns also contains mapping info:
1. Read name
2. Numeric id of this read
3. The position of this CpG in the forward strand of the read
4. The position of this CpG in the reverse strand of the read
5. The chromosome name this read mapped to
6. The chromosome position this CpG mapped to
7. The strand (`+` for forward and `-` for reverse) of this read mapped to the reference
8. Methylation probability

Step 3: call modification frequencies

``` shell
/data1/chenying/bs3/gcn5mc/src/scripts/gcn5mc-freq.sh \
	--src /data1/chenying/bs3/gcn5mc \
    --model_path /data1/chenying/bs3/gcn5mc/models/gcn5mc-freq-s11.ckpt \
    --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --num_threads 48 \
    --bam 5mc-call-grch38/5mc-call.bam \
    --out freq-call-grch38
```

#### results

``` shell
freq-call-grch38/freq-call.count.txt
freq-call-grch38/freq-call.model.txt
```

`freq-call-grch38/freq-call.count.txt` contains calling results computed with *counting mode*.
`freq-call-grch38/freq-call.model.txt` contains calling results computed with *model mode*.
The two files share the same format as BS-seq (it is actually a bed file):
```shell
chr1	1950011	1950012	1.0000	1	0
chr1	1950075	1950076	1.0000	1	0
chr1	1950137	1950138	1.0000	1	0
chr1	1950262	1950263	1.0000	1	0
chr1	1950276	1950277	1.0000	1	0
```

Using `freq-call-grch38/freq-call.count.txt` for subsequent analysis is recommanded.

### Case 3: Calling 5mc at single-molecule resolution and calling haploid and diploid modification frequencies both

Step 1: map the bam to a reference:
``` shell
pbmm2 align --preset CCS --sort GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
	/data1/chenying/bs3/gcn5mc/examples/chr1-reads.bam chr1-reads-grch38.bam
```

Step 2: call 5mc with the mapped bam:
```shell
/data1/chenying/bs3/gcn5mc/src/scripts/gcn5mc.sh \
	--src /data1/chenying/bs3/gcn5mc \
    --model_path /data1/chenying/bs3/gcn5mc/models/gcn5mc-mol-k400.ckpt \
    --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --input chr1-reads-grch38.bam \
    --out 5mc-call-grch38
```
#### Results

Same as Case 2.

Step 3: add haplotypes to the bam:

``` shell
$ samtools index -@8 5mc-call-grch38/5mc-call.bam
$ whatshap haplotag --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
	--ignore-read-groups \
    --o 5mc-call-grch38/5mc-call-tag.bam \
    HG002_NA24385_GRCh38_1_22_v4.2_benchmark_phased_StrandSeq_whatshap.vcf.gz \
    5mc-call-grch38/5mc-call.bam
```

Step 4: calling modification frequencies on the taged bam:

```shell
/data1/chenying/bs3/gcn5mc/src/scripts/gcn5mc-freq.sh \
	--src /data1/chenying/bs3/gcn5mc \
    --model_path /data1/chenying/bs3/gcn5mc/models/gcn5mc-freq-s11.ckpt \
    --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --num_threads 48 --bam 5mc-call-grch38/5mc-call-tag.bam \
    --out freq-call-grch38-tag \
    --phase
```

#### results

``` shell
freq-call-grch38-tag/freq-call.count.txt
freq-call-grch38-tag/freq-call.model.txt
freq-call-grch38-tag/freq-call.hp1.count.txt
freq-call-grch38-tag/freq-call.hp1.model.txt
freq-call-grch38-tag/freq-call.hp2.count.txt
freq-call-grch38-tag/freq-call.hp2.model.txt
```