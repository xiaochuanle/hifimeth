# HiFiMeth (v1.1.0)

![alt text](https://img.shields.io/github/v/release/xiaochuanle/hifimeth) ![alt text](https://img.shields.io/badge/Bioinformatics-PacBio%20HiFi-green) ![alt text](https://img.shields.io/badge/License-MIT-yellow.svg)

- [HiFiMeth (v1.1.0)](#hifimeth-v110)
- [Introduction](#introduction)
- [Installation](#installation)
    - [Option 1: Pre-compiled Binary (Most Recommended)](#option-1-pre-compiled-binary-most-recommended)
    - [Option 2: Build from source](#option-2-build-from-source)
    - [Option 3: Install the GPU version](#option-3-install-the-gpu-version)
- [Usage](#usage)
- [Model evaluation](#model-evaluation)

# Introduction

HiFiMeth is a deep-learning-based tool specifically designed for 5mC (CpG/CHG/CHH) methylation detection using PacBio HiFi sequencing data.

🌟 Key Highlights
* The First & Only All-Context Detector: HiFiMeth is the first and currently the only tool capable of detecting 5mC in all sequence contexts (CpG, CHG, and CHH) specifically from PacBio HiFi data.
* Superior Accuracy & Reliability:
  * Extensively Validated: Benchmarked across 11 diverse datasets to ensure robust performance.
  * High Correlation with BS-seq: Demonstrates exceptional consistency with Bisulfite Sequencing (the gold standard):
    * CHG: Genome-wide Pearson correlation reaches 0.900 – 0.973.
    * CHH: Genome-wide Pearson correlation reaches 0.755 – 0.800, setting a new benchmark for non-CpG detection.

* Zero Dependency & Out-of-the-box:
  * Provided as a single pre-compiled binary executable.
  * Zero Python dependency: The inference engine is written entirely in C/C++, eliminating complex environment configurations (no more conda or pip issues).
* CPU-Only High Efficiency:
  * No GPU required: Runs efficiently on standard CPU threads.
  * Lightning-Fast Inference: Processing a 30x Arabidopsis genome (all contexts) takes only ~2 hours using 48 CPU threads.
* Plant-Genome Friendly: Optimized for the complex methylation patterns found in plant genomes, providing high-resolution epigenetic insights where other tools fall short.

HiFiMeth is engineered for high performance and reliability, leveraging the following libraries:
* [htslib (v1.19.1)](https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2): For high-throughput and robust processing of BAM/SAM/VCF genomic data formats.
* [OpenVINO™ (2025.4.0)](https://github.com/openvinotoolkit/openvino/tree/2025.4.0) Toolkit: A high-performance inference engine specifically optimized for Intel/AMD CPUs, enabling lightning-fast deep learning predictions without requiring a GPU.
* [zlib (v1.3.1)](https://zlib.net/fossils/zlib-1.3.1.tar.gz): Used for efficient handling and reading of compressed .gz data streams.

***Note for Users***: Although HiFiMeth is built on these  libraries, you do not need to install them separately. Our pre-compiled binary comes with these dependencies handled (either statically linked or bundled), ensuring a true "plug-and-play" experience.

# Installation

We provide multiple ways to install HiFiMeth.

### Option 1: Pre-compiled Binary (Most Recommended)
This is the fastest way to get started. HiFiMeth provides a standalone binary that is statically or dynamically linked with all necessary dependencies (htslib, OpenVINO, zlib). No root privileges or complex environment setups (like Conda or Python) are required.
```shell
$ wget https://github.com/xiaochuanle/hifimeth/releases/download/v1.1.0/hifimeth_1.1.10_Linux-amd64.tar.bz2
$ tar -jxvf hifimeth_1.1.10_Linux-amd64.tar.bz2
$ cd hifimeth_1.1.10_Linux-amd64/bin
$ echo "export PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc && source ~/.bashrc
```
The last command above adds the folder of `hifimeth` to the system `PATH`, and harcodes this to the shell configuration file `~/.bashrc`.

Now we can call `hifimeth` from any directory:
```shell
hifimeth
```

### Option 2: Build from source

See [build_hifimeth](build_hifimeth.md)

### Option 3: Install the GPU version

For convinient, we also implement a GPU version of HiFiMeth with the LibTorch library. Details are found at [build_hifimeth_gpu](build_hifimeth_gpu.md)

# Usage

To help users rapidly familiarize themselves with HiFiMeth, we have provided a small-scale turorial dataset [P.Patens](https://github.com/xiaochuanle/hifimeth/releases/download/v1.1.0/P.patens.tar.bz2). This data originates from a 5.3 Mbp chromosome of *Physcomitrella patens*. The contents of this compressed archive encompass all inputs referenced in the command-line examples in this section and the subsequent benchmarking analysis.

```shell
$ wget https://github.com/xiaochuanle/hifimeth/releases/download/v1.1.0/P.patens.tar.bz2
$ tar -xjvf P.patens.tar.bz2
```

After decompressing, the directory `P.Patens` contains a shell script `run.sh`, which contains all of the commands listed here and the benchmark below:
```shell
$ cd P.patens
$ ./run.sh
```

We next explain the commnads in `run.sh`. HiFiMeth comprises two core modules: read-level 5mC methylation detection and genome-wide methylation quantification, which are typically executed sequentially in a standard epigenomic studies. Read-level prediction is performed via the `hifimeth call` subcommand:
```shell
$ hifimeth call m84070_250716_151350_s2.bam mod.bam
```
  
In the above command, leveraging its OpenVINO-powered inference engine, HiFiMeth automatically detects and utilizes all available system CPU threads for parallel prediction, eliminating the need for explicit user configuration of thread counts. Furthermore, HiFiMeth automatically identifies and loads the requisite model files from its local directory by default, precluding the need for manual path specification. This high level of automation, encompassing both self-configuring CPU threading and autonomous model loading, significantly streamlines the user experience and enhances usability. By default, HiFiMeth concurrently detects three 5mC contexts (CpG, CHG, and CHH) in one single run. The resulting 5mC modification probabilities are integrated into the output BAM file (`mod.bam` in this example) using standard MM and ML tags. The “Base Modifications” section of the official SAM tags specification (https://samtools.github.io/hts-specs/SAMtags.pdf) provides a detailed description of how modifications are encoded within these two tags, ensuing full compatibility with downstream analysis tools.

To perform genome-wide quantification, the BAM file annotated with MM/ML tags is first mapped back to the reference genome with [pbmm2](https://github.com/PacificBiosciences/pbmm2). This is achieved by first indexing the reference genome and subsequently aligning the modified BAM using the index. This step ensures that the modification probabilities are mapped to specific genomic coordinates, allowing the quantification module to aggregate single-molecule modifications into genomic methylation profiles.
```shell
$ pbmm2 index –preset CCS \
    GCA_000002425.3_Phypa_V5_genomic.fasta \
    GCA_000002425.3_Phypa_V5_genomic.fasta.mmi
$ pbmm2 align –preset CCS –sort -j48 \
    GCA_000002425.3_Phypa_V5_genomic.fasta.mmi \
    mod.bam mod.pbmm2.bam
```

BAM files annotated with standard MM and ML tags, along with detailed alignment information (the CIGAR field), are then fed into the quantification module via the `hifimeth pileup` command:
```shell
$ hifimeth pileup \
    GCA_000002425.3_Phypa_V5_genomic.fasta \
    mod.pbmm2.bam P.patens
```

This command executes two primary tasks. First, it determines an adaptive probability threshold for each methylation context (CpG, CHG, or CHH) by analyzing the global distribution of modification probabilities within the input BAM file (`mod.pbmm2.bam` in this case). This threshold is used to assign a binary methylation state to each HiFi read. Second, the module calculates the site-specific methylation frequency across the genome. Specifically, for each read overlapping a target genomic locus, it is classified as methylated (M) if its modification probability exceeds the adaptive threshold, and unmethylated (U) otherwise. The methylation frequency of that locus is then determined by the formula 100.0×M/(M+U). Finally, the results are exported into three standardized BED-formatted files, categorized by sequence context (CpG, CHG, and CHH). Each BED file contains six columns, detailing the chromosome name, 0-based inclusive start position, exclusive end position, methylation frequency (%), and the respective counts of methylated and unmethylated HiFi reads overlapping the site. Here is a snapshot of `P.patens.CHH.cov.bed` (the prefix `P.patens` is specified by the above command):
```shell
CM009342.2      14986   14987   69.2308 9       4
CM009342.2      14989   14990   7.69231 1       12
CM009342.2      14993   14994   15.3846 2       11
CM009342.2      15006   15007   7.69231 1       12
CM009342.2      15010   15011   0       0       13
```

Notably, for the strand-symmetric CpG and CHG contexts, methylated and unmethylated HiFi read counts from negative-strand guanines are aggregated with those of their corresponding positive-strand cytosines. Consequently, methylation frequencies are calculated and reported exclusively for the cytosine positions on the positive strand.

# Model evaluation

We assess the performance of HiFiMeth across three distinct dimensions using the eleven datasets. Each 5mC context is evaluated independently to facilitate a granular and exhaustive assessment of the model’s performance. In this section, we use the CHH context as a representative example to demonstrate the evaluation workflow. 

First, to benchmark computational efficiency, specifically inference speed and memory usage, we execute the following command:
```
$ hifimeth call -c chh m84070_250716_151350_s2.bam mod.bam
```

This command isolates the detection CHH methylation while bypassing the other two contexts. Upon completion of the call, HiFiMeth automatically logs the wall-clock time and peak memory consumption.

Second, to assess genome-wide quantification performance, we calculate the Pearson correlation coefficient (r) between the HiFiMeth-generated `P.patens.CHH.cov.bed` and the corresponding frequency results derived from BS-seq data. Since the `.cov` files output by Bismark utilize a 1-based coordinate system, we first convert them to the standard 0-based BED format:
```shell
$ hifimeth cov2bed \
    GCA_000002425.3_Phypa_V5_genomic.fasta \
    CHH P.patens.CHH.gz.bismark.cov chh.bed
```

For the CHH context, this conversion primarily involves adjusting the chromosomal coordinates from 1-based to 0-based. For the strand-symmetric CpG and CHG contexts, the counts of methylated and unmethylated reads from reverse-strand guanines are additionally aggregated to their corresponding forward-strand cytosines. We then compute the Pearson correlation (r) between the two BED files as follows:
```
$ hifimeth corr P.patens.CHH.cov.bed chh.bed
```

Finally, we evaluate the performance of HiFiMeth at the single-molecule level. Ground-truth labels are established based on methylation frequencies derived from BS-seq data, where genomic loci with a read depth of at least 5x are selected. Specifically, sites with a methylation frequency of 0% are assigned a binary label of 0 (unmethylated), while those with 100% frequency are assigned a label of 1 (methylated). Using the BAM files containing modification probabilities (MM and ML tags) and alignment metadata (CIGAR), we calculate an adaptive probability threshold as previously described. Read-level sites mapping to the labeled genomic positions are predicted as methylated if their predicted probability is greater than or equal to this threshold, and unmethylated otherwise. 
To ensure a robust and balanced evaluation, we randomly sample 100,000 instances for each ground-truth label, repeating this process five times to generate five independent evaluation subsets:
```shell
$ hifimeth eval \
    GCA_000002425.3_Phypa_V5_genomic.fasta \
    chh.bed mod.pbmm2.bam read-level.eval
```

Performance is assessed through a confusion matrx, where True Positives (TP) and True Negatives (TN) represent correctly predicted methylated and unmethylated sites, respectively, while False Positives (FP) and False Negatives (FN) denote unmethylated sites incorrectly predicted as methylated and methylated sites incorrectly predicted as unmethylated, respectively. Based on these counts, we calculate the Precision (TP/[TP+FP]), Recall (TP/[TP+FN]), and the F1-score (2*Precision * Recall/[Precision + Recall]):
```shell
$ python read_level_eval.py read-level.eval 5
```
  
The python script read_level_eval.py is provided in the torutiral data.
