# Introduction

This example computes the modification frequencies on the first 100 kbp of the first sequence in reference genome [GCF_001433935.1_IRGSP-1.0_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz). The correlation between between directory counting results of `hifimeth` and the BS-seq (100X) results is 0.987.

This example contains four steps.

# Step 1: Prepare data

### 1.1 Download read data

```shell
$ wget https://github.com/xiaochuanle/hifimeth/releases/download/Examples/oryza-examples.tar.gz
$ tar xzvf oryza-examples.tar.gz
$ cd oryza-examples/
```
`oryza-examples.tar.gz` contains five files:
* `oryza.hifi_kinetics.bam` Read file.
* `Ory.CpG.gz.bismark.cov` Modification frequencies computed with BS-seq (100X)
* `SITE_CORR.py` A python script to compute correlation
* `normalize_bismark_freq_call.py` A python script to normalize `Ory.CpG.gz.bismark.cov`.
* `run_oryza_example.sh` Commands to run the whole example.

### 1.2 Download reference genome
```
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
$ gunzip GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
$ mv GCF_001433935.1_IRGSP-1.0_genomic.fna GCF_001433935.1_IRGSP-1.0_genomic.fasta
$ pbmm2 index --preset CCS GCF_001433935.1_IRGSP-1.0_genomic.fasta GCF_001433935.1_IRGSP-1.0_genomic.fasta.mmi
$ samtools faidx GCF_001433935.1_IRGSP-1.0_genomic.fasta
```

# Step 2: fill data paths

Edit the shell script `run_oryza_example.sh` and fill the following four paths:
```shell
hifimeth=/data1/chenying/hifimeth
reference=/data1/chenying/6ma/oryza-examples/GCF_001433935.1_IRGSP-1.0_genomic.fasta
bam=/data1/chenying/6ma/oryza-examples/oryza.hifi_kinetics.bam
bismark=/data1/chenying/6ma/oryza-examples/Ory.CpG.gz.bismark.cov
```

# Step 3: Run the example
``` shell
$ ./run_oryza_example.sh
```

# Step 4: Catch the results

After successfully runing step 3, we can see the correlations finally:
``` shell
bismark.normalized.tsv	hifimeth.count.tsv	0.9872284419143924
bismark.normalized.tsv	hifimeth.model.tsv	0.9875577156519046
```