# Introduction

This example computes the modification frequencies on the first 100 kbp of the first sequence in reference genome [GCF_000001735.4_TAIR10.1_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz). The correlation between between directory counting results of `hifimeth` and the BS-seq (100X) results is 0.956.

This example contains four steps.

# Step 1: Prepare data

### 1.1 Download read data

```shell
$ wget https://github.com/xiaochuanle/hifimeth/releases/download/Examples/ara-examples.tar.gz
$ tar xzvf ara-examples.tar.gz
$ cd ara-examples/
```
`ara-examples.tar.gz` contains five files:
* `arabidopsis.hifi_kinetics.bam` Read file.
* `Ara.CpG.gz.bismark.cov` Modification frequencies computed with BS-seq (100X)
* `SITE_CORR.py` A python script to compute correlation
* `normalize_bismark_freq_call.py` A python script to normalize `Ory.CpG.gz.bismark.cov`.
* `run_arab_example.sh` Commands to run the whole example.

### 1.2 Download reference genome

```shell
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
$ gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz 
$ mv GCF_000001735.4_TAIR10.1_genomic.fna GCF_000001735.4_TAIR10.1_genomic.fasta
$ pbmm2 index --preset CCS GCF_000001735.4_TAIR10.1_genomic.fasta GCF_000001735.4_TAIR10.1_genomic.fasta.mmi
$ samtools faidx GCF_000001735.4_TAIR10.1_genomic.fasta
```

# Step 2: fill data paths

Edit the shell script `run_arab_example.sh` and fill the following four paths:
```shell
hifimeth=/data1/chenying/hifimeth
reference=/data1/chenying/ara/ara-examples/GCF_000001735.4_TAIR10.1_genomic.fasta
bam=/data1/chenying/ara/ara-examples/arabidopsis.hifi_kinetics.bam
bismark=/data1/chenying/ara/ara-examples/Ara.CpG.gz.bismark.cov
```

# Step 3: Run the example
``` shell
$ ./run_arab_example.sh
```

# Step 4: Catch the results

After successfully runing step 3, we can see the correlations finally:
``` shell
bismark.normalized.tsv	hifimeth.count.tsv	0.956113042156851
bismark.normalized.tsv	hifimeth.model.tsv	0.9622230657523669
```

