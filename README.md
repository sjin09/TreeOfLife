# Code accompanying *Germline and somatic mutational processes across the Tree of Life*

Preprint available is at bioRxiv, .

This repository contains all the code and usage instructions required to perform the analysis described in the manuscript. Specifically, the repository provides code for the following:

- Counting the number of trinucleotides (3 mers) where the middle base is a pyrimidine base (cytosine and thymine) from a reference FASTA file.
- Generating SBS52 and SBS96 counts and plots from VCF files containing germline and somatic mutations, respectively.
    1. Himut calculates the expected number of somatic mutations based on the callable positions in the reference genome and the callable bases from Pacific Biosciences CCS reads. Additionally, himut generates a bar plot of the expected number of somatic mutations following the SBS96 classification system. Please refer to the methods section of the manuscript for a detailed description.
- Generating normalised SBS52 and SBS96 counts and plots, ensuring each trinucleotide contributes an equal proportion to germline and somatic mutations. 
- Transforming SBS96 counts into SBS52 counts for the comparison of germline and somatic mutational spectra and mutational signatures.
- R code to perform somatic mutational signature extraction using HDP, which requires a matrix where the rows are (normalised) SBS96 counts and the columns are SBS96 classification.
- R code to perform germline mutational signature extraction using HDP, which requires a matrix where the rows are (normalised) SBS52 counts and the columns are SBS 52 classification.
- Phylogenetic analysis of germline and somatic mutational signatures.

Python scripts and R code are used downstream of germline and somatic mutation detection using [deepvariant](https://github.com/google/deepvariant) and [himut](https://github.com/sjin09/himut), respectively. Please follow the user's guide in the relevant repository for detailed instructions on usage and implementation.

## User's Guide

### Trinucletide counts

```
python get_reference_tricounts.py -i sample.fasta --target sample.target -o sample.tri
```

### SBS52 counts

```
```

### Trinucleotide normalised SBS52 counts

```
```

### SBS96 counts

Please note that the script here is used to retrieve and plot raw SBS96 counts and not the expected number of somatic mutations based on the callalbe positions in the reference genome and the callable bases from Pacific Biosciences CCS reads.

```
python get_sbs96_counts.py -i sample.himut.vcf.bgz --ref-fasta sample.fasta -o sample.himut.sbs96.tsv
python get_sbs96_barplot.py -i sample.himut.sbs96.tsv --sample sample -o sample.himut.sbs96.pdf
```

### Trinucleotide normalised SBS96 counts

```
```

### R code for germline mutational signature extraction

```
```

### R code for somatic mutational signature extraction

```
```

### Julia script for phylogenetic signal analysis

```
```
