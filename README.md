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

### Reference FASTA file trinucletide (3-mer) counts

```
python get_reference_tricounts.py -i sample.fasta --target sample.target -o sample.tri
```

### SBS52 classification counts

Unlike somatic mutations, germline mutations are classified into 52 categories (SBS52 classification). Please refer to the methods section of the manuscript for a detailed description of how the SBS96 classification is collapsed into the SBS52 classification.

```
## The file sample.vcf.gz contains germline mutations.
python get_sbs52_counts.py -i sample.vcf.bgz --ref-fasta sample.fasta -o sample.sbs52.tsv
python get_sbs52_barplot.py -i sample.sbs52.tsv --sample sample -o sample.sbs52.pdf
```

#### SBS52 counts, where each trinucleotide contributes equally.


```
```

### SBS96 classification counts

Somatic mutations can be classified into 96 categories (SBS96 classification), depending on the 6 different classes of base substitution and 16 combinations of the bases immediately 5’ and 3’ to the mutation. 

Please note that the script here is used to retrieve and plot raw SBS96 counts and not the expected number of somatic mutations based on the callalbe positions in the reference genome and the callable bases from Pacific Biosciences CCS reads.

```
## The file sample.vcf.gz contains somatic mutations.
python get_sbs96_counts.py -i sample.vcf.bgz --ref-fasta sample.fasta -o sample.sbs96.tsv
python get_sbs96_barplot.py -i sample.sbs96.tsv --sample sample -o sample.sbs96.pdf
```

#### SBS96 counts, where each trinucleotide contributes equally.

```
```

#### Collapse SBS96 classification counts to SBS52 classification counts

```
## sbs96_to_sbs52.tsv can be found under the scripts directory
python sbs96_to_sbs52.py -i sample.sbs96.tsv --sbs96-to-sbs52 sbs96_to_sbs52.tsv -o sample.sbs96_to_sbs52.tsv
python sbs96_to_sbs52.py -i sample.tri_equal_weights.sbs96.tsv --sbs96-to-sbs52 sbs96_to_sbs52.tsv -o sample.tri_equal_weights.sbs96_to_sbs52.tsv
```

### R code for mutational signature extraction

#### Germline mutational signature extraction
```
```

#### Somatic mutational signature extraction

```
```

### Julia script for phylogenetic signal analysis

```
```
