# Code accompanying *Germline and somatic mutational processes across the tree of life*

Preprint available is at bioRxiv, .

This repository contains all the code and usage instructions required to perform the analysis described in the manuscript. Specifically, the repository provides code for the following:

- Counting the number of trinucleotides (3 mers) where the middle base is a pyrimidine base (cytosine and thymine) from a reference FASTA file.
- Generating SBS52 and SBS96 counts and plots from VCF files containing germline and somatic mutations, respectively.
    1. Himut calculates the observed number of somatic mutations based on the callable positions in the reference genome and the callable bases from Pacific Biosciences CCS reads. Additionally, himut generates a bar plot of the observed number of somatic mutations following the SBS96 classification system. Please refer to the methods section of the manuscript for a detailed description.
- Generating normalised SBS52 and SBS96 counts and plots, ensuring each trinucleotide contributes an equal proportion to germline and somatic mutations. 
- Transforming SBS96 counts into SBS52 counts for the comparison of germline and somatic mutational spectra and mutational signatures.
- R code to perform somatic mutational signature extraction using HDP, which requires a matrix where the rows are (normalised) SBS96 counts and the columns are SBS96 classification.
- R code to perform germline mutational signature extraction using HDP, which requires a matrix where the rows are (normalised) SBS52 counts and the columns are SBS 52 classification.
- Phylogenetic analysis of germline and somatic mutational signatures.

Python scripts and R code are used downstream of germline and somatic mutation detection using [deepvariant](https://github.com/google/deepvariant) and [himut](https://github.com/sjin09/himut), respectively. Please follow the user's guide in the relevant repository for detailed instructions on usage and implementation.

## User's Guide

### Create a conda environment install dependencies to run python scripts

```
conda env create --file environment.yml  # create conda environment
conda activate treeoflife
poetry install  # install the package to the virtual environment
```

### Reference FASTA file trinucleotide (3-mer) counts

There are 64 possible trinucleotides (4 ** 3). The script below counts the occurrences of all 64 trinucleotides. However, if the middle base is a purine (adenine or guanine), the count is incremented using the reverse complement of the trinucleotide. Hence, the counts of trinucleotides where the middle base is a pyrimidine are returned. 

```
python get_reference_tricounts.py \
    -i sample.fasta \
    --target sample.target \
    -o sample.tri
```

### SBS52 classification counts

Unlike somatic mutations, germline mutations are classified into 52 categories (SBS52 classification). Please refer to the methods section of the manuscript for a detailed description of how the SBS96 classification is collapsed into the SBS52 classification.

```
## The file sample.vcf.gz contains germline mutations.

python get_sbs52_counts.py \
    -i sample.vcf.bgz \
    --ref-fasta sample.fasta \
    -o sample.sbs52.tsv

python get_sbs52_barplot.py \
    -i sample.sbs52.tsv \
    --sample sample \
    -o sample.sbs52.pdf
```

#### SBS52 counts, where each trinucleotide contributes equally.

```
python get_tri_equal_weight_sbs52_counts.py \
    -i sample.vcf.bgz \
    --ref-fasta sample.fasta \
    --target sample.target \
    --tri sample.tri \
    -o sample.tri_equal_weight.sbs52.tsv \
    --is-sample-reference-sample

python get_sbs52_barplot.py \
    -i sample.tri_equal_weight.sbs52.tsv \
    --sample sample \
    -o sample.tri_equal_weight.sbs52.pdf
```

### SBS96 classification counts

Somatic mutations can be classified into 96 categories (SBS96 classification), depending on the 6 different classes of base substitution and 16 combinations of the bases immediately 5’ and 3’ to the mutation. 

Please note that the script here is used to retrieve and plot raw SBS96 counts and not the observed number of somatic mutations based on the callable positions in the reference genome and the callable bases from Pacific Biosciences CCS reads.

```
## The file sample.vcf.gz contains somatic mutations.

python get_sbs96_counts.py \
    -i sample.vcf.bgz \
    --ref-fasta sample.fasta \
    -o sample.sbs96.tsv

python get_sbs96_barplot.py \
    -i sample.sbs96.tsv \
    --sample sample \
    -o sample.sbs96.pdf
```

#### SBS96 counts, where each trinucleotide contributes equally.

```
python get_tri_equal_weight_sbs96_counts.py \
    -i sample.sbs96.tsv \
    --tri sample.tri \
    -o sample.tri_equal_weight.sbs96.tsv
```

#### Collapse SBS96 classification counts to SBS52 classification counts

Please refer to the methods section of the manuscript or the `get_sbs96_to_sbs52_lookup_table.py` script for a detailed description of how the SBS96 classification is mapped to the SBS52 classification.

```
## sbs96_to_sbs52_lookup_table.tsv can be found under the scripts directory
## get_sbs96_to_sbs52_lookup_table.py is used to generate the sbs96_to_sbs52_lookup_table.tsv

python sbs96_to_sbs52.py \
    -i sample.sbs96.tsv \
    --sbs96-to-sbs52 sbs96_to_sbs52.tsv \
    -o sample.sbs96_to_sbs52.tsv

python sbs96_to_sbs52.py \
    -i sample.tri_equal_weight.sbs96.tsv \
    --sbs96-to-sbs52 sbs96_to_sbs52.tsv \
    -o sample.tri_equal_weight.sbs96_to_sbs52.tsv
```
### Calculate mutation burden per cell

Please note that the calculation here does not account for the number of mutations attributed to technical artefacts, which must be considered to obtain the correct mutation burden per cell for the sample. The number of mutations attributed to technical artefacts can be obtained through either mutational signature extraction or mutational signature attribution.

```
## Himut normcounts returns sample.observed_sbs96.tsv
## Here, sample.target needs to have both the autosomes and sex chromosomes to 
## calculate the mutation burden per genome and the mutation burden per cell.

python get_mutation_burden_per_cell.py \
    -i sample.observed_sbs96.tsv \
    --ref-fasta sample.fasta \
    --target sample.target \
    --ploidy 2 \
    -o sample.burden
```

- Avian sex chromosomes: chrW and chrZ
- Mammalian sex chromosomes: chrX and chrY
- Sex determination in insects of the Hymenoptera order is based on ploidy. Males possess a haploid genome, while females have a diploid genome.

### R code for mutational signature extraction

In the file `hdp_input.mat`, the rows represent samples, and the columns correspond to either the SBS52 or the SBS96 classification. Each element in the matrix indicates the counts of mutations.

```
## Germline mutational signature extraction

Rscript hdp_noprior_SBS52.R ${hdp_input.mat} ${chain_index} ${hdp_output_prefix} # repeat this ten times with chain index=1-10
Rscript hdp_extraction_SBS52.R ${hdp_output_prefix} ${hdp_input.mat} ${output_directory} ${output_prefix}

## Somatic mutational signature extraction
Rscript hdp_noprior_SBS96.R ${hdp_input.mat} ${chain_index} ${hdp_output_prefix} # repeat this ten times with chain index_1=10
Rscript hdp_extraction_SBS96.R ${hdp_output_prefix} ${hdp_input.mat} ${output_directory} ${output_prefix}
```
`hdp_extraction_SBS52/SBS96.R` returns two excel spreadsheets `${output_prefix}_HDP_sigs.csv` and `${output_prefix}_HDP_exposure.csv`. 
- In the file `${output_prefix}_HDP_sigs.csv`, the rows represent SBS52/SBS96 classifications, and the columns correspond to the mutational signatures. Each element in the matrix indicates the probability of a mutation occuring in a specific sequence context under a given mutational signature.
- In the file `${output_prefix}_HDP_exposure.csv`, the rows represent samples, and the columns correspond to the mutational signatures. Each element in the matrix reflects the contribution of each mutational signature to the sample's mutation burden.

### Calculate Abouheif's Cmean for each germline and somatic mutational signature

Given a table containing the taxonomic classification of each species and another table with their respective mutational signature exposures, the script `get_abouheif_C_mean.py` constructs a phylogenetic tree based on these classifications. It then calculates the phylogenetic proximity between species and assesses the phylogenetic signal of each mutational signature.


```
### Phylogenetic signal analysis of germline mutational signatures

python get_abouheif_C_mean.py \
    --taxonomic-classification dtol_sample_taxonomic_classification.csv \
    --signature-exposures gtol_exposure.csv \
    --pdf gtol_abouheif_cmean.pdf \
    -o gtol_abouheif_cmean.csv
```

```
### Phylogenetic signal analysis of somatic mutational signatures

python get_abouheif_C_mean.py \
    --taxonomic-classification dtol_sample_taxonomic_classification.csv \
    --signature-exposures stol_exposure.csv \
    --pdf stol_abouheif_cmean.pdf \
    -o stol_abouheif_cmean.csv \
    --is-somatic
```
