# Code accompanying *Somatic and germline mutational processes across the Tree of Life*

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
* [Installation instructions](docs/install.md)
* [Count trinucleotides in a reference FASTA file](docs/trinucleotide.md)
* [Counting somatic mutations using the SBS96 classification](docs/sbs96.md)
* [Counting germline mutations using the SBS52 classification](docs/sbs52.md)
* [Map SBS96 counts to SBS52 counts](docs/sbs96_to_sbs52.md)
* [Calculate mutation burden per cell](docs/burden.md)
* [Perform *de novo* mutational signature extraction using HDP](docs/mutational_signature_extraction.md)
* [Phylogenetic signal analysis of mutational signatures](docs/abouheifs_cmean.md)

## Citation
If you use any of the scripts provided here, please cite our bioRxiv preprint.

## Help
If you encounter bugs or have further questions or requests, you can raise an issue at the issue page. I have left academia and may not be able to provide active support, but I will do my best to respond when possible.
