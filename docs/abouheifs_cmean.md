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
