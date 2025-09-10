## Calculate mutation burden per cell
Please note that the calculation here does not account for the number of mutations attributed to technical artefacts, which must be considered to obtain the correct mutation burden per cell for the sample. The number of mutations attributed to technical artefacts can be obtained through either mutational signature extraction or mutational signature attribution.

```
## Himut normcounts adjusted returns sample.observed_sbs96.tsv
## Here, sample.target needs to have both the appropriate autosomes and sex chromosomes to reflect the genome of the sample.

python get_mutation_burden_per_cell.py \
    -i sample.observed_sbs96.tsv \
    --ref-fasta sample.fasta \
    --target sample.target \
    -o sample.burden
```

- Avian sex chromosomes: chrW and chrZ
- Mammalian sex chromosomes: chrX and chrY
