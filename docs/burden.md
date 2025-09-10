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
