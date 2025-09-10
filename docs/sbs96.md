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

### SBS96 counts, where each trinucleotide contributes equally.

```
python get_tri_equal_weight_sbs96_counts.py \
    -i sample.sbs96.tsv \
    --tri sample.tri \
    -o sample.tri_equal_weight.sbs96.tsv
```
