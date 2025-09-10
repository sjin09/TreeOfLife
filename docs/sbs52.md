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

### SBS52 counts, where each trinucleotide contributes equally.

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
