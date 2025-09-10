### Reference FASTA file trinucleotide (3-mer) counts

There are 64 possible trinucleotides (4 ** 3). The script below counts the occurrences of all 64 trinucleotides. However, if the middle base is a purine (adenine or guanine), the count is incremented using the reverse complement of the trinucleotide. Hence, the counts of trinucleotides where the middle base is a pyrimidine are returned.

```
python get_reference_tricounts.py \
    -i sample.fasta \
    --target sample.target \
    -o sample.tri
```
