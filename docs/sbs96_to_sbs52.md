## Mapping SBS96 classification to SBS52 classification
Please refer to the methods section of the manuscript or the `get_sbs96_to_sbs52_lookup_table.py` script for a detailed description of how the SBS96 classification is mapped to the SBS52 classification.

### Collapse SBS96 classification counts to SBS52 classification counts

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
