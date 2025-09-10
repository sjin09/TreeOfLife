## Mutational signature extraction
Given a matrix with germline and somatic mutations categorized under SBS52 and SBS96 classifications, respectively, HDP can *de novo* extract mutational signatures.

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
