# SNP_AIMeR
This R script assesses the diagnostic power of SNP combinations using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within [adegenet](https://github.com/thibautjombart/adegenet) as well as [pegas](https://github.com/emmanuelparadis/pegas). 
Its utility is in identifying ancestry informative markers (AIMs) to include in marker panels designed to identify the source population of unknown samples.

The user provides a range of panel sizes to test, and within that range, SNP_AIMeR tests every marker combination for how well it can assign the correct population to samples with known origin.

Due to computational size, we recommend testing no more than 15 markers. (add runtime example (48 cores, 192 GB RAM))


**Requires**
1. .str (STRUCTURE) formatted genotype file. Individuals must have population assignments.
2. Config file 
```
min_range: 1                                                  # Minimum panel size
max_range: 15                                                 # Maximum panel size
assignment_rate_threshold: 0.9                                # Value from 0 to 1
cross_validation_replicates: 100                              # We recommend no less than 100 replicates

structure_file: "SNP_check_toy_dataset_176inds_15SNPs.str"
number_of_individuals: 176                                    # Same as adegenet's "n.ind"
number_of_loci: 15                                            # Same as adegenet's "n.loc"
one_data_row_per_individual: FALSE                            # TRUE or FALSE
column_sample_IDs: 1                                          # Column number with individual sample names
column_population_assignments: 2                              # Column number with population information
column_other_info:                                            # Column number
row_markernames: 1                                            # Row number with marker names
no_genotype_character: -9                                     # Default is "-9"
optional_population_info:                                     # Optional
genotype_character_separator:                                 # Optional
```

**To Run**
library(SNP_AIMeR)
SNP_AIMeR(config_file)

**Output**



