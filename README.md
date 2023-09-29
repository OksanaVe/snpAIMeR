# SNP_AIMeR
This R script assesses the diagnostic power of SNP combinations using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within [adegenet](https://github.com/thibautjombart/adegenet) as well as [pegas](https://github.com/emmanuelparadis/pegas). 
Its utility is in (1) identifying ancestry informative markers (AIMs) and (2) evaluating marker combinations for how well they determine the source population of unknown samples.

The user provides candidate markers, genotypes from known populations, and a range of combination sizes to test. Within that size range, SNP_AIMeR tests every marker combination for its ability to assign the correct population to samples with known origin.

Due to the factorial growth rate, we recommend testing no more than 15 markers. For example, testing 15 markers with combination sizes 1 to 15 (15! = 1.308 E+12 combinations total) on a system with 48 cores, 192 GB RAM took xx runtime. To help with this, SNP_AIMeR automatically uses n - 1 the number of available processor cores. 

**Requires**
1. A .stru (STRUCTURE) formatted genotype file. Individuals must have population assignments.
2. SNP_AIMeR config file (YAML format)
```
min_range: 1                                                  # Minimum combination size
max_range: 15                                                 # Maximum combination size
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

**To Run**<br>
library(SNP_AIMeR)<br>
SNP_AIMeR(config_file)

**Output**



