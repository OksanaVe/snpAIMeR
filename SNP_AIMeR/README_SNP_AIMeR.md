# SNP_AIMeR
This R function assesses the diagnostic power of SNP combinations using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within [adegenet](https://github.com/thibautjombart/adegenet) as well as [pegas](https://github.com/emmanuelparadis/pegas). 

Its value is in (1) identifying ancestry informative markers (AIMs) and (2) evaluating how well different marker combinations can predict an unknown sample's population of origin. 

The user provides candidate markers, SNP genotypes from indivivduals of known origin, and a range of panel sizes. 

SNP_AIMeR tests every marker combination within the specified minimum and maximum panel sizes. To assess the performance of a SNP combination, 20% of the individuals in the dataset are withheld as test data. The predicted population for each test individual is then compared with the known population label from the input file.

Because of the number of possible combinations, we recommend testing no more than 15 markers. For example, testing 15 markers in panel sizes of 1 to 15 (32,767 total combinations) with 1,000 cross-validation replicates on a system with 48 processor cores took about 5 hours and 20 GB RAM. To mitigate run time, SNP_AIMeR automatically uses n - 1 the number of available processor cores. 


## Requirements
.stru (STRUCTURE) formatted genotype file. Individuals must have population assignments.


## Run interactively (user-friendly)
```
>SNP_AIMeR("interactive")
```

Upon executing the function, the user is prompted with the following:
```
Enter path to working directory: 
Enter path to STRUCTURE file: 
```
Then, the user is prompted (by adegenet) for information about the SNP genotype file:
```
How many genotypes are there? 
How many markers are there? 
Which column contains labels for genotypes ('0' if absent)? 
Which column contains the population factor ('0' if absent)? 
Which other optional columns should be read (press 'return' when done)? 
Which row contains the marker names ('0' if absent)? 
Are genotypes coded by a single row (y/n)? 
```
Finally, after a few messages about the data (again from adegenet), the user is prompted for the following:
```
Minimum number of markers in combination:
Maximum number of markers in combination: 
Enter assignment rate threshold (minimum rate of successful assignments): 
```

## Run without interaction
```
>SNP_AIMeR("non-interactive", "config_file")
```

Non-interactive mode requires a config file in YAML format. Example here
```
min_range: 1                                                  # Minimum combination size
max_range: 2                                                  # Maximum combination size
assignment_rate_threshold: 0.9                                # Value from 0 to 1
cross_validation_replicates: 1000                             # We recommend no less than 100 replicates

structure_file: "SNP_check_toy_dataset_176inds_15SNPs.str"    # File path is in quotes
number_of_individuals: 176                                    # Same as adegenet's "n.ind"
number_of_loci: 15                                            # Same as adegenet's "n.loc"
one_data_row_per_individual: FALSE                            # TRUE or FALSE
column_sample_IDs: 1                                          # Column number with individual sample names
column_population_assignments: 2                              # Column number with individual population of origin
column_other_info:                                            # Column number
row_markernames: 1                                            # Row number with marker names
no_genotype_character: -9                                     # Default is "-9"
optional_population_info:                                     # Optional
genotype_character_separator:                                 # Optional
```

## Output
During the analysis, after each panel size is evaluated, replicate cross-validation data for the last combination tested is displayed in a histogram. Cross-validation data is not saved, however, mean values are available in the output files "All_combinations_assignment_rate.csv" and "Combination_assignment_rate_means.csv". The output file "Above_threshold_assignment_rate.csv" lists the combinations with a mean correct assignment rate above the user-specifed threshold.

"Single_marker_assignment_rate.pdf" is each candidate marker's individual assignment rate. This is the same as running min_range=1, max_range=1.
<br clear="left"/>
<img src="https://github.com/OksanaVe/SNP_check/assets/131922755/0526f289-f2c7-45ff-95f0-0214b5d4a328" align="left" width="25%" height="25%" />
<br clear="left"/>

<br clear="left"/>
"Combination_assignment_rate_means.pdf" summarizes the average assignment rate for each tested panel size.
<br clear="left"/>
<img align="left" src="https://github.com/OksanaVe/SNP_check/assets/131922755/8d89dade-ca46-42e0-a9fd-703e0f2a38e7" width="25%" height="25%" />
<br clear="left"/>


## Toy dataset
A toy dataset of 15 SNPs and 176 individuals is provided [here](https://github.com/OksanaVe/SNP_check/blob/main/SNP_check_toy_dataset_176inds_15SNPs.str). The example YAML file is already setup for this dataset. For interactive mode, use
the following prompt responses.
```
library(SNP_AIMeR)
SNP_AIMeR("interactive")
Enter path to working directory: ./
Enter path to STRUCTURE file: SNP_check_toy_dataset_176inds_15SNPs.str

How many genotypes are there? 
176

How many markers are there? 
15

Which column contains labels for genotypes ('0' if absent)? 
1

Which column contains the population factor ('0' if absent)? 
2

Which other optional columns should be read (press 'return' when done)? 
1: 

Which row contains the marker names ('0' if absent)? 
1

Are genotypes coded by a single row (y/n)? 
n

Converting data from a STRUCTURE .stru file to a genind object... 

Data file contains  15  markers
File contains the following group definitions:
group_1 group_2 
    120      56 
Minimum number of markers in combination: 10
Maximum number of markers in combination: 15
Enter assignment rate threshold (minimum rate of successful assignments): 0.9


