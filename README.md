# SNP_AIMeR
This R script assesses the diagnostic power of SNP combinations using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within [adegenet](https://github.com/thibautjombart/adegenet) as well as [pegas](https://github.com/emmanuelparadis/pegas). 

Its value is in (1) identifying ancestry informative markers (AIMs) and (2) evaluating marker combinations for how well they can predict an unknown sample's population of origin. 

The user provides candidate markers, genotypes from known populations, and a range of marker combination sizes to test. Within that size range, SNP_AIMeR tests every marker combination.

Due to the number of possible combinations, we recommend testing no more than 15 markers. For example, testing 15 markers in combination sizes from 1 to 15 (32,767 total combinations) and 1000 cross-validation replicates on a system with 48 processor cores took about 5 hours and used about 20 GB of memory. To mitigate run time, SNP_AIMeR automatically uses n - 1 the number of available processor cores. 


## Requirements
.stru (STRUCTURE) formatted genotype file. Individuals must have population assignments.


## Run interactively
SNP_AIMeR("interactive")

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
Finally, after a few messages (again from adegenet) about the data, the user will be prompted for the following:
```
Minimum number of markers in combination:
Maximum number of markers in combination: 
Enter assignment rate threshold (minimum rate of successful assignments): 
```

## Run without interaction
SNP_AIMeR("non-interactive", config_file)

To run SNP_AIMeR non-interactively, a config file in YAML format is required. Example <here>
```
min_range: 1                                                  # Minimum combination size
max_range: 15                                                 # Maximum combination size
assignment_rate_threshold: 0.9                                # Value from 0 to 1
cross_validation_replicates: 1000                             # We recommend no less than 100 replicates

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

With those inputs, it will randomize combinations of SNPs between the specified minimum and maximum numbers and use a priori population delimitations (via the .str file) to assess whether those SNP combinations lead to successful reassignment of a test portion of the data (20% of individuals in the dataset). 

**Output**
During the analysis, each cluster size's replicate test data assignment rate data is displayed after the analysis for that size is complete. Plots are saved but the date is in the output csvs.


Each candidate marker's individual assignment rate. This is the same as running min_range=1, max_range=1. Each data point is the test data assignment rate for a single replicate analysis.

Summary of average assignment rate for each cluster size tested


**Toy dataset**
A toy dataset of 15 SNPs and 176 individuals is provided here. The following commands and prompt responses will facilitate analyzing this dataset.
```
> library(SNP_AIMeR)
> SNP_AIMeR("interactive")
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







