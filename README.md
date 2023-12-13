<img align="left" src="https://user-images.githubusercontent.com/131922755/278390582-0652ccaf-a4f4-41ae-8ba9-f8d81e7a115f.png" width="40%" height="40%" />
<br clear="left"/>
<br clear="left"/>

This R package assesses the diagnostic power of SNP combinations using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within the [adegenet](https://github.com/thibautjombart/adegenet) R package.
Its value is in (1) identifying ancestry informative markers (AIMs) and (2) evaluating how well different marker combinations can predict an unknown sample's population of origin. 

The user provides candidate markers, SNP genotypes from individuals of known origin, a range of panel sizes, and a threshold value for an acceptable rate of correct sample identification.

snpAIMeR tests every marker combination within the specified minimum and maximum panel sizes. For each cross-validation replicate, individuals are randomly divided with 80% for the DAPC and 20% withheld as test samples. Results from the DAPC are used to predict the population of origin for each test individual, which is then compared with the known population label from the input file.

Because of the number of possible combinations, we recommend testing no more than 15 markers. For example, testing 15 markers in panel sizes of 1 to 15 (32,767 total combinations) with 1,000 cross-validation replicates on a system with 48 processor cores took about 5 hours and 20 GB RAM. To mitigate run time, snpAIMeR automatically uses n - 1 the number of available processor cores. Reducing the number of cross-validation replicates also reduces run time, however, we recommend no less than 100 replicates.


## Requirements
.stru (STRUCTURE) formatted genotype file. Individuals must have population assignments.


## Usage
```
snpAIMeR(run_mode, config_file = NULL, verbose = TRUE)
```

## Run interactively (user-friendly)
```
>snpAIMeR("interactive")
```
Upon executing the function, the user is prompted with the following (do not quote paths):
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
Finally, after a few messages about the data (again from adegenet), the user is prompted for the following (we recommend no less than 100 cross-validation replicates):
```
Minimum number of markers in combination:
Maximum number of markers in combination:
Assignment rate threshold (minimum rate of successful assignments):
Number of cross-validation replicates:
```

## Run without interaction
```
>snpAIMeR("non-interactive", "config_file")
```
Non-interactive mode requires a config file in YAML format. Example [here](https://github.com/OksanaVe/snpAIMeR/blob/main/snpAIMeR_config.yml)
```
min_range: 1                                                  # Minimum combination size
max_range: 5                                                  # Maximum combination size
assignment_rate_threshold: 0.9                                # Value from 0 to 1
cross_validation_replicates: 100                              # We recommend no less than 100 replicates
working_directory: "./"                                       # Path name in quotes; use "./" for current directory

structure_file: "snpAIMeR_toy_dataset_176inds_5SNPs.str"      # Path name in quotes
number_of_individuals: 176                                    # Same as adegenet's "n.ind"
number_of_loci: 5                                             # Same as adegenet's "n.loc"
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
For each panel size, when all the combinations have been evaluated, replicate cross-validation data for the last combination is displayed as a histogram.
* "All_combinations_assign_rate.csv" has the mean correct assignment rate for each combination tested (average of all cross-validation replicates)
* "Panel_size_assign_rate.csv" has the mean correct assignment rate for each panel size tested (average of all combinations)
* "Above_threshold_assign_rate.csv" lists the combinations with a mean correct assignment rate above the user-specified threshold.

"Single_marker_assignment_rate.pdf" is each candidate marker's individual assignment rate. This is the same as running min_range=1, max_range=1.
<br clear="left"/>
<img src="https://github.com/OksanaVe/SNP_check/assets/131922755/0526f289-f2c7-45ff-95f0-0214b5d4a328" align="left" width="25%" height="25%" />
<br clear="left"/>

"Panel_size_assign_rate.pdf" is a visualization of "Panel_size_assign_rate.csv"
<br clear="left"/>
<img align="left" src="https://github.com/OksanaVe/SNP_check/assets/131922755/8d89dade-ca46-42e0-a9fd-703e0f2a38e7" width="25%" height="25%" />
<br clear="left"/>


## Toy dataset
A toy dataset of 5 SNPs and 176 individuals is provided [here](https://github.com/OksanaVe/snpAIMeR/blob/main/snpAIMeR_toy_dataset_176inds_5SNPs.str). The example YAML file is already setup for this dataset. For interactive mode, use the following prompt responses.
```
snpAIMeR("interactive")
Enter path to working directory: ./
Enter path to STRUCTURE file: snpAIMeR_toy_dataset_176inds_15SNPs.str

How many genotypes are there? 
176

How many markers are there? 
5

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

Data file contains  5  markers
File contains the following group definitions:
group_1 group_2 
    120      56 
Minimum number of markers in combination: 1
Maximum number of markers in combination: 5
Assignment rate threshold (minimum rate of successful assignments): 0.9
Number of cross-validation replicates: 100

