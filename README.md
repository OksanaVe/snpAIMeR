# SNP_check
This R script assesses the diagnostic power of SNP combinations from a .str formatted genotype file using leave-one-out style cross-validation. To do so, it uses [Discriminant Analysis of Principal Components](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94) within [adegenet](https://github.com/thibautjombart/adegenet) as well as [pegas](https://github.com/emmanuelparadis/pegas). 

Upon running the script, the user will be prompted with the following:

```
Enter path to working directory: 
Enter path to STRUCTURE file: 
Minimum number of markers in combination:
Maximum number of markers in combination: 
Enter assignment rate threshold (minimum rate of successful assignments): 
```

With those inputs, it will randomize combinations of SNPs between the specified minimum and maximum numbers and use a priori population delimitations (via the .str file) to assess whether those SNP combinations lead to successful reassignment of a test portion of the data (generally 25 % of individuals in the dataset). 

![assignment_success](Rate_vs_NumOfMarkers_SET15_wild_rearing_356inds.png)

Output files include csv files with marker assignment rates and summarized results. 
