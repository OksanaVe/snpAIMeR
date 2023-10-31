snpAIMeR <- function(run_mode, config_file=NULL) {
  library(yaml)
  require(doParallel)
  require(parallel)
  require(foreach)
  require(adegenet)
  require(pegas)
  #require(tidyverse)
  require(dplyr)
  require(tidyr)
  require(readr)
  require(forcats)
  require(ggplot2)
  
  # Setup backend to use many processors
  print(parallel::detectCores())
  n.cores <- parallel::detectCores() - 1
  # Create the cluster
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  print(my.cluster)
  suppressWarnings(doParallel::registerDoParallel(cl = my.cluster))
  
  run_mode <- tolower(gsub('[ -]', '', run_mode))
  
  if((run_mode == "noninteractive") & (is.null(config_file))) {
    print("Error: non-interactive mode requires config file")
  }
  if ((run_mode == "noninteractive") & (!is.null(config_file))) {
    config = yaml.load_file(config_file)

    # SNP_AIMeR settings
    setwd(config$working_directory)
    min_loc_num <- as.integer(config$min_range)
    max_loc_num <- as.integer(config$max_range)
    assignment_rate_threshold <- as.double(config$assignment_rate_threshold)
    cv_replicates <- as.integer(config$cross_validation_replicates)
    # adegenet::read.structure settings
    pop_file <- config$structure_file
    number_of_individuals <- as.integer(config$number_of_individuals)
    number_of_loci <- as.integer(config$number_of_loci)
    column_sample_IDs <- as.integer(config$column_sample_IDs)
    column_population_assignments <- as.integer(config$column_population_assignments)
    row_markernames <- as.integer(config$row_markernames)
    column_other_info <- as.integer(config$column_other_info)
    pop_data <- read.structure(pop_file, n.ind=number_of_individuals, n.loc=number_of_loci, onerowperind=config$one_data_row_per_individual, col.lab=column_sample_IDs, col.pop=column_population_assignments, col.others=column_other_info, row.marknames=row_markernames, NA.char=config$no_genotype_character, pop=config$optional_population_info, sep=config$genotype_character_separator, ask=FALSE, quiet=FALSE)

    loci=locNames(pop_data)
    cat("Data file contains", length(loci), "markers\n") 
    cat("File contains the following group definitions:")
    print(table(pop_data@pop))
  }
  
  if(run_mode == "interactive") {
    working_dir <- readline(prompt="Enter path to working directory: ")
    setwd(working_dir)
    str_fn <- readline(prompt="Enter path to STRUCTURE file: ")
    pop_data <- read.structure(str_fn)
    loci=locNames(pop_data)
    cat("Data file contains ", length(loci), " markers\n") 
    cat("File contains the following group definitions:")
    print(table(pop_data@pop))
  
    min_loc_num <- as.integer(readline(prompt="Minimum number of markers in combination: "))
    max_loc_num <- as.integer(readline(prompt="Maximum number of markers in combination: "))
    assignment_rate_threshold <- as.double(readline(prompt="Threshold value for the minimum rate of successful assignments: "))
	cv_replicates <- as.integer(readline(prompt="Number of cross-validation replicates: "))
  }

	df = data.frame(marker=character(), avg_success_rate=double())
  write.table(df,"Combination_assignment_rate_means.csv", row.names=FALSE, sep = ",")
	write.table(df,"All_combinations_assignment_rate.csv", row.names=FALSE, sep = ",")
	write.table(df,"Above_threshold_assignment_rate.csv", row.names=FALSE, sep = ",")

	mn_rate = data.frame(markers=integer(), mean_rate=double())
	
	if(min_loc_num == 1) {
		range = c(min_loc_num:max_loc_num)
	} else if(min_loc_num > 1) {
		range = c(1, min_loc_num:max_loc_num)
	} else if(min_loc_num < 1) {
		print("Incorrect number of markers provided")
	}

	for(n in range) {
		rate_ls <- list()
		loc_comb = combn(loci,n)
		all_loci = data.frame(markers=character(), avg_success_rate=double())
		good_loci = data.frame(markers=character(), avg_success_rate=double())
		
		for(i in 1:ncol(loc_comb)){
			test = pop_data[,loc=loc_comb[,i]]

			results <- foreach(j=1:cv_replicates, .combine='c') %dopar% {
				tryCatch({
					kept.id <- unlist(tapply(1:adegenet::nInd(test), adegenet::pop(test), function(e) sample(e, round(0.75*length(e),0), replace=FALSE)))
					x <- test[kept.id]
					x.sup <- test[-kept.id]
					dapc4 <- adegenet::dapc(x,n.pca=50,n.da=20)
					pred.sup <- adegenet::predict.dapc(dapc4, newdata=x.sup)
					a <- mean(as.character(pred.sup$assign)==as.character(adegenet::pop(x.sup)))
				}, error=function(e){
				return(0)
			})
		}
			rt = list(results)
			rate_ls <- append(rate_ls, rt)
			
			markers = toString(loc_comb[,i])
			mean_rate <- sum(results)/length(results)
			print(paste0("Combination ", i))
			cat(markers, "\n")
			cat(mean_rate, "\n")
		  # Write data to output csv files
			all_loci[nrow(all_loci)+1,] = c(markers, mean_rate)
			if(mean_rate >= assignment_rate_threshold) {
				good_loci[nrow(good_loci)+1,] = c(markers, mean_rate)
			}
		}
		# Make box plot of each marker's cross-validation replicates
		if(n == 1) {
		  #Make the assignment rate tidy
		  loc_comb_list <- as.list(loc_comb)
		  loc_comb_df <- as.data.frame(do.call(rbind, loc_comb_list))
		  names(loc_comb_df)[1] <- "marker_name"
		  rate_df <- as.data.frame(do.call(rbind, rate_ls))
		  rate_df <- cbind(loc_comb_df, rate_df)
		  rate_df <- pivot_longer(rate_df, -marker_name, names_prefix="V", names_to="replicate", values_to="rate")
		  # Group the data by marker
		  rate_df_group_marker <- rate_df %>% group_by(marker_name)
		  # Get mean of all markers/replicates
		  mean_results <- round(mean(results), 4)
		  # Make box plot
		  box_plot_title <- paste0("Mean for all markers = ", mean_results)
		  each_marker_plot <- ggplot(rate_df_group_marker, aes(x=fct_reorder(marker_name, parse_number(marker_name)), y=rate)) +
		    geom_boxplot(fill="gray", outlier.size=0.5) +
	      labs(title=box_plot_title, x="Marker", y="Assignment rate") +
		    theme_light() +
		    theme(axis.text.x = element_text(angle = 45, hjust = 1))
	    ggsave("Single_marker_assignment_rate.pdf")
	    print(each_marker_plot)
	   }
    # Histogram of rate distribution. Only prints to screen.
		else if(n > 1) {hist(results, main=paste0("Combinations of ", n, " markers, " , cv_replicates, " cross-validation replicates"), xlab="Assignment rate", ylab="Frequency")
			abline(v = mean(results), col = "red", lwd = 2)
			mtext(paste("Mean =", round(mean(results), 4)), side=3, col="red")
		}
		# Get the mean assignment rate for each size group
		mn_rate[nrow(mn_rate)+1,] = c(n, mean(results))
		combination_mean_matrix <- cbind(mn_rate$markers, mn_rate$mean_rate)
		colnames(combination_mean_matrix) <- c("combination","mean")
		combination_mean_df <- as.data.frame(combination_mean_matrix)
		
		write.table(combination_mean_df, "Combination_assignment_rate_means.csv", row.names=FALSE, col.names=FALSE, append=FALSE, sep = ",")
		write.table(all_loci, "All_combinations_assignment_rate.csv", row.names=FALSE, col.names=FALSE, append=TRUE, sep = ",")
		write.table(good_loci, "Above_threshold_assignment_rate.csv", row.names=FALSE, col.names=FALSE, append=TRUE, sep = ",")
	}
  # Make line plot of combination assignment rate means
	if (nrow(combination_mean_df) > 1) {
	  combination_means <- ggplot(combination_mean_df, aes(x=combination, y=mean)) +
      geom_line() +
	    geom_point() +
	    scale_x_continuous(breaks=seq(round(max(combination_mean_df$combination),0))) +
	    labs(title="Number of markers vs average assignment rate", x="Number of markers in combination", y="Average assignment rate") +
	    theme_light()
	  ggsave("Combination_assignment_rate_means.pdf")
	  print(combination_means)
	}

	cat(nrow(good_loci), "Marker combinations that passed the assignment rate threshold","\n")

	suppressWarnings(parallel::stopCluster(cl = my.cluster))
}
