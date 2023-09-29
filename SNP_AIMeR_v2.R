SNP_AIMeR <- function(config_file) {
  library(yaml)
  require(doParallel)
  require(parallel)
  require(foreach)
  require(adegenet)
  require(pegas)
  require(tidyverse)
  
  # Setup backend to use many processors
  print(parallel::detectCores())
  n.cores <- parallel::detectCores() - 1
  # Create the cluster
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  print(my.cluster)
  suppressWarnings(doParallel::registerDoParallel(cl = my.cluster))
  
  config = yaml.load_file(config_file)
  
  # SNP_AIMeR settings
  min_loc_num <- as.integer(config$min_range)
  max_loc_num <- as.integer(config$max_range)
  assignment_rate_threshold <- as.double(config$assignment_rate_threshold)
  cv_replicates <- as.integer(config$cross_validation_replicates)
  # adegenet::read.structure settings
  file <- config$structure_file
  number_of_individuals <- as.integer(config$number_of_individuals)
  number_of_loci <- as.integer(config$number_of_loci)
  column_sample_IDs <- as.integer(config$column_sample_IDs)
  column_population_assignments <- as.integer(config$column_population_assignments)
  row_markernames <- as.integer(config$row_markernames)
  column_other_info <- as.integer(config$column_other_info)
  
	ludens <- read.structure(file, n.ind=number_of_individuals, n.loc=number_of_loci, onerowperind=config$one_data_row_per_individual, col.lab=column_sample_IDs, col.pop=column_population_assignments, col.others=column_other_info, row.marknames=row_markernames, NA.char=config$no_genotype_character, pop=config$optional_population_info, sep=config$genotype_character_separator, ask=FALSE, quiet=FALSE)

	loci=locNames(ludens)
	cat("Data file contains", length(loci), "markers\n") 
	cat("File contains the following group definitions:")
	print(table(ludens@pop))

	df = data.frame(marker=character(), avg_success_rate=double())
  write.table(df,"Group_assignment_rate_means.csv", row.names=FALSE, sep = ",")
	write.table(df,"All_results_marker_assignment_rate.csv", row.names=FALSE, sep = ",")
	write.table(df,"Above_threshold_marker_assignment_rate.csv", row.names=FALSE, sep = ",")

	mn_rate = data.frame(markers=integer(), mean_rate=double())
	
	if(min_loc_num == 1) {
		range = c(min_loc_num:max_loc_num)
	} 
	else if(min_loc_num > 1) {
		range = c(1, min_loc_num:max_loc_num)
	} 
	else if(min_loc_num < 1) {
		print("Incorrect number of markers provided")
		}
	
	for(n in range) {
		rate_ls <- list()
		loc_comb = combn(loci,n)
		all_loci = data.frame(markers=character(), avg_success_rate=double())
		good_loci = data.frame(markers=character(), avg_success_rate=double())
		
		for(i in 1:ncol(loc_comb)){
			tryCatch({
				test = ludens[,loc=loc_comb[,i]]

				results <- foreach(j=1:cv_replicates, .combine='c') %dopar% {
					tryCatch({
						kept.id <- unlist(tapply(1:adegenet::nInd(test), adegenet::pop(test), function(e) sample(e, round(0.75*length(e),0), replace=FALSE)))
						x <- test[kept.id]
						x.sup <- test[-kept.id]
						dapc4 <- adegenet::dapc(x,n.pca=50,n.da=20)
						pred.sup <- adegenet::predict.dapc(dapc4, newdata=x.sup)
						a <- mean(as.character(pred.sup$assign)==as.character(adegenet::pop(x.sup)))
					}, error=function(e){})
					return(a)
				}
				rt = list(results)
				rate_ls <- append(rate_ls, rt)
				mean_rate <- sum(results)/length(results)
				print(paste0("Combination ", i))
				markers = toString(loc_comb[,i])
				cat(markers, "\n")
				cat(mean_rate, "\n")
				all_loci[nrow(all_loci)+1,] = c(markers, mean_rate)
				if(mean_rate >= assignment_rate_threshold) {
					good_loci[nrow(good_loci)+1,] = c(markers, mean_rate)
				}
			}, error=function(e){})
		}
		# Make box plot of each marker's cross-validation replicates
		if(n == 1) {
		  #Make the assignment rate tidy
		  loc_comb_list <- as.list(loc_comb)
		  loc_comb_df <- as.data.frame(do.call(rbind, loc_comb_list))
		  #print("LOC_COMB_DF")
		  #print(loc_comb_df)
		  names(loc_comb_df)[1] <- "markerID"
		  rate_df <- as.data.frame(do.call(rbind, rate_ls))
		  #print("RATE_DF")
		  #print(rate_df)
		  rate_df <- cbind(loc_comb_df, rate_df)
		  rate_df <- pivot_longer(rate_df, -markerID, names_prefix="V", names_to="replicate", values_to="rate")
		  # Group the data by marker
		  rate_df_group_by_markerID <- rate_df %>% group_by(markerID)
		  # Get mean of all markers/replicates
		  mean_results <- round(mean(results), 4)
		  # Make box plot
		  title <- paste0("Mean = ", mean_results)
		  each_marker_plot <- ggplot(rate_df_group_by_markerID, aes(x=fct_reorder(markerID, parse_number(markerID)), y=rate)) +
		  geom_boxplot(fill="gray", outlier.size=0.5) +
	    labs(title=title, x="Marker", y="Assignment rate") +
		  theme_linedraw() +
		  theme(axis.text.x = element_text(angle = 45, hjust = 1))
		  
	    ggsave("Assignment_rate_for_each_marker.pdf")
	    print(each_marker_plot)
	   }
    # Histogram of rate distribution. Only prints to screen (intentional)
		else if(n > 1) {hist(results, main=title) 
			abline(v = mean(results), col = "red", lwd = 2)
			mtext(paste("Mean =", round(mean(results), 4)), side=3, col="red")
		}
		# Get the mean assignment rate for each size group
		mn_rate[nrow(mn_rate)+1,] = c(n, mean(results))
		df1 <- cbind(mn_rate$markers, mn_rate$mean_rate)
		colnames(df1) <- c("combination","mean")
		df1 <- as.data.frame(df1)
		
		write.table(df1, "Group_assignment_rate_means.csv", row.names=FALSE, col.names=FALSE, append=TRUE, sep = ",")
		write.table(all_loci, "All_results_marker_assignment_rate.csv", row.names=FALSE, col.names=FALSE, append=TRUE, sep = ",")
		write.table(good_loci, "Above_threshold_marker_assignment_rate.csv", row.names=FALSE, col.names=FALSE, append=TRUE, sep = ",")
	}
  # Make line plot of combination assignment rate means
	if (nrow(df1) > 1) {
	  combination_means <- ggplot(df1, aes(x=combination, y=mean)) +
      geom_line() +
	    geom_point() +
	    labs(title="Average assignment rate vs number of markers", x="Number of markers in combination", y="Average assignment rate") +
	    theme(panel.border=element_rect(color="black", fill=NA, linewidth=1))
	  ggsave("Combination_assignment_rate_means.pdf")
	  print(combination_means)
	}

	cat(nrow(good_loci), "marker combinations passed threshold","\n")

	suppressWarnings(parallel::stopCluster(cl = my.cluster))
}
