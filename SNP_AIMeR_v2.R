SNP_AIMeR <- function(config_file) {
  require(adegenet)
  require(pegas)
  require(foreach)
  require(doParallel)
  require(parallel)
  library(yaml)
  
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
  # adegenet::read.structure settings
  file <- config$structure_file
  number_of_individuals <- as.integer(config$number_of_individuals)
  number_of_loci <- as.integer(config$number_of_loci)
  column_sample_IDs <- as.integer(config$column_sample_IDs)
  column_population_assignments <- as.integer(config$column_population_assignments)
  row_markernames <- as.integer(config$row_markernames)
  column_other_info <- as.integer(config$column_other_info)
  
	ludens <- read.structure(file, n.ind=number_of_individuals, n.loc=number_of_loci, onerowperind=config$one_data_row_per_individual, col.lab=column_sample_IDs, col.pop=column_population_assignments, col.others=column_other_info, NA.char=config$no_genotype_character, pop=config$optional_population_info, sep=config$genotype_character_separator, ask=FALSE, quiet=FALSE)

	loci=locNames(ludens)
	cat("Data file contains ", length(loci), " markers\n") 
	cat("File contains the following group definitions:")
	print(table(ludens@pop))

	all_loci = data.frame(markers=character(), avg_success_rate=double())
	good_loci = data.frame(markers=character(), avg_success_rate=double())
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

		for(i in 1:ncol(loc_comb)){
			tryCatch({
				test = ludens[,loc=loc_comb[,i]]

				results <- foreach(j=1:100, .combine='c') %dopar% {
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
				if(mean_rate >= rate_threshold) {
					good_loci[nrow(good_loci)+1,] = c(markers, mean_rate)
				}
			}, error=function(e){})
		}
		title <- paste0("Assignment rate for ", n, " marker(s)")
		if(n == 1 && length(rate_ls) < 16) {
			boxplot(rate_ls, col="steelblue", main=title) 
			mtext(paste("Mean =", round(mean(results), 4)), side=3, col="red")
			dev.copy(pdf, paste0(title, ".pdf"))
			dev.off()
		} 
		else if(n == 1 && length(rate_ls) >= 16) {
			y = split(rate_ls, rep(1:ceiling(length(rate_ls)/15), each=15)[1:length(rate_ls)])
			for(w in 1:length(y)){
				boxplot(y[[w]], col="steelblue", main=title)
				boxplot(y[[w]], col="steelblue", main=title)
				dev.copy(pdf, paste0(title, ".pdf"))
				dev.off()
			}
		}
		else if(n > 1) {hist(results, main=title) 
			abline(v = mean(results), col = "red", lwd = 2)
			mtext(paste("Mean =", round(mean(results), 4)), side=3, col="red")
		}
		mn_rate[nrow(mn_rate)+1,] = c(n,mean(results))
		
		df <- cbind(mn_rate$markers, mn_rate$mean_rate)
		colnames(df) <- c("group","mean")
		write.csv(df, "Combination_assignment_rate_means.csv", row.names=FALSE)
	
		write.csv(all_loci, "All_results_marker_assignment_rate.csv")
		write.csv(good_loci, "Above_threshold_marker_assignment_rate.csv")

	}
  plot(mn_rate$markers, mn_rate$mean_rate, col="red", xlab="number of markers in combination", ylab="avg assignment rate", pch=16, type="b",lwd=1.5,lty=3, main="Avg assignment rate vs number of markers")
	dev.copy(pdf, paste("Combinations.pdf"))
	dev.off()

	cat(nrow(good_loci), " marker combinations passed threshold")
	
	suppressWarnings(parallel::stopCluster(cl = my.cluster))
}
