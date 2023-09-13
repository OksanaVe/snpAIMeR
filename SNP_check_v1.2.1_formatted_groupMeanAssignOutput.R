SNP_check <- function() {

	require(adegenet)
	require(pegas)
	require(foreach)
	require(doParallel)
	require(parallel)

	#Setup backend to use many processors
	print(parallel::detectCores())
	n.cores <- parallel::detectCores() - 1
	#create the cluster
	my.cluster <- parallel::makeCluster(
		n.cores, 
		type = "PSOCK"
	)
	print(my.cluster)

	suppressWarnings(doParallel::registerDoParallel(cl = my.cluster))

	working_dir <- readline(prompt="Enter path to working directory: ")
	setwd(working_dir)
	
	str_fn <- readline(prompt="Enter path to STRUCTURE file: ")
	ludens <- read.structure(str_fn)
	loci=locNames(ludens)
	cat("Data file contains ", length(loci), " markers\n") 
	cat("File contains the following group definitions:")
	print(table(ludens@pop))

	min_loc_num <- readline(prompt="Minimum number of markers in combination: ")
	min_loc_num <- as.integer(min_loc_num)
	max_loc_num <- readline(prompt="Maximum number of markers in combination: ")
	max_loc_num <- as.integer(max_loc_num)
	rate_threshold <- readline(prompt="Enter assignment rate threshold (minimum rate of successful assignments): ")
	rate_threshold <- as.double(rate_threshold)

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

	}
	df <- cbind(mn_rate$markers, mn_rate$mean_rate)
	colnames(df) <- c("group","mean")
	write.csv(df, "Combination_assignment_rate_means.csv", row.names=FALSE)

	plot(mn_rate$markers, mn_rate$mean_rate, col="red", xlab="number of markers in combination", ylab="avg assignment rate", pch=16, type="b",lwd=1.5,lty=3, main="Avg assignment rate vs number of markers")
	dev.copy(pdf, paste("Combinations.pdf"))
	dev.off()

	write.csv(all_loci, "All_results_marker_assignment_rate.csv")
	write.csv(good_loci, "Above_threshold_marker_assignment_rate.csv")
	cat(nrow(good_loci), " marker combinations passed threshold")

	suppressWarnings(parallel::stopCluster(cl = my.cluster))
}

