SNP_check <- function() {

require(adegenet)
require(pegas)

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

for(n in min_loc_num:max_loc_num) {

loc_comb = combn(loci,n)

for(i in 1:ncol(loc_comb)){
tryCatch({
rate <- c()
test = ludens[,loc=loc_comb[,i]]
for(j in 1:1000) {
	tryCatch({
    pop_counts <- table(ludens@pop)
    counts = as.vector(pop_counts)
    pop_sample = sort(round(0.75*counts, 0))
    #kept.id <- runif(100, min=1, max=120)
    kept.id <- unlist(tapply(1:nInd(test), pop(test), function(e) sample(e, pop_sample,replace=FALSE)))
    x <- test[kept.id]
    x.sup <- test[-kept.id]
    dapc4 <- dapc(x,n.pca=50,n.da=20)
    pred.sup <- predict.dapc(dapc4, newdata=x.sup)
    a <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
    rate <- c(rate, a)}, error=function(e){})
}
mean_rate <- sum(rate)/length(rate)
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
hist(rate, main=title)
abline(v = mean(rate), col = "red", lwd = 2)
mtext(paste("Mean =", round(mean(rate), 4)), side=3, col="red")
mn_rate[nrow(mn_rate)+1,] = c(n,mean(rate))
}
plot(mn_rate$markers, mn_rate$mean_rate, col="red", xlab="number of markers in combination", ylab="avg assignment rate", pch=16, type="b",lwd=1.5,lty=3, main="Avg assignment rate vs number of markers")
write.csv(all_loci, "All_results_marker_assignment_rate.csv")
write.csv(good_loci, "Above_threshold_marker_assignment_rate.csv")
cat(nrow(good_loci), " marker combinations passed threshold")
}
