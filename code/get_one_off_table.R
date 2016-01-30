getRegionNonRandom <- function(region){
	file_name <- paste0("data/mothur_pool/", region, ".mock.screen.unique.good.filter.unique.error.summary")

	error <- read.table(file=file_name, header=T, row.names=1)
	error.nochim <- error[error$numparents==1,]

	one.off <- error.nochim[error.nochim$mismatches==1,]

	one.off.n <- sum(one.off$weight)
	one.off.table <- table(one.off$weight)
	nseqs <- sum(error.nochim$weight)
	one.off.table
}

regions <- c("V1V3", "V1V5", "V1V6", "V1V9", "V3V5", "V4")
non_random <- lapply(regions, getRegionNonRandom)
names(non_random) <- regions

max_freq <- max(sapply(non_random, function(x){max(as.numeric(names(x)))}))
one_offs <- matrix(0, nrow=max_freq, ncol=length(regions))
rownames(one_offs) <- 1:max_freq
colnames(one_offs) <- regions

for(r in regions){
	data <- unlist(unname(non_random[r]))
	one_offs[as.numeric(names(data)),r] <- data
	one_offs[,r]
}

write.table(one_offs, "data/process/non_random_analysis.tsv", sep="\t", quote=F)

