regions <- c("V1V3", "V1V5", "V1V6", "V1V9", "V3V5", "V4")
pretty_region <- c("V1V3"="V1-V3", "V1V5"="V1-V5", "V1V6"="V1-V6", "V1V9"="V1-V9", "V3V5"="V3-V5", "V4"="V4")

getAveSobs <- function(region, sample){
	file_name <- paste0("data/mothur_pool/", region, ".", sample, ".screen.unique.good.filter.unique.precluster.pick.an.ave-std.summary")

	data_sobs <- NA

	if(file.size(file_name) != 0){
		data <- read.table(file=file_name, header=T)
		data_sobs <- data[data$method=="ave", "sobs"]
	}

	data_sobs
}

getOTUResults <- function(region){
 	perfect <- read.table(file=paste0("data/mothur_pool/HMP_MOCK.", region, ".pick.phylip.an.summary"), header=T)
 	nochims <- read.table(file=paste0("data/mothur_pool/", region, ".mock.screen.unique.good.filter.unique.precluster.perfect.an.ave-std.summary"), header=T)

 	obs.mock <- getAveSobs(region, "mock")
 	obs.mouse <- getAveSobs(region, "mouse")
 	obs.human <- getAveSobs(region, "human")
 	obs.soil <- getAveSobs(region, "soil")

	nseqs <- nochims[nochims$method=="ave", "nseqs"]

	return(c(perfect[1,2], nochims[1,"sobs"], obs.mock, obs.soil, obs.mouse, obs.human, nseqs))
}

otu_table <- t(sapply(regions, getOTUResults))
nseqs <- otu_table[[1,7]]

colnames(otu_table) <- c("no_error","no_chimeras", "mock", "soil", "mouse", "human", "nseqs")

write.table(otu_table, file="data/process/sobs_table.tsv")
