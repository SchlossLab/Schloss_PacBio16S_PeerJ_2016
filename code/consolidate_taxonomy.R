#	read in a bunch of taxonomy files and figure out how deeply each dataset
#	was classified. outputs a summary file into data/process/

count_levels <- function(tax_string){
	length(unlist(strsplit(tax_string, split=";")))
}

count_tax_levels <- function(combo){
	taxonomy_file_name <- paste0("data/mothur_pool/",combo["region"],".",combo["dataset"],".screen.unique.good.filter.unique.precluster.pick.", combo["db"], ".wang.taxonomy")
	taxonomy <- read.table(file=taxonomy_file_name, stringsAsFactors=FALSE)
	classified <- gsub("unclassified;", "", taxonomy$V2)
	depth_count <- table(factor(sapply(classified, count_levels), levels=0:7))
	depth_count / sum(depth_count)
}

tax_levels <- 0:7
regions <- c("V4", "V1V3", "V3V5", "V1V5", "V1V6", "V1V9")
databases <- c("pds", "gg", "nr")
datasets <- c("mock", "mouse", "human", "soil")

combos <- cbind(region=rep(regions, each=length(databases)*length(datasets)),
				dataset=rep(datasets, length(regions)*length(databases)),
				db=rep(rep(databases, each=length(datasets)), length(regions)))

tax_level_counts <- t(sapply(1:nrow(combos), function(x){count_tax_levels(combos[x,])}))
genus_species <- tax_level_counts[,"6"] + tax_level_counts[,"7"]
composite <- cbind(combos, tax_level_counts, genus_species)

write.table(file="data/process/taxonomy.depth.analysis", composite, quote=F, row.names=F)
