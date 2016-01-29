regions <- c("V1V3", "V1V5", "V1V6", "V1V9", "V3V5", "V4")

mock_error <- read.table(file="data/process/mock.error.report")
mock_error_nochim <- mock_error[mock_error$numparents==1,]
mock_error_good <- mock_error_nochim[mock_error_nochim$error <= 0.10,]




## Basic error correction analysis
get_raw_error <- function(region){
	region_data <- mock_error_good[mock_error_good$region == region,]
	return(c(100*sum(region_data$mismatches)/sum(region_data$total), nrow(region_data)))
}
raw_rates <- sapply(regions, get_raw_error)
rownames(raw_rates) <- c("error", "nseqs")


get_basic_error <- function(region){
	region_data <- mock_error_good[mock_error_good$region == region & mock_error_good$good_coords & mock_error_good$good_homop,]
	return(c(100*sum(region_data$mismatches)/sum(region_data$total), nrow(region_data)))
}
basic_rates <- sapply(regions, get_basic_error)
rownames(basic_rates) <- c("error", "nseqs")



mock_error_basic <- mock_error_good[mock_error_good$good_coords & mock_error_good$good_homop,]






#primer and barcode analysis
mock_error_basic$pbdiffs <- mock_error_basic$fbdiffs + mock_error_basic$rbdiffs + mock_error_basic$fpdiffs + mock_error_basic$rpdiffs



# want to get the error rate and fraction of sequences for each region with a
# threshold of mismatches or fewer

pb_by_region_mm <- aggregate(mock_error_basic$error, by=list(mock_error_basic$pbdiffs, mock_error_basic$region), function(x){c(mean=mean(100*x), nseqs=length(x))})

aggregate_pb_data <- function(region){
	pb_region <- pb_by_region_mm[pb_by_region_mm$Group.2 == region,]
	product <- pb_region$x[,"mean"] * pb_region$x[,"nseqs"]
	pb_region_nseqs <- cumsum(pb_region$x[,"nseqs"])
	pb_region_error <- cumsum(product) / pb_region_nseqs
	c(error=pb_region_error, nseqs=pb_region_nseqs)
}
pb_data <- sapply(regions, aggregate_pb_data)

pb_error <- pb_data[grep("error", rownames(pb_data)),]
rownames(pb_error) <- 0:6

pb_nseqs <- pb_data[grep("nseqs", rownames(pb_data)),]
rownames(pb_nseqs) <- 0:6

pb_error_reduction_1 <- 100 * (1 - pb_error["1",]/basic_rates["error",])
pb_nseqs_reduction_1 <- 100 * (1 - pb_nseqs["1",]/basic_rates["nseqs",])




#coverage_analysis - need to get the error rate for coverage values between 3 and 40 for each region

#error_by_coverage <- function(region){
#	coverage_region <- mock_error_basic[mock_error_basic$region == region,]
#	coverage_error_agg <- aggregate(100 * coverage_region$error, by=list(coverage_region$coverage), mean)
#
#	coverage_error <- rep(NA, max(mock_error_basic$coverage))
#	coverage_error[coverage_error_agg$Group.1] <- coverage_error_agg[,"x"]
#	coverage_error
#}
#
#coverage_error <- sapply(regions, error_by_coverage)
#seems to plateau right around 10


error_at_below_coverage_cutoff <- function(region, cutoff=10){
	region_data <- mock_error_basic[mock_error_basic$region == region,]
	deep <- region_data$coverage >= cutoff
	coverage_data <- region_data[deep,]
	error <- 100 * sum(coverage_data$mismatches)/sum(coverage_data$total)
	c(error=error, nseqs=nrow(coverage_data))
}

coverage_error <- sapply(regions, error_at_below_coverage_cutoff)

coverage_error_reduction <- 100 * (1 - coverage_error["error",]/basic_rates["error",])
coverage_nseqs_reduction <- 100 * (1 - coverage_error["nseqs",]/basic_rates["nseqs",])



#screen by predicted error rate
error_at_above_pred_error_cutoff <- function(region, threshold=0.999){
	region_data <- mock_error_basic[mock_error_basic$region == region,]
	above_threshold <- region_data[region_data$pred_error >= threshold,]
	c(error=100*sum(above_threshold$mismatches)/sum(above_threshold$total), nseqs=nrow(above_threshold))
}

pred_error_data <- sapply(regions, error_at_above_pred_error_cutoff, threshold=0.999)

pred_error_reduction <- 100 * (1 - pred_error_data["error",]/basic_rates["error",])
pred_nseqs_reduction <- 100 * (1 - pred_error_data["nseqs",]/basic_rates["nseqs",])




#oligos & coverage analysis
get_oligos_coverage <- function(region){
	region_data <- mock_error_basic[mock_error_basic$region == region,]
	good1 <- region_data$pbdiffs <= 1 & region_data$coverage >= 10
	region_data_good1 <- region_data[good1,]
	c(	error = 100*mean(region_data_good1[,"error"]), nseqs = nrow(region_data_good1))
}

oligos_coverage_error <- matrix(unlist(lapply(regions, get_oligos_coverage)), ncol=2, byrow=T)
rownames(oligos_coverage_error) <- regions
colnames(oligos_coverage_error) <- c("error", "nseqs")

#     0.10.error 0.10.frac 1.10.error 1.10.frac
#V1V3  0.4515344 0.5898447  0.4632910 0.7968639
#V1V5  0.4187311 0.5557715  0.4270938 0.7384810
#V1V6  0.4399351 0.5410259  0.4460726 0.7167325
#V1V9  0.4317629 0.4444776  0.4420337 0.5752849
#V3V5  0.6421004 0.6854427  0.6555240 0.8102832
#V4    0.6825316 0.5512332  0.7023460 0.8245533

oligos_coverage_error_reduction <- 100 * (1 - oligos_coverage_error[,"error"]/basic_rates["error",])
oligos_coverage_nseqs_reduction <- 100 * (1 - oligos_coverage_error[,"nseqs"]/basic_rates["nseqs",])



#Oligos & QScore Analysis
get_oligos_pred_error <- function(region){
	region_data <- mock_error_basic[mock_error_basic$region == region,]
	good1 <- region_data$pbdiffs <= 1 & region_data$pred_error >= 0.999
	region_data_good1 <- region_data[good1,]
	c(error = 100 * mean(region_data_good1[,"error"]), nseqs = nrow(region_data_good1) )
}

oligos_pred_error <- matrix(unlist(lapply(regions, get_oligos_pred_error)), ncol=2, byrow=T)
rownames(oligos_pred_error) <- regions
colnames(oligos_pred_error) <- c("error", "nseqs")

#     0.999.error 0.999.frac 1.999.error 1.999.frac
#V1V3   0.2806049  0.4749822   0.2806049  0.4749822
#V1V5   0.2177965  0.4266279   0.2177965  0.4266279
#V1V6   0.2168503  0.4296550   0.2168503  0.4296550
#V1V9   0.2177073  0.3438720   0.2177073  0.3438720
#V3V5   0.2229757  0.5131411   0.2229757  0.5131411
#V4     0.3485654  0.3961328   0.3485654  0.3961328

oligos_pred_error_reduction <- 100 * (1 - oligos_pred_error[,"error"]/basic_rates["error",])
oligos_pred_nseqs_reduction <- 100 * (1 - oligos_pred_error[,"nseqs"]/basic_rates["nseqs",])


# coverage & pred error analysis
get_coverage_pred <- function(region){
	region_data <- mock_error_basic[mock_error_basic$region == region,]
	good <- region_data$coverage >= 10 & region_data$pred_error >= 0.999
	region_data_good <- region_data[good,]
	c(100 * mean(region_data_good[,"error"]), nrow(region_data_good) )
}

coverage_pred_error <- matrix(unlist(lapply(regions, get_coverage_pred)), ncol=2, byrow=T)
rownames(coverage_pred_error) <- regions
colnames(coverage_pred_error) <- c("error", "nseqs")

#    10.999.error 10.999.frac
#V1V3   0.2843963    65.86436
#V1V5   0.2188792    56.55108
#V1V6   0.2178326    57.33756
#V1V9   0.2191284    42.66675
#V3V5   0.2236951    60.74254
#V4     0.3527900    61.79213

coverage_pred_error_reduction <- 100 * (1 - coverage_pred_error[,"error"]/basic_rates["error",])
coverage_pred_nseqs_reduction <- 100 * (1 - coverage_pred_error[,"nseqs"]/basic_rates["nseqs",])


# combine all filters analysis
get_all_filters <- function(region){
	region_data <- mock_error_basic[mock_error_basic$region == region,]
	good1 <- region_data$pbdiffs <= 1 & region_data$pred_error >= 0.999 & region_data$coverage >= 10
	region_data_good1 <- region_data[good1,]
	c(error = 100 * mean(region_data_good1[,"error"]), nseqs =nrow(region_data_good1))
}

all_filters_error <- matrix(unlist(lapply(regions, get_all_filters)), ncol=2, byrow=T)
rownames(all_filters_error) <- regions
colnames(all_filters_error) <- c("error", "nseqs")

#     0.10.999.error 0.10.999.frac 1.10.999.error 1.10.999.frac
#V1V3      0.2808218     0.4696170      0.2840387     0.6183602
#V1V5      0.2184805     0.4192252      0.2182989     0.5374066
#V1V6      0.2170989     0.4182253      0.2174090     0.5421862
#V1V9      0.2176989     0.3213802      0.2191334     0.4016057
#V3V5      0.2241187     0.5038275      0.2240779     0.5841414
#V4        0.3492870     0.3932393      0.3526755     0.5665274

all_filters_error_reduction <- 100 * (1 - all_filters_error[,"error"]/basic_rates["error",])
all_filters_nseqs_reduction <- 100 * (1 - all_filters_error[,"nseqs"]/basic_rates["nseqs",])




#Unique and pre-cluster analysis
get_weighted_error_rate <- function(file_name){
	error_summary <- read.table(file=file_name, header=T, stringsAsFactors=FALSE)
	nochim <- error_summary[error_summary$numparents==1,]
    c(error = 100 * sum(nochim$weight * nochim$mismatches) / sum(nochim$weight * nochim$total), nseqs = sum(nochim$weight))
}

precluster_file_names <- paste0("data/mothur_pool/", regions,".mock.screen.unique.good.filter.unique.precluster.error.summary")

precluster_error <- sapply(precluster_file_names, get_weighted_error_rate)
colnames(precluster_error) <- gsub("data/mothur_pool/(.*).mock.screen.unique.good.filter.unique.precluster.error.summary", "\\1", colnames(precluster_error))

precluster_error_reduction <- 100 * (1 - precluster_error["error",]/basic_rates["error",])
precluster_nseqs_reduction <- 100 * (1 - precluster_error["nseqs",]/basic_rates["nseqs",])







#pool results
base_error_reduction <- rbind(	raw = NA,
								basic = NA,
								pb = pb_error_reduction_1,
								pred = pred_error_reduction,
								coverage = coverage_error_reduction,
								coverage_pred = coverage_pred_error_reduction,
								pb_coverage = oligos_coverage_error_reduction,
								pb_pred = oligos_pred_error_reduction,
								all = all_filters_error_reduction,
								precluster = precluster_error_reduction
							)

base_nseqs_reduction <- rbind(	raw = NA,
								basic = NA,
								pb = pb_nseqs_reduction_1,
								pred = pred_nseqs_reduction,
								coverage = coverage_nseqs_reduction,
								coverage_pred = coverage_pred_nseqs_reduction,
								pb_coverage = oligos_coverage_nseqs_reduction,
								pb_pred = oligos_pred_nseqs_reduction,
								all = all_filters_nseqs_reduction,
								precluster = precluster_nseqs_reduction
							)

error <- rbind(	raw = raw_rates["error",],
				basic = basic_rates["error",],
				pb = pb_error["1",],
				pred = pred_error_data["error",],
				coverage = coverage_error["error",],
				coverage_pred = coverage_pred_error[,"error"],
				pb_coverage = oligos_coverage_error[,"error"],
				pb_pred = oligos_pred_error[,"error"],
				all = all_filters_error[,"error"],
				precluster = precluster_error["error",]
			)

nseqs <- rbind(	raw = raw_rates["nseqs",],
				basic = basic_rates["nseqs",],
				pb = pb_nseqs["1",],
				pred = pred_error_data["nseqs",],
				coverage = coverage_error["nseqs",],
				coverage_pred = coverage_pred_error[,"nseqs"],
				pb_coverage = oligos_coverage_error[,"nseqs"],
				pb_pred = oligos_pred_error[,"nseqs"],
				all = all_filters_error[, "nseqs"],
				precluster = precluster_error["nseqs",]
			)

summary_data <- data.frame(
				cbind(	error = as.vector(error),
						nseqs = as.vector(nseqs),
						error_red = as.vector(base_error_reduction),
						nseqs_red = as.vector(base_nseqs_reduction)
					))
summary_data$region <- rep(colnames(error), each=nrow(error))
summary_data$method <- rep(rownames(error), ncol(error))

write.table(summary_data, file="data/process/error_summary.tsv", quote=F, row.names=F, sep='\t')
