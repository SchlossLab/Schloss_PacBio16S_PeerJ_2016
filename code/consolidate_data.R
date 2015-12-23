# This command will read in the *error.summary file to get a table containing
# various parameters in the mock community data

get_error_summary <- function(error_summary_file){

	read.table(file=error_summary_file, header=T, row.names=1)

}


# This command will read in the data from a *.ccs_stats file to get the fold
# sequence coverage and expected error rate

get_coverage_pred_error <- function(ccs_stats_file, sequences){

	ccs_stats <- read.table(file=ccs_stats_file, row.names=1)[sequences, ]
	colnames(ccs_stats) <- c("coverage", "pred_error")
	return(ccs_stats)

}


# This command will take the *.mismatches file and reformat it to look right

get_mismatches <- function(mismatches_file){

  mismatches <- read.table(file=mismatches_file, row.names=1)
  colnames(mismatches) <- c("fbdiffs", "rbdiffs", "fpdiffs", "rpdiffs")

  mismatches$fbdiffs <- gsub("fbdiffs=(\\d*)\\(.*", "\\1", mismatches$fbdiffs)
  mismatches$rbdiffs <- gsub("rbdiffs=(\\d*)\\(.*", "\\1", mismatches$rbdiffs)
  mismatches$fpdiffs <- gsub("fpdiffs=(\\d*)\\(.*", "\\1", mismatches$fpdiffs)
  mismatches$rpdiffs <- gsub("rpdiffs=(\\d*)\\(.*", "\\1", mismatches$rpdiffs)

  return(mismatches)

}



# We need to calculate the mode to find the best start and end positions in the
# reference alignment

mode <- function(x){
  return(as.numeric(names(sort(table(x), decreasing=T)[1])))
}


# Here we'll process the sequence alignment data to see whether things start/end
# at the correct locations, have the proper homopolymer length, and have any
# ambiguous base calls

get_alignment_quality <- function(alignment_summary_file, sequences){

	alignment_summary <- read.table(file=alignment_summary_file, header=T, row.names=1)

	mode_start <- mode(alignment_summary$start)
	mode_end <- mode(alignment_summary$end)

	#find sequences that start and end at correct location in alignment
	good_coords <- alignment_summary$start == mode_start & alignment_summary$end == mode_end

	#find sequences with less than or equal to 8 nt
	good_homop <- alignment_summary$polymer <= 8

	#find sequences with no ambiguous base calls
	good_ambig <- alignment_summary$ambig == 0

	alignment_quality <- cbind(good_coords, good_homop, good_ambig)
	rownames(alignment_quality) <- rownames(alignment_summary)
	colnames(alignment_quality) <- c("good_coords", "good_homop", "good_ambig")
	return(alignment_quality[sequences,])

}


# Now we want to calculate various parameters related to the quality scores
# for the sequences. We'll calculate...
#	* The mean quality score (w/o accounting for Q->P->Q relationship)
#	* The mean quality score (w accounting for Q->P->Q relationship)
#	* The minimum quality score across the sequence
#	* The expected number of errors across the sequence


# read in the qual file and output scores that overlap with those that were in
# the mock community. scores are outputted as a list with each sequence having a
# vector of quality scores

get_quality_scores <- function(quality_score_file, sequences){

	quality_score_data <- scan(file=quality_score_file, what="", sep="\n", quiet=T)

	seq_names <- quality_score_data[1:length(quality_score_data) %% 2 == 1]
	seq_scores <- quality_score_data[1:length(quality_score_data) %% 2 == 0]

	seq_names <- gsub(">(.*)\t.*", "\\1", seq_names)
	names(seq_scores) <- seq_names

	d <- lapply(seq_scores[sequences], function(x){as.numeric(unlist(strsplit(x, " ")))})

}


# convert Phred scores into probabilities: P = 10^(-Q/10)

get_probabilities <- function(scores){
	10^(-0.1*scores)
}


# convert probabilities into a Phred score: Q = -log10(P)

get_phred <- function(probabilities){
	-10 * log10(probabilities)
}


# get the average Quality score (note this is done *without* converting to
# probabilities first)

get_mean_score <- function(scores){
	mean(scores)
}


# get the average Quality score (note this is done *with* converting to
# probabilities first)

get_log_corrected_mean_score <- function(scores){
	get_phred(mean(get_probabilities(scores)))
}


# get the smallest Quality score

get_min_score <- function(scores){
	min(scores)
}


# sum up the probabilities across the length of the sequence
# http://bioinformatics.oxfordjournals.org/content/31/21/3476

get_expected_errors <- function(scores){
	sum(get_probabilities(scores))
}


# calculate various parameters for each sequence and output as a table

get_quality_score_report <- function(quality_score_file, sequences){

	quality_scores <- get_quality_scores(quality_score_file, sequences)

	mean_score <- sapply(quality_scores, get_mean_score)
	log_corrected_mean_score <- sapply(quality_scores, get_log_corrected_mean_score)
	min_score <- sapply(quality_scores, get_min_score)
	expected_errors <- sapply(quality_scores, get_expected_errors)

	cbind(mean_score, log_corrected_mean_score, min_score, expected_errors)

}


generate_report <- function(report_file) {

	error_summary_file <- gsub("mock.report", "mock.filter.error.summary", report_file)
	error_summary <- get_error_summary(error_summary_file)

	sequence_names <- rownames(error_summary)

	ccs_stats_file <- gsub("mock.report", "ccs_stats", report_file)
	coverage_pred_error <- get_coverage_pred_error(ccs_stats_file, sequence_names)

	mismatches_file <- gsub("mock.report", "mock.mismatches", report_file)
	mismatches <- get_mismatches(mismatches_file)

	alignment_summary_file <- gsub("mock.report", "mock.filter.summary", report_file)
	alignment_quality <- get_alignment_quality(alignment_summary_file, sequence_names)

	qual_file <- gsub("mock.report", "mock.qual", report_file)
	quality_score_report <- get_quality_score_report(qual_file)


	#make sure the tables are all in the same order...
	stopifnot(sum(rownames(error_summary) != rownames(coverage_pred_error)) == 0)
	stopifnot(sum(rownames(error_summary) != rownames(mismatches)) == 0)
	stopifnot(sum(rownames(error_summary) != rownames(alignment_quality)) == 0)
	stopifnot(sum(rownames(error_summary) != rownames(quality_score_report)) == 0)

	batch <- rep(gsub("data/mothur_(.*)/.*", "\\1", report_file), nrow(error_summary))
	

	report_data <- cbind(
		error_summary,
		coverage_pred_error,
		mismatches,
		alignment_quality,
		quality_score_report,
		batch
	)

	write.table(report_data, file=report_file, quote=F)
}

# R -e "source('code/consolidate_data.R');generate_report('data/mothur_june/V1V3.report')"
