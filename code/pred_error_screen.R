get_sequences <- function(fasta_file_name){
	fasta <- scan(fasta_file_name, quiet=T, what="character", sep="\n")
	sequences <- fasta[(1:length(fasta) %% 2)==0]
	seq_names <- fasta[(1:length(fasta) %% 2)==1]
	names(sequences) <- gsub(">(.*ccs).*", "\\1", seq_names)
	sequences
}

get_ccs_data <- function(ccs_file_name, threshold = 0.999){
	ccs_data <- read.table(ccs_file_name)
	ccs_data[ccs_data$V3>=threshold, "V1"]
}

get_good_seqs <- function(sequences, good_accnos){
	good_seqs <- names(sequences) %in% good_accnos
	sequences[good_seqs]
}

make_fasta <- function(sequences){
	names(sequences) <- gsub("^(.*)", ">\\1\n", names(sequences))
	sequences
}

pred_error_screen <- function(fasta_file_name, ccs_file_name){
	output_file_name <- gsub("fasta", "screen.fasta", fasta_file_name)

	seq_data <- get_sequences(fasta_file_name)
	ccs_data <- get_ccs_data(ccs_file_name)
	good_seq_data <- get_good_seqs(seq_data, ccs_data)
	output_fasta <- make_fasta(good_seq_data)

	write.table(output_fasta, output_file_name, sep="", quote=F, col.names=F)
}

