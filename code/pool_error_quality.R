pool <- function(file_string){
	files <- unlist(strsplit(file_string, " "))

	pooled_data <- read.table(file=files[1], header=T)
	files <- files[-1]

	qscores <- pooled_data$qscore

	for(f in files){
		temp_data <- read.table(file=files[1], header=T)
		stopifnot(temp_data$qscore == qscores)
	}
	pooled_data$qscore <- qscores
	write.table(pooled_data, "data/process/mock.quality.report", quote=F, row.names=F, sep='\t')
}
