error_summary <- read.table(file="data/process/error_summary.tsv", header=T)

getJitter <- function(x){
	jitter <- 0.15
	runif(6, x-jitter, x+jitter)
}


methods <- c("raw", "basic", "pb", "pred", "coverage", "pb_coverage", "pb_pred",
				"coverage_pred", "all", "precluster")



cairo_pdf(file="submission/figure_2.pdf", width=6.5, height=6.0)

l <- layout(matrix(c(1,2,3), nrow=3), heights=c(1,1,0.6))

par(mar=c(0.5, 7, 0.5, 0.5))
plot(1, type="n", xlim=c(1, length(methods)), ylim=c(0,1), axes=F, xlab="", ylab="")

plotRegionMethod <- function(index, data){
	method <- methods[index]
	subset <- error_summary[error_summary$method==method,]
	regions <- gsub("V", "", subset$region)
	text(x=getJitter(index), y=subset[,data], label=regions, cex=1.25)
}

lapply(1:length(methods), plotRegionMethod, "error")
text(x=4, y=0.8, label="*", cex=3, font=2)

mtext(text="Error rate (%)", side=2, line=4, cex.lab=1.2)
box()
axis(2, las=2, at=seq(0,1, 0.2), label=format(seq(0, 1, 0.2),nsmall=1),
		 cex.axis=1.5)

abline(v=seq(3.5,9.5, 1), col="grey")
abline(v=c(1.5, 2.5, 9.5), col="black", lwd=2)

text(x=1, y=0.95, label="A", cex=2, font=2)

par(mar=c(0.5, 7, 0.5, 0.5))
plot(1, type="n", xlim=c(1, 10), ylim=c(0,100), axes=F,
		 ylab="Reads remaining\nfrom basic criteria (%)", xlab="", cex.lab=1.5)

lapply(1:length(methods), plotRegionMethod, "nseqs_red")


box()
abline(v=seq(3.5,9.5, 1), col="grey")
abline(v=c(1.5, 2.5, 9.5), col="black", lwd=2)

text(x=1, y=95, label="B", cex=2, font=2, las=1)
axis(2, las=2, at=seq(0,100, 25), label=seq(0, 100, 25), cex.axis=1.5)
text(1,0, label="NA")
text(2,0, label="NA")

plot(1, type="n", xlim=c(1, 10), ylim=c(0,1), axes=F, ylab="", xlab="")


methods <- c("raw", "basic", "pb", "pred", "coverage", "pb_coverage", "pb_pred",
				"coverage_pred", "all", "precluster")

text(x=c(1:10)+0.1, y=rep(1, 11), label=c("Initial data", "Basic criteria",
		"\u22641 Mismatch", "Predicted", "Coverage",
		"\u22641 Mismatch & Coverage",
		"\u22641 Mismatch & Predicted",
		"Coverage & Predicted",
		"\u22641 Mismatch,\nCoverage & Predicted", "After Pre-clustering"),
		srt=90, pos=2, cex=1)

layout(1)
dev.off()
