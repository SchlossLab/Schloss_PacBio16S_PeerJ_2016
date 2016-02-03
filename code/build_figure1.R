regions <- c("V1V3", "V1V5", "V1V6", "V1V9", "V3V5", "V4")
pretty_region <- c("V1V3"="V1-V3", "V1V5"="V1-V5", "V1V6"="V1-V6", "V1V9"="V1-V9", "V3V5"="V3-V5", "V4"="V4")

clrs <- rainbow(length(regions))
names(clrs) <- regions



mock_error <- read.table(file="data/process/mock.error.report")
mock_error_nochim <- mock_error[mock_error$numparents==1,]
mock_error_good <- mock_error_nochim[mock_error_nochim$error <= 0.1,]
mock_error_basic <- mock_error_good[mock_error_good$good_coords & mock_error_good$good_homop,]


subset <- mock_error_basic[sample(1:nrow(mock_error_good), 0.05*nrow(mock_error_good)),]


pdf(file="results/figures/figure_1.pdf", width=6.5, height=3.0)
layout(matrix(c(1,2,3), nrow=1), heights=1, widths=c(1,1,1))

#Panel A depicting relationship between predicted and observed error rates
par(mar=c(5, 4, 0.5, 0.5))

plot(subset$error, 1-subset$pred_error, pch=19, cex=0.15, xlim=c(0,0.1), ylim=c(0,0.1), xlab="", ylab="", axes=F)
axis(1, at=seq(0,0.1,0.02), label=seq(0,10,2))
axis(2, at=seq(0,0.1,0.02), label=seq(0,10,2), las=2)
box()
text(x=0.0050, y=0.0975, label="A", cex=2, font=2)
mtext(side=2, line=2, text="Predicted Error (%)", cex=0.7)
mtext(side=1, line=2.5, text="Observed Error (%)", cex=0.7)


#Panel B depicting relationship between errors in barcodes/primers and the
#rest of the sequence.
par(mar=c(5,3,0.5,0))

bcprimer <- apply(mock_error_basic[, c("fbdiffs", "rbdiffs", "fpdiffs", "rpdiffs")], 1, sum)
error_bcprimer <- aggregate(mock_error_basic[,"error"], by=list(bcprimer), function(x){c(mean(x), quantile(x, probs=c(0.025, 0.975)))})

plot(error_bcprimer$x[1:6,1]~error_bcprimer$Group.1[1:6], pch=19,
		 ylab="", xlab="", axes=F, ylim=c(0,0.04))
points(error_bcprimer$x[1:6,1]~error_bcprimer$Group.1[1:6], type="l")
#arrows(x0=error_bcprimer$Group.1[1:6], y0=error_bcprimer$x[1:6,1], y1=error_bcprimer$x[1:6,2], angle=90, length=0.1)
#arrows(x0=error_bcprimer$Group.1[1:6], y0=error_bcprimer$x[1:6,1], y1=error_bcprimer$x[1:6,3], angle=90, length=0.1)

axis(1)
axis(2, at=seq(0,0.04,0.01), label=seq(0,4, 1), las=2)
box()
text(x=0.25, y=0.039, label="B", cex=2, font=2, las=2)
mtext(side=2, line=2, text="Observed Error rate (%)", cex=0.7)
mtext(side=1, line=3, text="Total mismatches to\nbarcodes and primers", cex=0.7)

#Panel C depicting relationship between sequencing coverage and error rate for
#the different regions that were sequence
error_by_coverage <- function(region){
	coverage_region <- mock_error_basic[mock_error_basic$region == region,]
	coverage_error_agg <- aggregate(coverage_region$error, by=list(coverage_region$coverage), mean)

	coverage_error <- rep(NA, max(mock_error_basic$coverage))
	coverage_error[coverage_error_agg$Group.1] <- coverage_error_agg[,"x"]
	coverage_error
}

coverage_fold <- 1:40
coverage_error <- sapply(regions, error_by_coverage)[coverage_fold,]

plotLines <- function(region){
	points(coverage_fold, coverage_error[,region], type="l", lwd=2, col=clrs[region])
}


par(mar=c(5,1,0.5, 2))
plot(1, xlim=range(coverage_fold), ylim=c(0,0.04), type="n", xlab="", ylab="", axes=F)
axis(1)
axis(2, at=seq(0,0.04, 0.01), label=FALSE, tick=TRUE, las=2)

lapply(regions, plotLines)

legend(x=25, y=0.04, legend=pretty_region[regions], lty=1, lwd=2, col=clrs[regions], cex=0.8)
box()
text(x=2, y=0.039, label="C", cex=2, font=2, las=2)
mtext(side=1, line=2.5, text="Coverage", cex=0.7)

dev.off()
