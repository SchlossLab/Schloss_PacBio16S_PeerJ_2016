pretty_region <- c("V4"="V4", "V3V5"="V3-V5", "V1V3"="V1-V3", "V1V5"="V1-V5", "V1V6"="V1-V6", "V1V9"="V1-V9")
regions <- names(pretty_region)

clrs <- rainbow(length(regions))
names(clrs) <- regions

#plot one_offs
one_offs <- read.table(file="data/process/non_random_analysis.tsv", header=T)
z <- prop.table(as.matrix(one_offs), 2)[1:10,]

pdf(file="submission/figure_4.pdf", width=4, height=4)
par(mar=c(4,4,0.5,0.5))
plot(NA, xlim=c(1,10), ylim=c(0,1), xlab="Number of times 1-nt\nvariants were observed", ylab="Percentage of 1-nt variants", axes=F)
for(r in regions){
	points(1:10, z[,r], type="l", col=clrs[r], lwd=3)
}
axis(1, at=1:10, label=1:10)
axis(2, at=seq(0,1,0.2), label=seq(0, 100, 20), las=2)
box()
legend(x=6, y=0.9, legend=pretty_region, lwd=3, col=clrs, cex=0.8)
dev.off()
