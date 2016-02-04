data <- read.table(file="data/process/taxonomy_depth_analysis.tsv", header=T)

rdp <- 100 * cbind("Mock"=data[data$db=="pds" & data$dataset=="mock", "genus_species"],
				 "Human"=data[data$db=="pds" & data$dataset=="human", "genus_species"],
				 "Mouse"=data[data$db=="pds" & data$dataset=="mouse", "genus_species"],
				 "Soil"=data[data$db=="pds" & data$dataset=="soil", "genus_species"])
rownames(rdp) <- c("V4", "V3-V5", "V1-V3", "V1-V5", "V1-V6", "V1-V9")

gg <- 100 * cbind("Mock"=data[data$db=="gg" & data$dataset=="mock", "genus_species"],
				"Human"=data[data$db=="gg" & data$dataset=="human", "genus_species"],
				"Mouse"=data[data$db=="gg" & data$dataset=="mouse", "genus_species"],
				"Soil"=data[data$db=="gg" & data$dataset=="soil", "genus_species"])
rownames(gg) <- c("V4", "V3-V5", "V1-V3", "V1-V5", "V1-V6", "V1-V9")

silva <- 100 * cbind("Mock"=data[data$db=="nr" & data$dataset=="mock", "genus_species"],
					 "Human"=data[data$db=="nr" & data$dataset=="human", "genus_species"],
					 "Mouse"=data[data$db=="nr" & data$dataset=="mouse", "genus_species"],
					 "Soil"=data[data$db=="nr" & data$dataset=="soil", "genus_species"])
rownames(silva) <- c("V4", "V3-V5", "V1-V3", "V1-V5", "V1-V6", "V1-V9")

gg.sp <- 100 * cbind("Mock"=data[data$db=="gg" & data$dataset=="mock", "X7"],
					 "Human"=data[data$db=="gg" & data$dataset=="human", "X7"],
					 "Mouse"=data[data$db=="gg" & data$dataset=="mouse", "X7"],
					 "Soil"=data[data$db=="gg" & data$dataset=="soil", "X7"])
rownames(gg.sp) <- c("V4", "V3-V5", "V1-V3", "V1-V5", "V1-V6", "V1-V9")



pdf(file="submission/figure_3.pdf", width=5, height=7)


par(mar=c(5, 3, 0.5, 0.5))
dotchart(gg, xlim=c(0,100), col="black", xlab="Unique reads that classified\nto genus or species level (%)", pch=19)

points(x=rdp[,"Mock"], y=25:30, pch=15, col="black")
points(x=rdp[,"Human"], y=17:22, pch=15, col="black")
points(x=rdp[,"Mouse"], y=9:14, pch=15, col="black")
points(x=rdp[,"Soil"], y=1:6, pch=15, col="black")

points(x=silva[,"Mock"], y=25:30, pch=17, col="black")
points(x=silva[,"Human"], y=17:22, pch=17, col="black")
points(x=silva[,"Mouse"], y=9:14, pch=17, col="black")
points(x=silva[,"Soil"], y=1:6, pch=17, col="black")

points(x=gg.sp[,"Mock"], y=25:30, pch=21, bg="gray")
points(x=gg.sp[,"Human"], y=17:22, pch=21, bg="gray")
points(x=gg.sp[,"Mouse"], y=9:14, pch=21, bg="gray")
points(x=gg.sp[,"Soil"], y=1:6, pch=21, bg="gray")

legend(x=55, y=8, legend=c("RDP (gen.)", "SILVA (gen.)", "greengenes (gen.+sp.)", "greengenes (sp.)"), pch=c(15, 17, 19, 21), col="black", pt.bg=c("black", "black", "black", "gray"), bg="white", cex=0.8)

dev.off()
