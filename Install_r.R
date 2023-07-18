packageInstall = c("R.utils", 
	"optparse",
	"data.table",
	"BiocManager",
	"plot.matrix",
	"UpSetR",
	"svMisc",
	"matrixStats",
	"crayon",
	"VennDiagram",
	"gridExtra")

install.packages(packageInstall,repos = "http://cran.us.r-project.org")

library(BiocManager)
BiocManager::install(c("GenomicRanges","GenomicAlignments","rtracklayer","BSgenome"))
