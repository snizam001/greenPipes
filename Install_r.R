dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))

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

for(i in c(1:length(packageInstall))){
if(!require(packageInstall[i],character.only = TRUE)){
	install.packages(packageInstall[i],repos = "http://cran.us.r-project.org")
	}
library(packageInstall[i],character.only = TRUE)
}

library(BiocManager)
BiocManager::install(c("GenomicRanges","GenomicAlignments","rtracklayer","BSgenome"))
