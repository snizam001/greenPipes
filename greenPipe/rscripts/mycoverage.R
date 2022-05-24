#--------------------------------
#---- mycoverage.R
#--------------------------------


#!/usr/bin/env Rscript

if(!require(optparse)){
        install.packages("optparse")}
library(optparse)
#--- 
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
#--- 
if(!require(GenomicAlignments)){
	print ("GenomicAlignments is not avaiable in your system")
        BiocManager::install("GenomicAlignments")
        }

library(GenomicAlignments)
#--- 
if(!require(rtracklayer)){
	print ("rtracklayer is not avaiable in your system")
        BiocManager::install("rtracklayer")
        }

library(rtracklayer)
#--- 
if(!require(BSgenome)){
	print ("BSgenome is not avaiable in your system")
        BiocManager::install("BSgenome")
        }
#----
if(!require(plot.matrix)){
        install.packages("plot.matrix")}
library(plot.matrix)

#--- 
option_list = list(
make_option(c("-m", "--myMaxGap"), type="integer", default=200,help="maximum gap between reads and region", metavar="character"),
make_option(c("-i", "--inputregions"), type="character", default=NA,help="Input regions in bed format", metavar="character"),
make_option(c("-r", "--regionName"), type="character", default=NA,help="name of the region", metavar="character"),
make_option(c("-b", "--bamList"), type="character", default=NA,help="List of the bamfiles", metavar="character"),
make_option(c("-n", "--names"), type="character", default=NA,help="Name of the experiments or bamfiles", metavar="character"),
make_option(c("-o", "--output"), type="character", default=NA,help="Prefix of the output files", metavar="character"),
make_option(c("-g", "--genome"), type="character", default=NA,help="genome version hg38/hg19", metavar="character"),
make_option(c("-a", "--alignment"), type="character", default=NA,help="pair or single, comma seperated values", metavar="character")

);



parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
distance=opt$myMaxGap
myregions=opt$inputregions
l_names=opt$names
l_bamfiles=opt$bamList
output=opt$output
alignment=opt$alignment
regionName=opt$regionName
genome=opt$genome
#---
if(genome=='hg38') {
	if(!require(BSgenome.Hsapiens.UCSC.hg38)){
		print ("BSgenome.Hsapiens.UCSC.hg38 is not avaiable in your system")
       	 BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
        	}

	library(BSgenome.Hsapiens.UCSC.hg38)
}
if(genome=='hg19'){
	if(!require(BSgenome.Hsapiens.UCSC.hg19)){
		print ("BSgenome.Hsapiens.UCSC.hg19 is not avaiable in your system")
       	 BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
        	}

	library(BSgenome.Hsapiens.UCSC.hg19)
}
#---
mynames=unlist(strsplit(l_bamfiles, ","))
l_names=unlist(strsplit(l_names, ","))
alignment=unlist(strsplit(alignment, ","))
myregions=unlist(strsplit(myregions, ","))
regionName=unlist(strsplit(regionName, ","))
#---
myoutputs=list()
for (j in c(1:length(myregions))){ 
	myregion=myregions[j]
	bedfile <- import.bed(myregion)
	start(bedfile) - distance -> start(bedfile) 
	end(bedfile) + distance -> end(bedfile) 
	if(genome=='hg19'){ seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, bedfile)}
	if(genome=='hg38'){ seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, bedfile)}
	elementMetadata(bedfile) [[ 'GC_percentage' ]] <- as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE))
	#--
	for (i in c(1:length(l_names))) {
		if(alignment[i]!="pair") {
			mybamfile=mynames[i]
			bamfile <- readGAlignments(mybamfile) 
			elementMetadata(bedfile) [[ l_names[i] ]] <- (10^6*(countOverlaps(bedfile,bamfile)/length(bamfile)))/bedfile$GC_percentage
			print(paste('->>>>>> Finished: ',l_names[i] ,sep=''))
#			print(bedfile)	
		}
		if(alignment[i]=='pair') {
			mybamfile=mynames[i]
			bamfile <- readGAlignmentPairs(mybamfile) 
			elementMetadata(bedfile) [[ l_names[i] ]] <- (10^6*(countOverlaps(bedfile,bamfile)/length(bamfile)))/bedfile$GC_percentage
			print(paste('->>>>>> Finished: ',l_names[i] ,sep=''))
#			print(bedfile)
		}
	}
	myoutputs[[j]] = data.frame(bedfile)
}


jpeg(paste(output,"-mycoverage.jpeg",sep=""),unit="in",res=300,height=5*length(myregions),width=5)
opar=par(mfrow=c(length(myregions),1))
par(mar=c(10,5,5,10))

write.table("",paste(output,"-mycoverage.txt",sep=""),row.names=F,col.names=F,quote=F)

for (j in c(1:length(myregions))){ 
	myregions_matrix=as.matrix(myoutputs[[j]][, -which(names(myoutputs[[j]]) %in% c("seqnames","start", 'end', 'width', 'strand', 'name', 'score','experiment', 'GC_percentage')) ])
	myregions_distance=dist(myregions_matrix)
	hclust(myregions_distance) -> myregions_hclust
	myregions_matrix=myregions_matrix[myregions_hclust$order,]
	plot(log2(myregions_matrix+1),border=NA,las=2,axis.row=NULL,xlab='',ylab='',col=heat.colors(150), main = regionName[j], cex.axis = 0.6 );
	myoutputs[[j]]$experiment <- rep(regionName[j],nrow(myoutputs[[j]]))
	write.table(myoutputs[[j]],paste(output,"-mycoverage.txt",sep=""),row.names=F,col.names=F,quote=F, append=T)
}


dev.off() 


