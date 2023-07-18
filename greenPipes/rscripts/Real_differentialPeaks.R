#!/usr/bin/env Rscript
#Real_differentialPeaks.R --files  --dir

#--- install necessary packages
if(!require(R.utils)){
        install.packages("R.utils",repos = "http://cran.us.r-project.org")}
library(R.utils)

if(!require(optparse)){
	install.packages("optparse",repos = "http://cran.us.r-project.org")
    library(optparse)
}

#--------- main

option_list = list(
make_option(c("-f","--files"), type="character",default=NA,help="input file"),
make_option(c("-p","--pvalue"), default=0.0001,help="cutoff pvalue"),
make_option(c("-c","--cvalue"), default=4.0,help="cutoff foldchange"),
make_option(c("-d","--dir"), type="character",default=NA,help="directory"),
make_option(c("-x","--prefix"), type="character",default=NA,help="prefix of output image")
);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)

if(is.na(opt$files)) {
print('Missing input file; run Rscript Real_differentialPeaks.R --help');
} else {
	dir=opt$dir
#	out=paste(dir,'/',out,sep='')
	files=strsplit(opt$file,',')
	files=files[[1]]
#	jpeg(out,unit='in',res=300,height=3.5,width=3.5*length(files))
#	opar<-par(mfrow=c(1,length(files)))
	Totaldata <- list()
	max_fc =0; min_fc =0; max_p =0; min_p =0;
	for(i in c(1:length(files))){
		filename = files[i]
		#filename = paste(dir,'/',files[i],'.TotalPeaks.txt',sep='')
		Totaldata[[i]] <- read.table(filename)
		max_fc=max(max_fc,log2(Totaldata[[i]][,10]))
		min_fc=min(min_fc,log2(Totaldata[[i]][,10]))
		z=max(max_p,-log10(Totaldata[[i]][,11]))
		zz = -log10(Totaldata[[i]][,11])
		if(z!=Inf) {max_p=z} else { max_p = max(zz[!is.infinite(zz)])}
		min_p=min(min_p,-log10(Totaldata[[i]][,11]))
	}
	#---
	for(j in c(1:length(files))){
		data = Totaldata[[j]]
		out=paste(dir,'/',opt$prefix,'.jpeg',sep='')
		jpeg(out,unit='in',res=300,height=3.5,width=3.5)
			for(i in c(1:nrow(data))) {
				ifelse((-log10(data[i,11])>= -log10(opt$pvalue) && log2(data[i,10]) >= log2(opt$cvalue)),'#00008B64', '#FFA50064' ) -> data[i,12]
				}
			plot(log2(data$V10),-log10(data$V11),xlab='log2(Fold change)',ylab='-log10(pvalue)',col=as.vector(data$V12),pch=20,main=files[j],axes=F ,xlim=c(min_fc-4,max_fc+4),ylim=c(min_p-4,max_p+4))
			axis(1)
			axis(2)
			abline(v= c(log2(opt$cvalue),-log2(opt$cvalue)),lty=2)
			abline(h= -log10(opt$pvalue),lty=2)
			box()
			dim(data[(-log10(data[,11])>= -log10(opt$pvalue) & log2(data[,10]) >= log2(opt$cvalue)),])[1]-> value_upperRight
			dim(data[(-log10(data[,11])< -log10(opt$pvalue) & log2(data[,10]) < log2(opt$cvalue)),])[1]-> value_lowerLeft
			legend('topright',legend=value_upperRight,bty='n',cex=1)
			legend('bottomleft',legend=value_lowerLeft,bty='n',cex=1)
		dev.off()
	}
}
