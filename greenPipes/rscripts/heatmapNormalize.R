#!/usr/bin/env Rscript

#----- heatmapNormalize.R <- first create the computeMatrix file then use different normalization
# method to normalized the matrix.

#In the end, python will creadte new gz file and generate heatmap

#zcat Fig4EPeaks.gz | head -n1 > xyz && zcat Fig4EPeaks-SpikeNormalized.gz >> xyz &&
# gzip -c xyz > Fig4EPeaks-SpikeNormalized.gz && plotHeatmap -m Fig4EPeaks-SpikeNormalized.gz
# -o Fig4EPeaks-SpikeNormalized.jpeg --dpi 300 --heatmapWidth 1.5 --refPointLabel "" --regionsLabel
# ""  --heatmapHeight 20 && eog Fig4EPeaks-SpikeNormalized.jpeg
# incorporate the standarrrrrrd deviation in the coverage plot.
#Think: you can not use it find significant differences(?)
if(!require(R.utils)){
        install.packages("R.utils",repos = "http://cran.us.r-project.org")}
library(R.utils)
if(!require(optparse)){
        install.packages("optparse",repos = "http://cran.us.r-project.org")}
library(optparse)
#---
if(!require(data.table)){
        install.packages("data.table",repos = "http://cran.us.r-project.org")}
library(data.table)
#---
if(!require(matrixStats)){
        install.packages("matrixStats",repos = "http://cran.us.r-project.org")}
library(matrixStats)
#---
if(!require(crayon)){
        install.packages("crayon",repos = "http://cran.us.r-project.org")}
library(crayon)
#---
option_list = list(
make_option(c("-i", "--input"),
    type="character",
    default=NA,
    help="input computeMatrix output file (full path)",
    metavar="character"
           ),
make_option(c("-o", "--output"),
    type="character",
    default=NA,
    help="prefix of output file (full path)",
    metavar="character"
           ),
make_option(c("-l", "--sampleLabels"),
    type="character",
    default=NA,
    help="Sample labels",
    metavar="character"
           ),
make_option(c("-c", "--count"),
    type="character",
    default=NA,
    help="comma seperated values of the read count of the human or spikIn",
    metavar="character"
           )
);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)

    input=opt$input

    output=opt$output

    sampleLabels=opt$sampleLabels

	  sampleLabels=unlist(
        strsplit(
            sampleLabels,
            ","))

    normLabels="normalized"

    count=opt$count

	  count=as.numeric(unlist(
        strsplit(
            count,
            ",")))

    if (any(is.na(sampleLabels))) {
		sampleLabels = paste("Samples",seq(1,length(count)), sep = "-")
    }


#---
print("Install pv and pigz, else normalized files will not generate.")
#-----------------
if (!file.exists(input)) {
	print(paste(input, " :input file does not exist!", sep =""))
}

print(paste(
	"[] Total number of the samples in the input file are: ",
	length(count),"\n",
	sep=""
	))

#_________________
d = fread(input,skip=1)
print ("[] Read input file")

#_________________
h=data.frame(d[,-c(1:6)]);
c=0
bins = ncol(h) / length(count)

for(i in seq(1,ncol(d)-bins,bins)) {
	c=c+1;
	j = i+(bins - 1);
	h[,c(i:j)]=h[,c(i:j)]/count[c];
	print (paste("#",i,j,count[c],sep="-"))
}
print ("[] Dataset normalized")

#_________________
colMeans(as.matrix(h)) -> hMean

for(i in seq(1,ncol(d)-bins,bins)) {
	j = i+(bins - 1)
	hMean[i:j] = hMean[i:j] - min(hMean[i:j])
}

#__________________
jpeg(paste(output,'-',"normalized-curve",".jpeg",sep=""),
	unit = "in",
	res = 300,
	height = 3.5,
	width = max(3.5,(1*length(count)))
	)
par(mar=c(10,4,4,2))
plot(hMean, type = 'n', ylab = "Coverage", xlab = "", axes =F )

xx = 0
for(r in seq(1,
				bins * (length(count)),
				bins)
			)
		{
      if (xx == 0) {xx = 1} else {xx = 0}
      if (xx == 1) {
			rect(r,-500000000000000000, (r+bins)-1, 5000000000000, col = "cyan", border=FALSE)
      }
			}

lines(hMean, col = "black")
box() ; axis(2)

hAt = seq(
		round(bins/2),
		bins * (length(count)),
		bins
		)
axis(1, at = hAt, labels = sampleLabels, las = 2)
dev.off()

#_________________
oo = paste(output,'-',normLabels,".gz", sep = "")
fwrite(data.frame(d[,c(1:6)],h),
	paste(output,'-',normLabels,".gz", sep = ""),
	quote=F,
	sep="\t",
	row.names=F,
	col.names=F
	)

system(paste("pv ", input ," | pigz -dc | head -n1 | pigz -c > nnnnnn  && cat nnnnnn ", oo, " > nz && mv nz ", oo, " && rm nnnnnn", sep =""))
#_________________
fwrite(data.frame(sample=rep(sampleLabels, each = bins),coverage=hMean),
	paste(output,'-',normLabels,"-curve.gz", sep = ""),
	quote=F,
	sep="\t",
	row.names=F,
	col.names=F
	)
