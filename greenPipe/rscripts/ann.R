#!/usr/bin/env Rscript
#----- ann.R

if(!require(optparse)){
        install.packages("optparse")}
library(optparse)
#--- 
if(!require(data.table)){
        install.packages("data.table")}
library(data.table)
#--- 
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
#--- 
if(!require(GenomicRanges)){
    print ("GenomicRanges is not avaiable in your system")
    system('sudo apt-get install libcurl4-openssl-dev')
        BiocManager::install("GenomicRanges")
        }

library(GenomicRanges)
#--- 
if(!require(plot.matrix)){
        install.packages("plot.matrix")}
library(plot.matrix)
#--
#--- 
if(!require(UpSetR)){
        install.packages("UpSetR")}
library(UpSetR)
#--- 
option_list = list(
make_option(c("-m", "--myMaxGap"), 
    type="character", 
    default="200",
    help="maximum gap between peaks (comma seperated values)",
    metavar="character"
           ),

make_option(c("-i", "--input"), 
    type="character", 
    default=NA,
    help="Input bedfile", 
    metavar="character"),

make_option(c("-d", "--datalist"), 
    type="character", 
    default=NA,
    help="List of files with which input will be matched", 
    metavar="character"),

make_option(c("-c", "--core"), 
    type="integer", 
    default=NA,
    help="Number of threads", 
    metavar="character"),

make_option(c("-n", "--names"), 
    type="character", 
    default=NA,
    help="Name of the experiments (comma seperated values)", 
    metavar="character"),

make_option(c("-o", "--output"), 
    type="character", 
    default=NA,
    help="Prefix of the output files", 
    metavar="character")

);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
    myMaxGap=opt$myMaxGap
    core=opt$core
    mynames=opt$names
    output=opt$output
    input=opt$input
    datalist=opt$datalist
#--- 
mynames=unlist(
        strsplit(
            mynames, 
            ","))

datalist=unlist(
        strsplit(
            datalist, 
            ","))

myMaxGaps=unlist(
        strsplit(
            as.character(myMaxGap), 
            ","))

if (length(myMaxGaps) == 1){
    myMaxGaps=rep(myMaxGaps,length(datalist))
}

if (length(myMaxGaps) != length(datalist) ){
    stop("number of annSize is not equal to input annotation bed files. you can give single annSize for all input annotation bed fils\n----------------------------------")
    
}

mybed=fread(input,
        header=F,
        sep='\t')

mybed=mybed[,c(1:3)]

colnames(mybed) <- c('Chromosome',
            'Start',
            'End')

mybed=makeGRangesFromDataFrame(mybed)

#--- 
for(i in c(1:length(datalist))) {

            mydata = datalist[i]
            myname = mynames[i]
            mydata=fread(mydata,header=F,sep='\t')
            mydata=mydata[,c(1:3)]

            colnames(mydata) <- c('Chromosome',
                        'Start',
                        'End')


            mydata=makeGRangesFromDataFrame(mydata)
            match=unique(
                data.frame(
                    findOverlaps(
                        mybed, 
                        mydata, 
                        ignore.strand=TRUE, 
                        maxgap = as.integer(myMaxGaps[i])))[,1])

            up=mybed[match,]
            elementMetadata(up) [[ myname ]] <- rep(1,nrow(data.frame(up)))

            down=mybed[-match,]
            elementMetadata(down) [[ myname ]]  <- rep(0,nrow(data.frame(down)))

            total = append(up,down)
            total -> mybed
            print (paste('-> Finished (Myannotations): ', myname,sep=''))

}
#--- 

data.frame(total) -> total
ftable(total[,-c(1:5)]) -> totalOutput
outputTxt = paste(output,".totalOutput.txt",sep="")
outputMat = paste(output,".totalMatrix.txt",sep="")

write.ftable(totalOutput,
    outputTxt,
    quote=F,
    sep="\t")
write.table(total,
    outputMat,
    quote=F,
    sep="\t",
    row.names=F)

#--- 
heatmap(as.matrix(total[,-c(1:5)]),Colv=NA)-> d2
as.matrix(total[,-c(1:5)]) -> xyz
xyz[d2$rowInd,d2$colInd] -> my

jpeg(
        paste(output,'.Distribution.jpeg',sep=''),
        unit='in',
        res=300,
        height=5,
        width=5)

layout(
        matrix(c(1,1,1,1,1,1,2,2),
        nrow=2,
        byrow=F))

plot(my,
        col=c('white','#323896'),
        border=NA,
        las=2,
        axes=F,
        xlab='',
        ylab='',
        main="")

plot(rowSums(my),
        rev(c(1:nrow(my))),
        type='l',
        yaxs = "i",
        axes=F,
        xlab='',
        ylab=''); 

axis(1)

abline(v=1,lty=2,col='darkblue')
abline(v=0,lty=2,col='orange')

system('rm Rplots.pdf')

dev.off() 

jpeg(paste(output,'.UpSetPlot.jpeg',sep=''),unit='in',res=300,height=3.5,width=5)
upset(total[,-c(1:5)], mainbar.y.label = "Overlap", sets.x.label = "Individuals")
dev.off()

