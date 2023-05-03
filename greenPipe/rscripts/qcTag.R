#!/usr/bin/env Rscript

if(!require(optparse)){
        install.packages("optparse")
    library(optparse)
}

if(!require(data.table)){
        install.packages("data.table")
    library(data.table)
}


option_list = list(
    make_option(c("-a", "--file1"),
                type="character",
                default=NA,
                help="tagCountDistribution: experiment",
                metavar="character"),

    make_option(c("-b", "--file2"),
                type="character",
                default=NA,
                help="tagCountDistribution: control",
                metavar="character"),

    make_option(c("-c", "--file3"),
                type="character",
                default=NA,
                help="tagAutocorrelation: experiment",
                metavar="character"),

    make_option(c("-d", "--file4"),
                type="character",
                default=NA,
                help="tagAutocorrelation: control",
                metavar="character"),

    make_option(c("-e", "--file5"),
                type="character",
                default=NA,
                help="FragmentLength: experiment",
                metavar="character"),

    make_option(c("-f", "--file6"),
                type="character",
                default=NA,
                help="FragmentLength: control",
                metavar="character"),

    make_option(c("-o", "--output"),
                type="character",
                default=NA,
                help="Name and location of output JPEG file",
                metavar="character")
);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
    file1=opt$file1
    file2=opt$file2
    file3=opt$file3
    file4=opt$file4
    file5=opt$file5
    file6=opt$file6
    output=opt$output

tgCnt_expr=read.table(file1,skip=1)
tgCnt_ctrl=read.table(file2,skip=1)

tgCor_expr=read.table(file3,skip=1)
tgCor_ctrl=read.table(file4,skip=1)

if ( !is.na(file5) || !is.na(file6) ){

    tgFr_expr=fread(file5)
    tgFr_expr=data.frame(tgFr_expr)
    tgFr_ctrl=fread(file6)
    tgFr_ctrl=data.frame(tgFr_ctrl)
    tgFr_expr$diff = tgFr_expr[,3]-tgFr_expr[,2]
    tgFr_ctrl$diff = tgFr_ctrl[,3]-tgFr_ctrl[,2]
    density(tgFr_expr$diff) -> tgFr_expr_density
    density(tgFr_ctrl$diff) -> tgFr_ctrl_density

    max_y=max(c(tgFr_expr_density$y,tgFr_expr_density$y))
    max_x=max(c(tgFr_ctrl_density$x,tgFr_ctrl_density$x))
    min_y=min(c(tgFr_expr_density$y,tgFr_expr_density$y))
    min_x=min(c(tgFr_ctrl_density$x,tgFr_ctrl_density$x))
}

jpeg(
    output,
    unit='in',
    res=300,
    height=2.5*3,
    width=2.5*2)

opar<-par(
    mfrow=c(3,2)
)

plot(tgCnt_expr,
     ylim=c(0,1),
     xlab="Tag position",
     ylab="Frequency of tags per position",
     main='Experiment: frequency');
grid()

plot(tgCnt_ctrl,
     ylim=c(0,1),
     xlab="Tag position",
     ylab="Frequency of tags per position",
     main='Control: frequency');
grid()

#---
plot(tgCor_expr[,1],
     tgCor_expr[,2],
     xlab="Relative distance between reads (bps)",
     ylab="Total reads",main='Experiment: autocorrelation',
     type='l');
lines(tgCor_expr[,1],
      tgCor_expr[,3],
      col='grey72')
grid()

plot(tgCor_ctrl[,1],
     tgCor_ctrl[,2],
     xlab="Relative distance between reads (bps)",
     ylab="Total reads",
     main='Control: autocorrelation',
     type='l');

lines(tgCor_ctrl[,1],
      tgCor_ctrl[,3],
      col='grey72')
grid()
#---
if (!is.na(file5) || !is.na(file6)){

    plot(tgFr_expr_density$x,
         tgFr_expr_density$y,
         xlab="Fragment length in base-pairs",
         ylab="Density",
         main='Density of fragment length',
         xlim=c(min_x,max_x),
         ylim=c(min_y,max_y),
         type='l')

    lines(tgFr_ctrl_density$x,
          tgFr_ctrl_density$y,
          col='grey72')

    legend('topright',
           col=c('black','grey72'),
           lty=2,
           bty='n',
           legend=c('Experiment','Control')
          )
    grid()
}

dev.off()
