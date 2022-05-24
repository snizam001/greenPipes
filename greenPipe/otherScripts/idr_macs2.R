#' plot irreproducible discovery rate
#'
#' Make plots from the output of the function \code{calculateIDR}
#'
#' The output of \code{calculateIDR()} is a list of tuples of ChIP samples.
#' \code{plotIDR()} draws plots which compare the two members of a tuple.
#' For every tuple six plots are generated:
#' - number of peaks in common as a function of the number of _all_ significant peaks
#' - slope of the previous plot
#' - number of peaks in common as a function of the number of _matched_ significant peaks
#' - slope of the previous plot
#' - IDR as a function of the number of significant peaks
#' - scatter plot of ranked peaks of the ChIP samples
#' Refer to Qunhua Li et al (2011) for an explanation on how to interpret
#' the plots.
#'
#' Ref: Qunhua Li, James B. Brown, Haiyan Huang, and Peter J. Bickel: Measuring reproducibility of high-throughput experiments.
#'      Ann Appl Stat. 2011 October 13 
#'
#' @export
#' @param chipTuple one item of the return value of \calc{calculateIDR}
#' @param idrCutoff peaks with an IDR greater than this cut-off are colored red in the scatter plot
#' @examples
#' chipTuples <- calculateIDR(c("IP1.bam", "IP2.bam"), c("input1.bam", "input2.bam"),
#'                             "Hsapiens", "UCSC", "hg19")
#' for (chipTuple in chipTuples) {
#'     pdf(paste(basename(chipTuple$rep1), "_VS_", basename(chipTuple$rep2), ".pdf", sep=""), paper="a4r", width=11, height=8.5)
#'     plotIDR(chipTuple)
#'     dev.off()
#' }
plotIDR <- function(chipTuple, idrCutoff=0.01) {
	df.txt <- 10

	ez <- get.ez.tt.all(chipTuple$em.output, chipTuple$uri.output$data12.enrich$merge1,
		chipTuple$uri.output$data12.enrich$merge2) # reverse =T for error rate

	# URI for all peaks
	uri.all.peaks <- chipTuple$uri.output$uri.n
	# URI for matched peaks
	uri.match <- get.uri.matched(chipTuple$em.output$data.pruned, df=df.txt)
	uri.matched.peaks <- uri.match$uri.n

	# calculate IDR values
	idrs <- 1-chipTuple$em.output$em.fit$e.z
	idrs[order(idrs)] <- cumsum(idrs[order(idrs)])/(1:length(idrs))
	peaksBelowIDRCutoff <- idrs <= idrCutoff

	############# plot and report output
	# plot correspondence curve for each pair,
	# plot number of selected peaks vs IDR
	# plot all into 1 file
	par(mfcol=c(2,3), pty="s", mar=c(2,5,2,2))
	plot.uri.group(list(uri.all.peaks), NULL, file.name=NULL, 1, title.txt="all peaks")
	plot.uri.group(list(uri.matched.peaks), NULL, file.name=NULL, 1, title.txt="matched peaks")
	plot.ez.group(list(ez), plot.dir=NULL, file.name=NULL, legend.txt=1, y.lim=c(0, 0.6))
	plot(
		rank(-chipTuple$em.output$data.pruned$sample1$signal.value),
		rank(-chipTuple$em.output$data.pruned$sample2$signal.value),
		pch=19, cex=0.01,
		xlab="peak rank rep1", ylab="peak rank rep2", cex.lab=2,
		col=ifelse(idrs <= peaksBelowIDRCutoff, "black", "red")
	)
	title(paste(sum(peaksBelowIDRCutoff), "peaks with IDR <=", idrCutoff))
	legend(0, length(chipTuple$em.output$data.pruned$sample2$signal.value), c(paste("IDR <=", idrCutoff), paste("IDR >", idrCutoff)), col=c("black", "red"), pch=16)

	return(sum(peaksBelowIDRCutoff))
}

#' Calculate irreproducibly discovery rate
#'
#' Assess the consistency between ChIP-Seq replicates according to Qunhua Li et al.
#'
#' Assess the consistency between ChIP-Seq replicates according to Qunhua Li et al.:
#'
#' Ref: Qunhua Li, James B. Brown, Haiyan Huang, and Peter J. Bickel: Measuring reproducibility of high-throughput experiments.
#'      Ann Appl Stat. 2011 October 13 
#'
#' The irreproducible discovery rate (IDR) is a measure for the probability
#' that a peak is called in a ChIP-Seq sample, if it has been called in another
#' sample. Peaks that result from biological activity should be called
#' consistently between replicates. They are assigned a low IDR. In contrast,
#' peaks that are noise are typically not called in all replicates and are
#' assigned a high IDR.
#' The output of this function is meant to be used as input to the function
#' \code{plotIDR}, which generates plots for the assessment of the IDR of
#' ChIP-Seq replicates. Refer to Qunhua Li et al (2011) for an explanation on
#' how to interpret the plots.
#'
#' @export
#' @param chipSamples vector of file names of immuno-precipitated samples in \code{.bam} format or list of \code{GAlignments} objects
#' @param inputSamples vector of file names of control samples in \code{.bam} format or list of \code{GAlignments} objects
#' @param org organism as described in the \code{BSGenome} package: \code{Hsapiens}, \code{Mmusculus}, ...
#' @param assembly assembly name: \code{UCSC}, \code{NCBI}, \code{TAIR}, ...
#' @param version assembly version: \code{hg19}, \code{mm9}, ...
#' @param readLength length of reads (helps identify phantom peak), if \code{NULL} then it is determined automatically from the first 10000 reads in the BAM file
#' @param shiftRange strand shifts at which cross-correlation is evaluated
#' @param binSize step size for \code{shiftRange}
#' @param crossCorrelationPeakShift user-defined cross-correlation peak shift (when given, SPP does not try to detect the peak shift automatically)
#' @param invalidCrossCorrelationPeaks strand shifts to exclude (to avoid phantom peaks)
#' @param halfWidth a numerical value to truncate the peaks to; \code{NULL} means use peak width reported by SPP
#' @param overlapRatio a value between 0 and 1. It controls how much overlaps two peaks need to have to be considered as calling the same region. It is the ratio of overlap / short peak of the two. When set to 0, two peaks are deemed as calling the same region, if they overlap by at least 1 bp
#' @param isBroadPeak if broadpeak is used, set to \code{T}; if narrowpeak is used, set to \code{F}
#' @param cluster a \code{snow} cluster: \code{cluster <- snow::makeCluster(4)}
#' @return object that is suitable as input to \code{plotIDR}
#' @examples
#' chipTuples <- calculateIDR(c("IP1.bam", "IP2.bam"), c("input1.bam", "input2.bam"),
#'                             "Hsapiens", "UCSC", "hg19")
#' for (chipTuple in chipTuples) {
#'     pdf(paste(basename(chipTuple$rep1), "_VS_", basename(chipTuple$rep2), ".pdf", sep=""), paper="a4r", width=11, height=8.5)
#'     plotIDR(chipTuple)
#'     dev.off()
#' }
calculateIDR <- function(chipSamples, inputSamples, org, assembly, version, readLength=NULL, shiftRange=c(-500,1500), binSize=5, crossCorrelationPeakShift=NULL, invalidCrossCorrelationPeaks=c(10, readLength+10), halfWidth=NULL, overlapRatio=0, isBroadPeak=F, cluster=NULL) {

	# read the length of the chromosomes, which will be used to concatenate chr's
	require(paste("BSgenome",org,assembly,version,sep="."),character.only=T)
	chr.size <- data.frame(V1=names(seqlengths(get(org))), V2=seqlengths(get(org)), row.names=NULL)

	# Load SPP library
	require(spp)

	if (is.null(readLength)) {
		# if no read length is given, automatically determine it
		# by taking the median of the first 10000 alignments
		if (class(chipSamples[[1]]) == "GAlignments") {
			require(GenomicAlignments)
			readLength <- round(median(qwidth(head(chipSamples[[1]], 10000))))
		} else if (is.character(chipSamples[[1]])) {
			require(Rsamtools)
			reads <- scanBam(BamFile(chipSamples[[1]], yieldSize=10000), param=ScanBamParam(what="seq"))
			readLength <- round(median(width(reads[[1]]$seq)))
		} else {
			stop("chipSamples must be a list of file paths or GAlignments objects")
		}
		cat(paste("auto-detected read length:", readLength), fill=T)
	}

	# load chip samples
	chip.data <- list()
	for (i in 1:length(chipSamples)) {
		chipSample <- chipSamples[i]
		if (class(chipSample[[1]]) == "GAlignments") {
			fileName <- paste("chipSample", i, sep="")
			alignments <- gAlignmentsToSppTags(chipSample[[1]])
		} else if (is.character(chipSample[[1]])) {
			fileName <- chipSample[[1]]
			alignments <- read.bam.tags(chipSample[[1]])
		} else {
			stop("chipSamples must be a list of file paths or GAlignments objects")
		}
		if (!is.null(names(chipSample)))
			fileName <- names(chipSample)
		chip.data[[length(chip.data)+1]] <- list(fileName=fileName, alignments=alignments)
	}

	# load input samples and merge them into a single pooled sample
	control.data <- list(tags=list(), quality=list())
	for (i in 1:length(inputSamples)) {
		inputSample <- inputSamples[i]
		if (class(inputSample[[1]]) == "GAlignments") {
			fileName <- paste("inputSample", i, sep="")
			bam.tags <- gAlignmentsToSppTags(inputSample[[1]])
		} else if (is.character(inputSample[[1]])) {
			fileName <- inputSample[[1]]
			bam.tags <- read.bam.tags(inputSample[[1]])
		} else {
			stop("inputSamples must be a list of file paths or GAlignments objects")
		}

		# merge samples
		keys <- unique(c(names(bam.tags$tags), names(control.data$tags)))
		control.data$tags <- setNames(mapply(c, bam.tags$tags[keys], control.data$tags[keys], SIMPLIFY=F), keys)
		keys <- unique(c(names(bam.tags$quality), names(control.data$quality)))
		control.data$quality <- setNames(mapply(c, bam.tags$quality[keys], control.data$quality[keys], SIMPLIFY=F), keys)
	}

	tuplesToCompare <- list()
	for (i in 1:(length(chip.data)-1))
		for (j in (i+1):length(chip.data))
			tuplesToCompare[[length(tuplesToCompare)+1]] <- list(rep1=i, rep2=j)

	# create pooled sample
	pooledSample <- list(fileName="pooledChipSamples", alignments=list(tags=list(), quality=list()))
	for (chip.data1 in chip.data) {
		# merge samples
		keys <- unique(c(names(chip.data1$alignments$tags), names(pooledSample$alignments$tags)))
		pooledSample$alignments$tags <- setNames(mapply(c, chip.data1$alignments$tags[keys], pooledSample$alignments$tags[keys], SIMPLIFY=F), keys)
		keys <- unique(c(names(chip.data1$alignments$quality), names(pooledSample$alignments$quality)))
		pooledSample$alignments$quality <- setNames(mapply(c, chip.data1$alignments$quality[keys], pooledSample$alignments$quality[keys], SIMPLIFY=F), keys)
	}
	chip.data[[length(chip.data)+1]] <- pooledSample

	# create pseudo-replicates
	set.seed(001) # ensures reproducibility
	for (chip.data1 in chip.data) {
		cat(paste("generating pseudo-replicates for", chip.data1$fileName), fill=T)
		pseudoRep1 <- list(fileName=paste(chip.data1$fileName, "PseudoRep1", sep="_"), alignments=list(tags=list(), quality=list()))
		pseudoRep2 <- list(fileName=paste(chip.data1$fileName, "PseudoRep2", sep="_"), alignments=list(tags=list(), quality=list()))

		for (i in 1:length(chip.data1$alignments$tags)) {
			# randomly decide which tags go in which pseudo-replicate
			subsample <- sample(c(T, F), length(chip.data1$alignments$tags[[i]]), replace=T)

			contig <- names(chip.data1$alignments$tags[i])

			# pseudo-replicate 1
			pseudoRep1$alignments$tags[[contig]] <- chip.data1$alignments$tags[[i]][subsample]
			pseudoRep1$alignments$quality[[contig]] <- chip.data1$alignments$quality[[i]][subsample]

			# pseudo-replicate 2 (complement of pseudo-replicate 1)
			pseudoRep2$alignments$tags[[contig]] <- chip.data1$alignments$tags[[i]][!subsample]
			pseudoRep2$alignments$quality[[contig]] <- chip.data1$alignments$quality[[i]][!subsample]
		}

		chip.data[[length(chip.data)+1]] <- pseudoRep1
		chip.data[[length(chip.data)+1]] <- pseudoRep2
		tuplesToCompare[[length(tuplesToCompare)+1]] <- list(rep1=length(chip.data)-1, rep2=length(chip.data))
	}

	# call peaks
	for (i in unique(as.integer(unlist(tuplesToCompare)))) {
		cat(paste("finding peaks for", chip.data[[i]]$fileName), fill=T)
		chip.data[[i]]$peaks <- find.peaks(chip.data[[i]]$alignments, control.data, readLength, shiftRange, binSize, crossCorrelationPeakShift, invalidCrossCorrelationPeaks, cluster)
		chip.data[[i]]$alignments <- list() # free memory
	}
	control.data <- list() # free memory

	############# process the data
	# process data, summit: the representation of the location of summit
	computeURIandEM <- function(tupleToCompare, chr.size, halfWidth, isBroadPeak, overlapRatio) {
		chip.tuple <- list(rep1=tupleToCompare$rep1$fileName, rep2=tupleToCompare$rep2$fileName)

		rep1 <- process.narrowpeak(tupleToCompare$rep1$peaks, chr.size, half.width=halfWidth, summit="offset", broadpeak=isBroadPeak)
		rep2 <- process.narrowpeak(tupleToCompare$rep2$peaks, chr.size, half.width=halfWidth, summit="offset", broadpeak=isBroadPeak)

		cat(paste("computing correspondence profile (URI) for", tupleToCompare$rep1$fileName, "and", tupleToCompare$rep2$fileName), fill=T)
		uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned, sig.value1="signal.value", sig.value2="signal.value", overlap.ratio=overlapRatio)
		chip.tuple$uri.output <- uri.output

		cat(paste("computing EM procedure for", tupleToCompare$rep1$fileName, "and", tupleToCompare$rep2$fileName), fill=T)
		em.output <- fit.em(uri.output$data12.enrich, fix.rho2=T)
		chip.tuple$em.output <- em.output

		return(chip.tuple)
	}
	tuplesToCompare <- lapply(tuplesToCompare, function(tupleToCompare, chip.data) { return(list(rep1=chip.data[[tupleToCompare$rep1]], rep2=chip.data[[tupleToCompare$rep2]])) }, chip.data)
	if (is.null(cluster)) {
		chip.tuples <- list()
		for (tupleToCompare in tuplesToCompare) {
			chip.tuples[[length(chip.tuples)+1]] <- computeURIandEM(tupleToCompare, chr.size, halfWidth, isBroadPeak, overlapRatio)
		}
	} else {
		clusterExport(cluster, as.vector(lsf.str(envir=.GlobalEnv))) # export all functions to cluster nodes
		chip.tuples <- clusterApplyLB(cluster, tuplesToCompare, computeURIandEM, chr.size, halfWidth, isBroadPeak, overlapRatio)
	}
	return(chip.tuples)
}

#' Get self-consistency ratios
#'
#' Calculate the ratio of peaks below a given IDR threshold between all pseudo-replicates
#'
#' The ChIP-Seq guidelines of the ENCODE consortium recommend that the number
#' of peaks with an IDR below 0.01 should not differ by more than a factor of 2
#' between any pair of pseudo-replicates. This function counts the peaks of all
#' pseudo-replicates with an IDR <= 0.01 and calculates the quotient of
#' the number of peaks for all possible pairs of pseudo-replicates. It produces
#' a N x N matrix, where N is the number of true replicates. The values of the
#' matrix are calculated as:
#'
#'   [peak count of replicate in column] / [peak count of replicate in row]
#'
#' All values should be within a range of 0.5 to 2. Deviations from this
#' range indicate low reproducibility.
#'
#' Ref: Stephen G. Landt, Georgi K. Marinov, Anshul Kundaje, Pouya Kheradpour,
#' Florencia Pauli, Serafim Batzoglou, Bradley E. Bernstein, Peter Bickel,
#' James B. Brown, Philip Cayting, Yiwen Chen, Gilberto DeSalvo,
#' Charles Epstein, Katherine I. Fisher-Aylor, Ghia Euskirchen, Mark Gerstein,
#' Jason Gertz, Alexander J. Hartemink, Michael M. Hoffman,
#' Vishwanath R. Iyer, Youngsook L. Jung, Subhradip Karmakar, Manolis Kellis,
#' Peter V. Kharchenko, Qunhua Li, Tao Liu, X. Shirley Liu, Lijia Ma,
#' Aleksandar Milosavljevic, Richard M. Myers, Peter J. Park,
#' Michael J. Pazin, Marc D. Perry, Debasish Raha, Timothy E. Reddy,
#' Joel Rozowsky, Noam Shoresh, Arend Sidow, Matthew Slattery,
#' John A. Stamatoyannopoulos, Michael Y. Tolstorukov, Kevin P. White,
#' Simon Xi, Peggy J. Farnham, Jason D. Lieb, Barbara J. Wold, Michael Snyder:
#' ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia
#' Genome Res. Sep 2012; 22(9): 1813–1831.
#'
#' @export
#' @param chipTuples output of \code{calculateIDR}
#' @param idrCutoff count only peaks below this IDR value (for mammalian cells a value of 0.01 is recommended)
#' @return matrix with quotients of peak counts for all possible combinations of replicates
#' @examples
#' chipTuples <- calculateIDR(c("IP1.bam", "IP2.bam"), c("input1.bam", "input2.bam"),
#'                             "Hsapiens", "UCSC", "hg19")
#' getSelfConsistencyRatios(chipTuples)
getSelfConsistencyRatios <- function(chipTuples, idrCutoff=0.01) {
	replicateNames <- c()
	peaksBelowIDRCutoff <- c()
	for (chip.tuple in chipTuples)
		if (length(grep("_PseudoRep.$", chip.tuple$rep1)) > 0 & length(grep("pooledChipSamples_PseudoRep.$", chip.tuple$rep1)) == 0) {
			replicateNames <- c(replicateNames, paste(chip.tuple$rep1, chip.tuple$rep2, sep="_VS_"))
			peaksBelowIDRCutoff <- c(peaksBelowIDRCutoff, sum(1-chip.tuple$em.output$em.fit$e.z <= idrCutoff))
		}
	result = sapply(peaksBelowIDRCutoff, function(x) { x/peaksBelowIDRCutoff })
	colnames(result) <- replicateNames
	rownames(result) <- replicateNames
	return(result)
}

#' Get cross-consistency ratios
#'
#' Calculate the ratio of peaks below a given IDR threshold between all true replicates and the pooled pseudo-replicates
#'
#' The ChIP-Seq guidelines of the ENCODE consortium recommend that the number
#' of peaks with an IDR below 0.01 should not differ by more than a factor of 2
#' between the pair of pooled pseudo-replicates and any pair of true replicates.
#' This function counts the peaks of all pairs of true replicates with an IDR
#' below 0.01 as well as the peaks of pseudo-replicates from a pooled sample.
#' The output is a list of quotients calculated as:
#'
#'   [peak count of pooled replicates] / [peak count of true replicate]
#'
#' All values should be lower than or equal to 2. Higher values indicate low
#' reproducibility.
#'
#' Ref: Stephen G. Landt, Georgi K. Marinov, Anshul Kundaje, Pouya Kheradpour,
#' Florencia Pauli, Serafim Batzoglou, Bradley E. Bernstein, Peter Bickel,
#' James B. Brown, Philip Cayting, Yiwen Chen, Gilberto DeSalvo,
#' Charles Epstein, Katherine I. Fisher-Aylor, Ghia Euskirchen, Mark Gerstein,
#' Jason Gertz, Alexander J. Hartemink, Michael M. Hoffman,
#' Vishwanath R. Iyer, Youngsook L. Jung, Subhradip Karmakar, Manolis Kellis,
#' Peter V. Kharchenko, Qunhua Li, Tao Liu, X. Shirley Liu, Lijia Ma,
#' Aleksandar Milosavljevic, Richard M. Myers, Peter J. Park,
#' Michael J. Pazin, Marc D. Perry, Debasish Raha, Timothy E. Reddy,
#' Joel Rozowsky, Noam Shoresh, Arend Sidow, Matthew Slattery,
#' John A. Stamatoyannopoulos, Michael Y. Tolstorukov, Kevin P. White,
#' Simon Xi, Peggy J. Farnham, Jason D. Lieb, Barbara J. Wold, Michael Snyder:
#' ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia
#' Genome Res. Sep 2012; 22(9): 1813–1831.
#'
#' @export
#' @param chipTuples output of \code{calculateIDR}
#' @param idrCutoff count only peaks below this IDR value (for mammalian cells a value of 0.01 is recommended)
#' @return list of quotients of peak counts (pooled pseudo-replicates / true replicates)
#' @examples
#' chipTuples <- calculateIDR(c("IP1.bam", "IP2.bam"), c("input1.bam", "input2.bam"),
#'                             "Hsapiens", "UCSC", "hg19")
#' getCrossConsistencyRatios(chipTuples)
getCrossConsistencyRatios <- function(chipTuples, idrCutoff=0.01) {
	replicateNames <- c()
	peaksBelowIDRCutoff <- c()
	for (chip.tuple in chipTuples)
		if (length(grep("_PseudoRep.$", chip.tuple$rep1)) == 0) {
			replicateNames <- c(replicateNames, paste(chip.tuple$rep1, chip.tuple$rep2, sep="_VS_"))
			peaksBelowIDRCutoff <- c(peaksBelowIDRCutoff, sum(1-chip.tuple$em.output$em.fit$e.z <= idrCutoff))
		} else if (chip.tuple$rep1 == "pooledChipSamples_PseudoRep1" & chip.tuple$rep2 == "pooledChipSamples_PseudoRep2") {
			pooledPeaksBelowIDRCutoff <- sum(1-chip.tuple$em.output$em.fit$e.z <= idrCutoff)
		}
	result = sapply(peaksBelowIDRCutoff, function(x) { pooledPeaksBelowIDRCutoff/x })
	names(result) <- replicateNames
	return(result)
}

##################### internal functions #####################

# use peak-caller SPP to find peaks
find.peaks <- function(chipData, pooledInputData, readLength, shiftRange=c(-500,1500), binSize=5, crossCorrelationPeakShift=NULL, invalidCrossCorrelationPeaks=NULL, cluster=NULL) {

	# Read ChIP tagAlign/BAM files
	chipData$num.tags <- sum(unlist(lapply(chipData$tags,function(d) length(d))))
	pooledInputData$num.tags <- sum(unlist(lapply(pooledInputData$tags,function(d) length(d))))

	# #################################    
	# Calculate cross-correlation for various strand shifts
	# #################################    
	# crosscorr
	# $cross.correlation : Cross-correlation profile as an $x/$y data.frame
	# $peak : Position ($x) and height ($y) of automatically detected cross-correlation peak.
	# $whs: Optimized window half-size for binding detection (based on the width of the cross-correlation peak) 
	crosscorr <- get.binding.characteristics(chipData,
	                remove.tag.anomalies=F,
			srange=shiftRange,
			bin=binSize,
			accept.all.tags=T,
			cluster=cluster)
	
	# Smooth the cross-correlation curve if required
	cc <- crosscorr$cross.correlation
	crosscorr$min.cc <- crosscorr$cross.correlation[ length(crosscorr$cross.correlation$y) , ] # minimum value and shift of cross-correlation
	sbw <- 2*floor(ceiling(5/binSize) / 2) + 1 # smoothing bandwidth
	cc$y <- runmean(cc$y,sbw,alg="fast")
	
	# Compute cross-correlation peak
	bw <- ceiling(2/binSize) # crosscorr[i] is compared to crosscorr[i+/-bw] to find peaks
	peakidx <- (diff(cc$y,bw)>=0) # cc[i] > cc[i-bw]
	peakidx <- diff(peakidx,bw)
	peakidx <- which(peakidx==-1) + bw        
	
	# exclude peaks within invalid region
	if ( is.null(invalidCrossCorrelationPeaks) ) {
		invalidCrossCorrelationPeaks <- c(10, readLength+10)
	}
	peakidx <- peakidx[(cc$x[peakidx] < invalidCrossCorrelationPeaks[1]) | (cc$x[peakidx] > invalidCrossCorrelationPeaks[2]) | (cc$x[peakidx] < 0) ]    
	cc <- cc[peakidx,]
	
	# Find max peak position and other peaks within 0.9*max_peakvalue that are further away from maxpeakposition   
	maxpeakidx <- which.max(cc$y)
	maxpeakshift <- cc$x[maxpeakidx]
	maxpeakval <- cc$y[maxpeakidx]
	peakidx <-which((cc$y >= 0.9*maxpeakval) & (cc$x >= maxpeakshift)) 
	cc <- cc[peakidx,]
	
	# pick best peak
	cc.peak <- cc[order(cc$y,decreasing=TRUE),]
	
	# Override peak shift if user supplies peak shift
	if (! is.null(crossCorrelationPeakShift)) {
		cc.peak <- approx(crosscorr$cross.correlation$x,crosscorr$cross.correlation$y,crossCorrelationPeakShift,rule=2)
	}
	
	# Reset values in crosscorr
	crosscorr$peak$x <- cc.peak$x[1]
	crosscorr$peak$y <- cc.peak$y[1]
	
	# Compute window half size
	whs.thresh <- crosscorr$min.cc$y + (crosscorr$peak$y - crosscorr$min.cc$y)/3
	crosscorr$whs <- max(crosscorr$cross.correlation$x[crosscorr$cross.correlation$y >= whs.thresh])
	
	# #################################    
	# Call peaks
	# #################################
	
	# Remove local tag anomalies
	chipData <- remove.local.tag.anomalies(chipData$tags)
	pooledInputData <- remove.local.tag.anomalies(pooledInputData$tags)
	
	# Find peaks
	narrow.peaks <- find.binding.positions(signal.data=chipData,control.data=pooledInputData,fdr=0.99,method=tag.lwcc,whs=crosscorr$whs,cluster=cluster)
	region.peaks <- add.broad.peak.regions(chipData,pooledInputData,narrow.peaks,window.size=max(50,round(crosscorr$whs/4)),z.thr=9)

	# convert peaks to data frame (as in write.narrowpeak.binding)
	margin <- round(crosscorr$whs/2)
	chrl <- names(region.peaks$npl)
	names(chrl) <- chrl
	md <- do.call(rbind, lapply(chrl, function(chr) {
		df <- region.peaks$npl[[chr]]
		x <- df$x
		rs <- df$rs
		if (is.null(rs)) {
			rs <- rep(NA, length(x))
		}
		re <- df$re
		if (is.null(re)) {
			re <- rep(NA, length(x))
		}
		ivi <- which(is.na(rs))
		if (any(ivi)) {
			rs[ivi] <- x[ivi] - margin
		}
		ivi <- which(is.na(re))
		if (any(ivi)) {
			re[ivi] <- x[ivi] + margin
		}
		cbind(chr, rs, re, df$y, -1, df$fdr, x - rs)
	}))
	md <- md[order(as.numeric(md[, 4]), decreasing = T), ]
        md <- md[1:min(nrow(md),300000),]
	return(md)
}

# revised on 2-20-10
# - fix error in pass.structure: reverse rank.combined, so that big sig.value
#  are ranked with small numbers (1, 2, ...)
# - fix error on get.ez.tt.all: get ez.cutoff from sorted e.z

#
# modified EM procedure to compute empirical CDF more precisely - 09/2009



# this file contains the functions for  
# 1. computing the correspondence profile (upper rank intersection and derivatives)
# 2. inference of copula mixture model
#
# It also has functions for
# 1. reading peak caller results
# 2. processing and matching called peaks
# 3. plotting results


################ read peak caller results

# process narrow peak format
# some peak callers may not report q-values, p-values or fold of enrichment
# need further process before comparison
#
# stop.exclusive: Is the basepair of peak.list$stop exclusive? In narrowpeak and broadpeak format they are exclusive.
# If it is exclusive, we need subtract peak.list$stop by 1 to avoid the same basepair being both a start and a stop of two 
# adjacent peaks, which creates trouble for finding correct intersect  
process.narrowpeak <- function(narrow.data, chr.size, half.width=NULL, summit="offset", stop.exclusive=T, broadpeak=F){

  if(broadpeak){
    bb.ori <- data.frame(chr=narrow.data[,1], start=as.numeric(narrow.data[,2]), stop=as.numeric(narrow.data[,3]), signal.value=as.numeric(narrow.data[,4]), p.value=as.numeric(narrow.data[,5]), q.value=as.numeric(narrow.data[,6]))
  }else{
    bb.ori <- data.frame(chr=narrow.data[,1], start=as.numeric(narrow.data[,2]), stop=as.numeric(narrow.data[,3]), signal.value=as.numeric(narrow.data[,4]), p.value=as.numeric(narrow.data[,5]), q.value=as.numeric(narrow.data[,6]), summit=as.numeric(narrow.data[,7]))
  }

  if(summit=="summit"){
    bb.ori$summit <- bb.ori$summit-bb.ori$start # change summit to offset to avoid error when concatenating chromosomes
  }
 
  bb <- concatenate.chr(bb.ori, chr.size)

  #bb <- bb.ori

  # remove the peaks that has the same start and stop value
  bb <- bb[bb$start != bb$stop,]

  if(stop.exclusive==T){
    bb$stop <- bb$stop-1
  }

  if(!is.null(half.width)){
    bb$start.ori <- bb$start    #Anshul changed this
    bb$stop.ori <- bb$stop 	#Anshul changed this

    # if peak is narrower than the specified window, stay with its width
    # otherwise chop wider peaks to specified width
    width <- bb$stop-bb$start +1
    is.wider <- width > 2*half.width

    if(summit=="offset" | summit=="summit"){ # if summit is offset from start
      bb$start[is.wider] <- bb$start.ori[is.wider] + bb$summit[is.wider]-half.width
      bb$stop[is.wider] <- bb$start.ori[is.wider] + bb$summit[is.wider]+half.width
    } else { 
      if(summit=="unknown"){
        bb$start[is.wider] <- bb$start.ori[is.wider]+round(width[is.wider]/2) - half.width
        bb$stop[is.wider] <- bb$start.ori[is.wider]+round(width[is.wider]/2) + half.width
      }
    }

    bb$start.ori <- bb.ori$start    #Anshul changed this
    bb$stop.ori <- bb.ori$stop      #Anshul changed this
  }

  bb <- clean.data(bb)
  invisible(list(data.ori=bb.ori, data.cleaned=bb))
}

# clean data 
# and concatenate chromosomes if needed
clean.data <- function(adata){

  # remove the peaks that has the same start and stop value
  adata <- adata[adata$start != adata$stop,]

  # if some stops and starts are the same, need fix them
  stop.in.start <- is.element(adata$stop, adata$start)
  n.fix <- sum(stop.in.start)
  if(n.fix >0){
    print(paste("Fix", n.fix, "stops\n"))
    adata$stop[stop.in.start] <- adata$stop[stop.in.start]-1 
  }  
 
  return(adata) 
}

# concatenate peaks
# peaks: the dataframe to have all the peaks
# chr.file: the file to keep the length of each chromosome 
# chr files should come from the species that the data is from
concatenate.chr <- function(peaks, chr.size){

 # chr.size <- read.table(chr.file)
  chr.o <- order(chr.size[,1])
  chr.size <- chr.size[chr.o,]

  chr.shift <- cumsum(c(0, chr.size[-nrow(chr.size),2]))
  chr.size.cum <- data.frame(chr=chr.size[,1], shift=chr.shift)  

  peaks$start.ori <- peaks$start
  peaks$stop.ori <- peaks$stop
  
  for(i in 1:nrow(chr.size)){
    is.in <- as.character(peaks$chr) == as.character(chr.size.cum$chr[i])
    if(sum(is.in)>0){
      peaks[is.in,]$start <- peaks[is.in,]$start + chr.size.cum$shift[i]
      peaks[is.in,]$stop <- peaks[is.in,]$stop + chr.size.cum$shift[i]
    }
  }

  invisible(peaks)
}


deconcatenate.chr <- function(peaks, chr.size){

  chr.o <- order(chr.size[,1])
  chr.size <- chr.size[chr.o,]

  chr.shift <- cumsum(c(0, chr.size[-nrow(chr.size),2]))
  chr.size.cum <- data.frame(chr=chr.size[,1], shift=chr.shift)  

  peaks$chr <- rep(NA, nrow(peaks))
  
  for(i in 1:(nrow(chr.size.cum)-1)){
    is.in <- peaks$start > chr.size.cum[i,2] & peaks$start <= chr.size.cum[i+1, 2]
    if(sum(is.in)>0){
      peaks[is.in,]$start <- peaks[is.in,]$start - chr.size.cum[i,2]
      peaks[is.in,]$stop <- peaks[is.in,]$stop - chr.size.cum[i,2]+1    
      peaks[is.in,]$chr <- chr.size[i,1]
    }
  }

  if(i == nrow(chr.size.cum)){
    is.in <- peaks$start > chr.size.cum[i, 2]
    if(sum(is.in)>0){
      peaks[is.in,]$start <- peaks[is.in,]$start - chr.size.cum[i,2]
      peaks[is.in,]$stop <- peaks[is.in,]$stop - chr.size.cum[i,2]+1    
      peaks[is.in,]$chr <- chr.size[i,1]
    }
  }
  
  invisible(peaks)
}

################ preprocessing peak calling output


# 
# read two calling results and sort by peak starting locations, 
# then find overlap between peaks
# INPUT:
#   rep1: the 1st replicate
#   rep2: the 2nd replicate
# OUTPUT:
#   id1, id2: the labels for the identified peaks on the replicates
find.overlap <- function(rep1, rep2){

  o1 <- order(rep1$start)
  rep1 <- rep1[o1,]
    
  o2 <- order(rep2$start)
  rep2 <- rep2[o2,]

  n1 <- length(o1)
  n2 <- length(o2)
  
  # assign common ID to peaks
  id1 <- rep(0, n1) # ID assigned on rep1
  id2 <- rep(0, n2) # ID assigned on rep2
  id <- 1 # keep track common id's
  
  # check if two replicates overlap with each other
  i <- 1
  j <- 1

  while(i <= n1|| j <= n2){

    # && (id1[n1] ==0 || id2[n2] ==0)
    
    # if one list runs out
    if(i > n1 && j < n2){
      
      j <- j+1
      id2[j] <- id
      id <- id +1
      next
    } else{
      if(j > n2 && i < n1){
        i <- i+1        
        id1[i] <- id
        id <- id +1
        next
      } else {
        if(i >= n1 && j >=n2)
          break
      }
    }

    # if not overlap

    if(!(rep1$start[i] <= rep2$stop[j] && rep2$start[j] <= rep1$stop[i])){

      # at the start of loop, when both are not assigned an ID
      # the one locates in front is assigned first
      if(id1[i] ==0 && id2[j]==0){
        if(rep1$stop[i] < rep2$stop[j]){
          id1[i] <- id
        } else {
          id2[j] <- id
        }
      } else { # in the middle of the loop, when one is already assigned
      # The one that has not assigned gets assigned
      #  if(id1[i] ==0){ # id1[i] is not assigned
      #    id1[i] <- id
      #  } else { # id2[i] is not assigned
      #    id2[j] <- id 
      #  }

        # order the id according to location
        if(rep1$stop[i] <= rep2$stop[j]){
          id1[i] <- max(id2[j], id1[i])
          id2[j] <- id  
        } else {
          if(rep1$stop[i] > rep2$stop[j]){
            id2[j] <- max(id1[i], id2[j])
            id1[i] <- id
          }
        }
          
      }
      
      id <- id +1
      
    } else { # if overlap
    
      if(id1[i] == 0 && id2[j] == 0){ # not assign label yet
        id1[i] <- id 
        id2[j] <- id
        id <- id +1
      } else { # one peak is already assigned label, the other is 0
        
        id1[i] <- max(id1[i], id2[j]) # this is a way to copy the label of the assigned peak without knowing which one is already assigned
        id2[j] <- id1[i] # syncronize the labels        
      }
      
    }
    
    if(rep1$stop[i] < rep2$stop[j]){
      i <- i+1
    } else {
      j <- j+1
    } 
    
  }

  invisible(list(id1=id1, id2=id2))
  
}

# Impute the missing significant value for the peaks called only on one replicate.
# value 
# INPUT:
#   rep1, rep2: the two peak calling output 
#   id1, id2: the IDs assigned by function find.overlap, vectors
#        If id1[i]==id2[j], peak i on rep1 overlaps with peak j on rep2
#   p.value.impute: the significant value to impute for the missing peaks 
# OUTPUT:   
#   rep1, rep2: peaks ordered by the start locations with imputed peaks
#   id1, id2: the IDs with imputed peaks
fill.missing.peaks <- function(rep1, rep2, id1, id2, p.value.impute){

#   rep1 <- data.frame(chr=rep1$chr, start=rep1$start, stop=rep1$stop, sig.value=rep1$sig.value)
#   rep2 <- data.frame(chr=rep2$chr, start=rep2$start, stop=rep2$stop, sig.value=rep2$sig.value)   
   
   o1 <- order(rep1$start)
   rep1 <- rep1[o1,]
    
   o2 <- order(rep2$start)
   rep2 <- rep2[o2,]  
     
   entry.in1.not2 <- !is.element(id1, id2)
   entry.in2.not1 <- !is.element(id2, id1)

   if(sum(entry.in1.not2) > 0){
   
     temp1 <- rep1[entry.in1.not2, ]

     # impute sig.value
     temp1$sig.value <- p.value.impute
     temp1$signal.value <- p.value.impute
     temp1$p.value <- p.value.impute
     temp1$q.value <- p.value.impute
     
     rep2.filled <- rbind(rep2, temp1)
     id2.filled <- c(id2, id1[entry.in1.not2])
   } else {
     id2.filled <- id2
     rep2.filled <- rep2
   }

   if(sum(entry.in2.not1) > 0){

     temp2 <- rep2[entry.in2.not1, ]

     # fill in p.values to 1
     temp2$sig.value <- p.value.impute
     temp2$signal.value <- p.value.impute
     temp2$p.value <- p.value.impute
     temp2$q.value <- p.value.impute
   

     # append to the end
     rep1.filled <- rbind(rep1, temp2)

     id1.filled <- c(id1, id2[entry.in2.not1])
   } else {
     id1.filled <- id1
     rep1.filled <- rep1
   }

   # sort rep1 and rep2 by the same id
   o1 <- order(id1.filled)
   rep1.ordered <- rep1.filled[o1, ]

   o2 <- order(id2.filled)
   rep2.ordered <- rep2.filled[o2, ]   
   
   invisible(list(rep1=rep1.ordered, rep2=rep2.ordered,
                  id1=id1.filled[o1], id2=id2.filled[o2]))
 }

# Merge peaks with same ID on the same replicates 
# (They are generated if two peaks on rep1 map to the same peak on rep2)
# need peak.list have 3 columns: start, stop and sig.value 
merge.peaks.best <- function(peak.list, id){

  i <- 1
  j <- 1
  dup.index <- c()
  sig.value <- c()
  start.new <- c()
  stop.new <- c()
  id.new <- c()

  # original data
  chr <- c()
  start.ori <- c()
  stop.ori <- c()
  
  signal.value <- c()
  p.value <- c()
  q.value <- c()

  while(i < length(id)){
    
    if(id[i] == id[i+1]){
      dup.index <- c(dup.index, i, i+1) # push on dup.index
    } else {
      if(length(dup.index)>0){ # pop from dup.index        
    #    sig.value[j] <- mean(peak.list$sig.value[unique(dup.index)]) # mean of -log(pvalue)
        sig.value[j] <- max(peak.list$sig.value[unique(dup.index)])
        start.new[j] <- peak.list$start[min(dup.index)]
        stop.new[j] <- peak.list$stop[max(dup.index)]
        id.new[j] <- id[max(dup.index)]
        
    #    signal.value[j] <- mean(peak.list$signal.value[unique(dup.index)])        #    p.value[j] <- mean(peak.list$p.value[unique(dup.index)]) # mean of -log(pvalue)
    #    q.value[j] <- mean(peak.list$q.value[unique(dup.index)]) # mean of -log(pvalue)     
        signal.value[j] <- max(peak.list$signal.value[unique(dup.index)]) 
        p.value[j] <- max(peak.list$p.value[unique(dup.index)]) 
        q.value[j] <- max(peak.list$q.value[unique(dup.index)]) 

        chr[j] <- as.character(peak.list$chr[min(dup.index)])
        start.ori[j] <- peak.list$start.ori[min(dup.index)]
        stop.ori[j] <- peak.list$stop.ori[max(dup.index)]
        
        dup.index <- c()
      } else { # nothing to pop
        sig.value[j] <- peak.list$sig.value[i]
        start.new[j] <- peak.list$start[i]
        stop.new[j] <- peak.list$stop[i]
        id.new[j] <- id[i]

        signal.value[j] <- peak.list$signal.value[i] 
        p.value[j] <- peak.list$p.value[i] 
        q.value[j] <- peak.list$q.value[i] 

        chr[j] <- as.character(peak.list$chr[i])
        start.ori[j] <- peak.list$start.ori[i]
        stop.ori[j] <- peak.list$stop.ori[i]
        
      }
      j <- j+1
    }
    i <- i+1
  }

  data.new <- data.frame(id=id.new, sig.value=sig.value, start=start.new, stop=stop.new, signal.value=signal.value, p.value=p.value, q.value=q.value, chr=chr, start.ori=start.ori, stop.ori=stop.ori)
  invisible(data.new)
}

# Merge peaks with same ID on the same replicates 
# (They are generated if two peaks on rep1 map to the same peak on rep2)
# need peak.list have 3 columns: start, stop and sig.value 
merge.peaks <- function(peak.list, id){

  i <- 1
  j <- 1
  dup.index <- c()
  sig.value <- c()
  start.new <- c()
  stop.new <- c()
  id.new <- c()

  # original data
  chr <- c()
  start.ori <- c()
  stop.ori <- c()
  
  signal.value <- c()
  p.value <- c()
  q.value <- c()

  while(i < length(id)){
    
    if(id[i] == id[i+1]){
      dup.index <- c(dup.index, i, i+1) # push on dup.index
    } else {
      if(length(dup.index)>0){ # pop from dup.index
        sig.value[j] <- mean(peak.list$sig.value[unique(dup.index)]) # mean of -log(pvalue)
        start.new[j] <- peak.list$start[min(dup.index)]
        stop.new[j] <- peak.list$stop[max(dup.index)]
        id.new[j] <- id[max(dup.index)]
        
        signal.value[j] <- mean(peak.list$signal.value[unique(dup.index)]) # mean of -log(pvalue)
        p.value[j] <- mean(peak.list$p.value[unique(dup.index)]) # mean of -log(pvalue)
        q.value[j] <- mean(peak.list$q.value[unique(dup.index)]) # mean of -log(pvalue)

        chr[j] <- as.character(peak.list$chr[min(dup.index)])
        start.ori[j] <- peak.list$start.ori[min(dup.index)]
        stop.ori[j] <- peak.list$stop.ori[max(dup.index)]
        
        dup.index <- c()
      } else { # nothing to pop
        sig.value[j] <- peak.list$sig.value[i]
        start.new[j] <- peak.list$start[i]
        stop.new[j] <- peak.list$stop[i]
        id.new[j] <- id[i]

        signal.value[j] <- peak.list$signal.value[i] 
        p.value[j] <- peak.list$p.value[i] 
        q.value[j] <- peak.list$q.value[i] 

        chr[j] <- as.character(peak.list$chr[i])
        start.ori[j] <- peak.list$start.ori[i]
        stop.ori[j] <- peak.list$stop.ori[i]
        
      }
      j <- j+1
    }
    i <- i+1
  }

  data.new <- data.frame(id=id.new, sig.value=sig.value, start=start.new, stop=stop.new, signal.value=signal.value, p.value=p.value, q.value=q.value, chr=chr, start.ori=start.ori, stop.ori=stop.ori)
  invisible(data.new)
}





# a wrap function to fill in missing peaks, merge peaks and impute significant values
# out1 and out2 are two peak calling outputs
pair.peaks <- function(out1, out2, p.value.impute=0){

  aa <- find.overlap(out1, out2)
  bb <- fill.missing.peaks(out1, out2, aa$id1, aa$id2, p.value.impute=0)

  cc1 <- merge.peaks(bb$rep1, bb$id1)
  cc2 <- merge.peaks(bb$rep2, bb$id2)

  invisible(list(merge1=cc1, merge2=cc2))
}



# overlap.ratio is a parameter to define the percentage of overlap
# if overlap.ratio =0, 1 basepair overlap is counted as overlap
# if overlap.ratio between 0 and 1, it is the minimum proportion of
# overlap required to be called as a match
# it is computed as the overlap part/min(peak1.length, peak2.length)
pair.peaks.filter <- function(out1, out2, p.value.impute=0, overlap.ratio=0){

  aa <- find.overlap(out1, out2)
  bb <- fill.missing.peaks(out1, out2, aa$id1, aa$id2, p.value.impute=0)

  cc1 <- merge.peaks(bb$rep1, bb$id1)
  cc2 <- merge.peaks(bb$rep2, bb$id2)

  frag12 <- cbind(cc1$start, cc1$stop, cc2$start, cc2$stop)
  
  frag.ratio <- apply(frag12, 1, overlap.middle)

  frag.ratio[cc1$sig.value==p.value.impute | cc2$sig.value==p.value.impute] <- 0

  cc1$frag.ratio <- frag.ratio
  cc2$frag.ratio <- frag.ratio

  merge1 <- cc1[cc1$frag.ratio >= overlap.ratio,]
  merge2 <- cc2[cc2$frag.ratio >= overlap.ratio,]
  
  invisible(list(merge1=merge1, merge2=merge2))
}

# x[1], x[2] are the start and end of the first fragment
# and x[3] and x[4] are the start and end of the 2nd fragment 
# If there are two fragments, we can find the overlap by ordering the
# start and stop of all the ends and find the difference between the middle two
overlap.middle  <- function(x){

  x.o <- x[order(x)]
  f1 <- x[2]-x[1]
  f2 <- x[4]-x[3]
  
  f.overlap <- abs(x.o[3]-x.o[2])
  f.overlap.ratio <- f.overlap/min(f1, f2)

  return(f.overlap.ratio)
}



#######
####### compute correspondence profile
#######

# compute upper rank intersection for one t
# tv: the upper percentile
# x is sorted by the order of paired variable
comp.uri <- function(tv, x){
  n <- length(x)
  qt <- quantile(x, prob=1-tv[1]) # tv[1] is t
#  sum(x[1:ceiling(n*tv[2])] >= qt)/n/tv[2]- tv[1]*tv[2] #tv[2] is v
  sum(x[1:ceiling(n*tv[2])] >= qt)/n

}

# compute the correspondence profile
# tt, vv: vector between (0, 1) for percentages
get.uri.2d <- function(x1, x2, tt, vv, spline.df=NULL){

  o <- order(x1, x2, decreasing=T)
  
  # sort x2 by the order of x1
  x2.ordered <- x2[o]
  
  tv <- cbind(tt, vv)
  ntotal <- length(x1) # number of peaks    

  uri <- apply(tv, 1, comp.uri, x=x2.ordered)

  # compute the derivative of URI vs t using small bins
  uri.binned <- uri[seq(1, length(uri), by=4)]
  tt.binned <- tt[seq(1, length(uri), by=4)]
  uri.slope <- (uri.binned[2:(length(uri.binned))] - uri.binned[1:(length(uri.binned)-1)])/(tt.binned[2:(length(uri.binned))] - tt.binned[1:(length(tt.binned)-1)])

  # smooth uri using spline
  # first find where the jump is and don't fit the jump
  # this is the index on the left
  # jump.left.old  <- which.max(uri[-1]-uri[-length(uri)])
  short.list.length <- min(sum(x1>0)/length(x1), sum(x2>0)/length(x2))

  if(short.list.length < max(tt)){
    jump.left <- which(tt>short.list.length)[1]-1
  } else {
    jump.left <- which.max(tt)
  }

#  reversed.index <- seq(length(tt), 1, by=-1)
#  nequal <- sum(uri[reversed.index]== tt[reversed.index])
#  temp  <- which(uri[reversed.index]== tt[reversed.index])[nequal]
#  jump.left <- length(tt)-temp
 
  if(jump.left < 6){
   jump.left <- length(tt)
  }
    
 
  if(is.null(spline.df))
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=6.4)
  else{
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=spline.df)
  }
  # predict the first derivative
  uri.der <- predict(uri.spl, tt[1:jump.left], deriv=1)

  invisible(list(tv=tv, uri=uri, 
                 uri.slope=uri.slope, t.binned=tt.binned[2:length(uri.binned)], 
                 uri.spl=uri.spl, uri.der=uri.der, jump.left=jump.left,
                 ntotal=ntotal))
 }


# change the scale of uri from based on t (percentage) to n (number of peaks or basepairs)
# this is for plotting multiple pairwise URI's on the same plot 
scale.t2n <- function(uri){

  ntotal <- uri$ntotal
  tv <- uri$tv*uri$ntotal
  uri.uri <- uri$uri*uri$ntotal
  jump.left <- uri$jump.left
  uri.spl <- uri$uri.spl
  uri.spl$x <- uri$uri.spl$x*uri$ntotal 
  uri.spl$y <- uri$uri.spl$y*uri$ntotal

  t.binned <- uri$t.binned*uri$ntotal
  uri.slope <- uri$uri.slope
  uri.der <- uri$uri.der
  uri.der$x <- uri$uri.der$x*uri$ntotal
  uri.der$y <- uri$uri.der$y

  uri.n <- list(tv=tv, uri=uri.uri, t.binned=t.binned, uri.slope=uri.slope, uri.spl=uri.spl, uri.der=uri.der, ntotal=ntotal, jump.left=jump.left)
  return(uri.n)
} 




# a wrapper for running URI for peaks from peak calling results
# both data1 and data2 are calling results in narrowpeak format
compute.pair.uri <- function(data.1, data.2, sig.value1="signal.value", sig.value2="signal.value", spline.df=NULL, overlap.ratio=0){

  tt <- seq(0.01, 1, by=0.01)
  vv <- tt

  if(sig.value1=="signal.value"){
    data.1.enrich <- data.frame(chr=data.1$chr, start.ori=data.1$start.ori, stop.ori=data.1$stop.ori, start=data.1$start, stop=data.1$stop, sig.value=data.1$signal.value, signal.value=data.1$signal.value, p.value=data.1$p.value, q.value=data.1$q.value)
  } else {
    if(sig.value1=="p.value"){ 
      data.1.enrich <- data.frame(chr=data.1$chr, start.ori=data.1$start.ori, stop.ori=data.1$stop.ori, start=data.1$start, stop=data.1$stop, sig.value=data.1$p.value, signal.value=data.1$signal.value, p.value=data.1$p.value, q.value=data.1$q.value)
    } else {
      if(sig.value1=="q.value"){
        data.1.enrich <- data.frame(chr=data.1$chr, start.ori=data.1$start.ori, stop.ori=data.1$stop.ori, start=data.1$start, stop=data.1$stop, sig.value=data.1$q.value, signal.value=data.1$signal.value, p.value=data.1$p.value, q.value=data.1$q.value)
      }
    }
  }

  if(sig.value2=="signal.value"){
    data.2.enrich <- data.frame(chr=data.2$chr, start.ori=data.2$start.ori, stop.ori=data.2$stop.ori, start=data.2$start, stop=data.2$stop, sig.value=data.2$signal.value, signal.value=data.2$signal.value, p.value=data.2$p.value, q.value=data.2$q.value)
  } else {
    if(sig.value2=="p.value"){ 
      data.2.enrich <- data.frame(chr=data.2$chr, start.ori=data.2$start.ori, stop.ori=data.2$stop.ori, start=data.2$start, stop=data.2$stop, sig.value=data.2$p.value, signal.value=data.2$signal.value, p.value=data.2$p.value, q.value=data.2$q.value)
    } else {
      if(sig.value2=="q.value"){
        data.2.enrich <- data.frame(chr=data.2$chr, start.ori=data.2$start.ori, stop.ori=data.2$stop.ori, start=data.2$start, stop=data.2$stop, sig.value=data.2$q.value, signal.value=data.2$signal.value, p.value=data.2$p.value, q.value=data.2$q.value)
      }
    }
  }

  ### by peaks
  # data12.enrich <- pair.peaks(data.1.enrich, data.2.enrich)
  data12.enrich <- pair.peaks.filter(data.1.enrich, data.2.enrich, p.value.impute=0, overlap.ratio)
  uri <- get.uri.2d(as.numeric(as.character(data12.enrich$merge1$sig.value)), as.numeric(as.character(data12.enrich$merge2$sig.value)), tt, vv, spline.df=spline.df)
  uri.n <- scale.t2n(uri)

  return(list(uri=uri, uri.n=uri.n, data12.enrich=data12.enrich, sig.value1=sig.value1, sig.value2=sig.value2))


}



# compute uri for matched sample
get.uri.matched <- function(data12, df=10){

  tt <- seq(0.01, 1, by=0.01)
  vv <- tt
  uri <- get.uri.2d(data12$sample1$sig.value, data12$sample2$sig.value, tt, vv, spline.df=df)

  # change scale from t to n
  uri.n <- scale.t2n(uri)

  return(list(uri=uri, uri.n=uri.n))
  
}

# map.uv is a pair of significant values corresponding to specified consistency FDR
# assuming values in map.uv and qvalue are linearly related
# data.set is the original data set
# sig.value is the name of the significant value in map.uv, say enrichment
# nominal.value is the one we want to map to, say q-value
# 
map.sig.value <- function(data.set, map.uv, nominal.value){

  index.nominal <- which(names(data.set$merge1)==nominal.value)
  nentry <- nrow(map.uv)  
  map.nominal <- rbind(map.uv[, c("sig.value1", "sig.value2")])

  for(i in 1:nentry){

    map.nominal[i, "sig.value1"] <- data.set$merge1[unique(which.min(abs(data.set$merge1$sig.value-map.uv[i, "sig.value1"]))), index.nominal]
    map.nominal[i, "sig.value2"] <- data.set$merge2[unique(which.min(abs(data.set$merge2$sig.value-map.uv[i, "sig.value2"]))), index.nominal]
  }

  invisible(map.nominal)
}


############### plot correspondence profile

# plot multiple comparison wrt one template
# uri.list contains the total number of peaks
# plot.missing=F: not plot the missing points on the right 
plot.uri.group <- function(uri.n.list, plot.dir, file.name=NULL, legend.txt, xlab.txt="num of significant peaks", ylab.txt="num of peaks in common", col.start=0, col.txt=NULL, plot.missing=F, title.txt=NULL){

  if(is.null(col.txt))
    col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")

  n <- length(uri.n.list)  

  ntotal <- c()
  for(i in 1:n)
    ntotal[i] <- uri.n.list[[i]]$ntotal

  jump.left <- c()
  jump.left.der <- c()
  ncommon <- c()
  for(i in 1:n){
#    jump.left[i]  <- which.max(uri.n.list[[i]]$uri[-1]-uri.n.list[[i]]$uri[-length(uri.n.list[[i]]$uri)])
#    if(jump.left[i] < 6)
#      jump.left[i] <- length(uri.n.list[[i]]$uri)

##  reversed.index <- seq(length(uri.n.list[[i]]$tv[,1]), 1, by=-1)
##  nequal <- sum(uri.n.list[[i]]$uri[reversed.index]== uri.n.list[[i]]$tv[reversed.index,1])
##  temp  <- which(uri.n.list[[i]]$uri[reversed.index]== uri.n.list[[i]]$tv[reversed.index,1])[nequal]
##  jump.left[i] <- length(uri.n.list[[i]]$tv[,1])-temp
##print(uri.n.list[[i]]$uri)
##print(uri.n.list[[i]]$tv[,1])
##   jump.left[i] <- uri.n.list[[i]]$jump.left

#    jump.left.der[i] <- sum(uri.n.list[[i]]$t.binned < uri.n.list[[i]]$uri.der$x[length(uri.n.list[[i]]$uri.der$x)])

    jump.left[i] <- uri.n.list[[i]]$jump.left
    jump.left.der[i] <- jump.left[i]
    ncommon[i] <- uri.n.list[[i]]$tv[jump.left[i],1]
  }


  if(plot.missing){
    max.peak <- max(ntotal)
  } else {
    max.peak <- max(ncommon)*1.05
  }

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "uri.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  }

  plot(uri.n.list[[1]]$tv[,1], uri.n.list[[1]]$uri, type="n", xlab=xlab.txt, ylab=ylab.txt, xlim=c(0, max.peak), ylim=c(0, max.peak), cex.lab=2)

  for(i in 1:n){

#    if(plot.missing){ 
#      points(uri.n.list[[i]]$tv[,1], uri.n.list[[i]]$uri, col=col.txt[i+col.start], cex=0.5 )
#    } else {
#      points(uri.n.list[[i]]$tv[1:jump.left[i],1], uri.n.list[[i]]$uri[1:jump.left[i]], col=col.txt[i+col.start], cex=0.5)
#    }
    lines(uri.n.list[[i]]$uri.spl, col=col.txt[i+col.start], lwd=4)
  }
  abline(coef=c(0,1), lty=3)
#  legend(0, max.peak, legend=legend.txt, col=col.txt[(col.start+1):length(col.txt)], lty=1, lwd=3, cex=2)
  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "duri.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  } 
  plot(uri.n.list[[1]]$t.binned, uri.n.list[[1]]$uri.slope, type="n", xlab=xlab.txt, ylab="slope", xlim=c(0, max.peak), ylim=c(0, 1.5), cex.lab=2)

  for(i in 1:n){
#    if(plot.missing){ 
#      points(uri.n.list[[i]]$t.binned, uri.n.list[[i]]$uri.slope, col=col.txt[i+col.start], cex=0.5)
#    } else {
#      points(uri.n.list[[i]]$t.binned[1:jump.left.der[i]], uri.n.list[[i]]$uri.slope[1:jump.left.der[i]], col=col.txt[i+col.start], cex=0.5)
#    }
    lines(uri.n.list[[i]]$uri.der, col=col.txt[i+col.start], lwd=4)
  }
  abline(h=1, lty=3)
#  legend(0.5*max.peak, 1.5, legend=legend.txt, col=col.txt[(col.start+1):length(col.txt)], lty=1, lwd=3, cex=2)

  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }
  
}



#######################
####################### copula fitting for matched peaks
#######################

# estimation from mixed copula model 

# 4-5-09
# A nonparametric estimation of mixed copula model


# updated

# c1, c2, f1, f2, g1, g2 are vectors
# c1*f1*g1 and c2*f2*g2 are copula densities for the two components
# xd1 and yd1 are the values of marginals for the first component
# xd2 and yd2 are the values of marginals for the 2nd component
#
# ez is the prob for being in the consistent group
get.ez <- function(p, c1, c2, xd1, yd1, xd2, yd2){

  return(p*c1*xd1*yd1/(p*c1*xd1*yd1 + (1-p)*c2*xd2*yd2))
}

# checked

# this is C_12 not the copula density function c=C_12 * f1* f2
# since nonparametric estimation is used here for f1 and f2, which
# are constant throughout the iterations, we don't need them for optimization
# 
# bivariate gaussian copula function
# t and s are vectors of same length, both are percentiles 
# return a vector
gaussian.cop.den <- function(t, s, rho){

  A <- qnorm(t)^2 + qnorm(s)^2
  B <- qnorm(t)*qnorm(s)

  loglik <-  -log(1-rho^2)/2 - rho/(2*(1-rho^2))*(rho*A-2*B)

  return(exp(loglik))
}

clayton.cop.den <- function(t, s, rho){

  if(rho > 0)
    return(exp(log(rho+1)-(rho+1)*(log(t)+log(s))-(2+1/rho)*log(t^(-rho) + s^(-rho)-1)))

  if(rho==0)
    return(1)

  if(rho<0)
    stop("Incorrect Clayton copula coefficient")
  
}


# checked
# estimate rho from Gaussian copula
mle.gaussian.copula <- function(t, s, e.z){

  # reparameterize to bound from rho=+-1
  l.c <- function(rho, t, s, e.z){
#    cat("rho=", rho, "\n")
    sum(e.z*log(gaussian.cop.den(t, s, rho)))}

  rho.max <- optimize(f=l.c, c(-0.998, 0.998), maximum=T, tol=0.00001, t=t, s=s, e.z=e.z)

#print(rho.max$m)

#cat("cor=", cor(qnorm(t)*e.z, qnorm(s)*e.z), "\t", "rho.max=", rho.max$m, "\n")
#  return(sign(rho.max$m)/(1+rho.max$m))
  return(rho.max$m)
}


# estimate mle from Clayton copula, 
mle.clayton.copula <- function(t, s, e.z){

  l.c <- function(rho, t, s, e.z){
    lc <- sum(e.z*log(clayton.cop.den(t, s, rho)))
#    cat("rho=", rho, "\t", "l.c=", lc, "\n")
    return(lc)
  }

  rho.max <- optimize(f=l.c, c(0.1, 20), maximum=T, tol=0.00001, t=t, s=s, e.z=e.z)

  return(rho.max$m)
}



# updated
# mixture likelihood of two gaussian copula
# nonparametric and ranked transformed
loglik.2gaussian.copula <- function(x, y, p, rho1, rho2, x.mar, y.mar){
 
  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
  c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)

  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}

loglik.2copula <- function(x, y, p, rho1, rho2, x.mar, y.mar, copula.txt){

  px.1 <- pdf.cdf$px.1
  px.2 <- pdf.cdf$px.2
  py.1 <- pdf.cdf$py.1
  py.2 <- pdf.cdf$py.2

  if(copula.txt=="gaussian"){
    c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
    c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  } else {
    if(copula.txt=="clayton"){
      c1 <- clayton.cop.den(px.1$cdf, py.1$cdf, rho1)
      c2 <- clayton.cop.den(px.2$cdf, py.2$cdf, rho2)
    }
  }  
  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}


# estimate the marginals of each component using histogram estimator in EM
# return the density, breaks, and cdf of the histogram estimator 
est.mar.hist <- function(x, e.z, breaks){

  binwidth <- c()
  nbin <- length(breaks)-1
  nx <- length(x) 

  # the histogram
  x1.pdf <- c()
  x2.pdf <- c()
  x1.cdf <- c()
  x2.cdf <- c()

  # the pdf for each point
  x1.pdf.value <- rep(NA, nx)
  x2.pdf.value <- rep(NA, nx)

  x1.cdf.value <- rep(NA, nx)
  x2.cdf.value <- rep(NA, nx) 

  for(i in 1:nbin){

    binwidth[i] <- breaks[i+1] - breaks[i]
    if(i < nbin)
      in.bin <- x>= breaks[i] & x < breaks[i+1]
    else    # last bin
      in.bin <- x>= breaks[i] & x <=breaks[i+1]

    # each bin add one observation to avoid empty bins
    # multiple (nx+nbin)/(nx+nbin+1) to avoid blowup when looking up for
    # quantiles 
    x1.pdf[i] <- (sum(e.z[in.bin])+1)/(sum(e.z)+nbin)/binwidth[i]*(nx+nbin)/(nx+nbin+1)        
    x2.pdf[i] <- (sum(1-e.z[in.bin])+1)/(sum(1-e.z)+nbin)/binwidth[i]*(nx+nbin)/(nx+nbin+1) 


#    x1.pdf[i] <- sum(e.z[in.bin])/sum(e.z)/binwidth[i]*nx/(nx+1)        
#    x2.pdf[i] <- sum(1-e.z[in.bin])/sum(1-e.z)/binwidth[i]*nx/(nx+1) 
    
# treat each bin as a value for a discrete variable    
#    x1.cdf[i] <- sum(x1.pdf[1:i]*binwidth[1:i])
#    x2.cdf[i] <- sum(x2.pdf[1:i]*binwidth[1:i])


    # cumulative density before reaching i
    if(i>1){
      x1.cdf[i] <- sum(x1.pdf[1:(i-1)]*binwidth[1:(i-1)])
      x2.cdf[i] <- sum(x2.pdf[1:(i-1)]*binwidth[1:(i-1)])    
    } else{
      x1.cdf[i] <- 0
      x2.cdf[i] <- 0
    }

    # make a vector of nx to store the values of pdf and cdf for each x
    # this will speed up the computation dramatically
    x1.pdf.value[in.bin] <- x1.pdf[i]
    x2.pdf.value[in.bin] <- x2.pdf[i]

    x1.cdf.value[in.bin] <- x1.cdf[i] + x1.pdf[i]*(x[in.bin]-breaks[i])
    x2.cdf.value[in.bin] <- x2.cdf[i] + x2.pdf[i]*(x[in.bin]-breaks[i])      
  }

#  x1.cdf <- cumsum(x1.pdf*binwidth)
#  x2.cdf <- cumsum(x2.pdf*binwidth)

  f1 <-list(breaks=breaks, density=x1.pdf, cdf=x1.cdf)
  f2 <-list(breaks=breaks, density=x2.pdf, cdf=x2.cdf)

  f1.value <- list(pdf=x1.pdf.value, cdf=x1.cdf.value)
  f2.value <- list(pdf=x2.pdf.value, cdf=x2.cdf.value)

  return(list(f1=f1, f2=f2, f1.value=f1.value, f2.value=f2.value))
}

# estimate the marginal cdf from rank
est.cdf.rank <- function(x, conf.z){

  # add 1 to prevent blow up
  x1.cdf <- rank(x[conf.z==1])/(length(x[conf.z==1])+1)

  x2.cdf <- rank(x[conf.z==0])/(length(x[conf.z==0])+1)

  return(list(cdf1=x1.cdf, cdf2=x2.cdf))
}

# df is a density function with fields: density, cdf and breaks, x is a scalar
get.pdf <- function(x, df){

  if(x < df$breaks[1])
    cat("x is out of the range of df\n")

  index <- which(df$breaks >= x)[1]

  if(index==1)
    index <- index +1
  return(df$density[index-1])  
}

# get cdf from histgram estimator for a single value
get.cdf <- function(x, df){

  index <- which(df$breaks >= x)[1]
  if(index==1)
    index <- index +1
  return(df$cdf[index-1])   
}

# df is a density function with fields: density, cdf and breaks
get.pdf.cdf <- function(x.vec, df){

  x.pdf <- sapply(x.vec, get.pdf, df=df)
  x.cdf <- sapply(x.vec, get.cdf, df=df) 
  return(list(cdf=x.cdf, pdf=x.pdf))
}

# E-step
# x and y are the original observations or ranks
# rho1 and rho2 are the parameters of each copula
# f1, f2, g1, g2 are functions, each is a histogram 
e.step.2gaussian <- function(x, y, p, rho1, rho2, x.mar, y.mar){

  # get pdf and cdf of each component from functions in the corresponding component 
  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
  c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  
  return(get.ez(p, c1, c2, px.1$pdf, py.1$pdf, px.2$pdf, py.2$pdf))
}

# E-step
# rho1 and rho2 are the parameters of each copula 
e.step.2copula <- function(x, y, p, rho1, rho2, x.mar, y.mar, copula.txt){

  # get pdf and cdf of each component from functions in the corresponding component 
  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  if(copula.txt=="gaussian"){
    c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
    c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  } else {
    if(copula.txt=="clayton"){
      c1 <- clayton.cop.den(px.1$cdf, py.1$cdf, rho1)
      c2 <- clayton.cop.den(px.2$cdf, py.2$cdf, rho2)
    } 
  }
  return(get.ez(p, c1, c2, px.1$pdf, py.1$pdf, px.2$pdf, py.2$pdf))
}




# M-step
m.step.2gaussian <- function(x, y, e.z, breaks){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  
  rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 

  p <- sum(e.z)/length(e.z) 

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar))
}

m.step.2copula <- function(x, y, e.z, breaks, copula.txt){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  if(copula.txt=="gaussian"){
    rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  
    rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
  } else {
    if(copula.txt=="clayton"){
      rho1 <- mle.clayton.copula(px.1$cdf, py.1$cdf, e.z)  
      rho2 <- mle.clayton.copula(px.2$cdf, py.2$cdf, 1-e.z)      
    }
  }
  
  p <- sum(e.z)/length(e.z) 

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar))
}



# E-step: pass values
# x and y are the original observations or ranks
# rho1 and rho2 are the parameters of each copula
# f1, f2, g1, g2 are functions, each is a histogram 
e.step.2gaussian.value <- function(x, y, p, rho1, rho2, pdf.cdf){

  c1 <- gaussian.cop.den(pdf.cdf$px.1$cdf, pdf.cdf$py.1$cdf, rho1)
  c2 <- gaussian.cop.den(pdf.cdf$px.2$cdf, pdf.cdf$py.2$cdf, rho2)
  
  e.z <- get.ez(p, c1, c2, pdf.cdf$px.1$pdf, pdf.cdf$py.1$pdf, 
               pdf.cdf$px.2$pdf, pdf.cdf$py.2$pdf)
  return(e.z)
}


e.step.2copula.value <- function(x, y, p, rho1, rho2, pdf.cdf, copula.txt){

  if(copula.txt =="gaussian"){
    c1 <- gaussian.cop.den(pdf.cdf$px.1$cdf, pdf.cdf$py.1$cdf, rho1)
    c2 <- gaussian.cop.den(pdf.cdf$px.2$cdf, pdf.cdf$py.2$cdf, rho2)
  } else {
    if(copula.txt =="clayton"){
      c1 <- clayton.cop.den(pdf.cdf$px.1$cdf, pdf.cdf$py.1$cdf, rho1)
      c2 <- clayton.cop.den(pdf.cdf$px.2$cdf, pdf.cdf$py.2$cdf, rho2)      
    }
  }
  
  e.z <- get.ez(p, c1, c2, pdf.cdf$px.1$pdf, pdf.cdf$py.1$pdf, 
               pdf.cdf$px.2$pdf, pdf.cdf$py.2$pdf)
  return(e.z)
}


# M-step: pass values
m.step.2gaussian.value <- function(x, y, e.z, breaks, fix.rho2){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

#  px.1 <- get.pdf.cdf(x, x.mar$f1)
#  px.2 <- get.pdf.cdf(x, x.mar$f2)
#  py.1 <- get.pdf.cdf(y, y.mar$f1)
#  py.2 <- get.pdf.cdf(y, y.mar$f2)

  px.1 <- x.mar$f1.value
  px.2 <- x.mar$f2.value
  py.1 <- y.mar$f1.value
  py.2 <- y.mar$f2.value

  rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  

  if(!fix.rho2)
    rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
  else
    rho2 <- 0

  p <- sum(e.z)/length(e.z) 

  pdf.cdf <- list(px.1=px.1, px.2=px.2, py.1=py.1, py.2=py.2)

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar,
              pdf.cdf=pdf.cdf))
}

m.step.2gaussian.value2 <- function(x, y, e.z, breaks, fix.rho2, x.mar, y.mar){

  # compute f1, f2, g1 and g2
#  x.mar <- est.mar.hist(x, e.z, breaks)
#  y.mar <- est.mar.hist(y, e.z, breaks)  

#  px.1 <- get.pdf.cdf(x, x.mar$f1)
#  px.2 <- get.pdf.cdf(x, x.mar$f2)
#  py.1 <- get.pdf.cdf(y, y.mar$f1)
#  py.2 <- get.pdf.cdf(y, y.mar$f2)

  px.1 <- x.mar$f1.value
  px.2 <- x.mar$f2.value
  py.1 <- y.mar$f1.value
  py.2 <- y.mar$f2.value

  rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  

  if(!fix.rho2)
    rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
  else
    rho2 <- 0

  p <- sum(e.z)/length(e.z) 

  pdf.cdf <- list(px.1=px.1, px.2=px.2, py.1=py.1, py.2=py.2)

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar,
              pdf.cdf=pdf.cdf))
}



m.step.2copula.value <- function(x, y, e.z, breaks, fix.rho2, copula.txt){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

#  px.1 <- get.pdf.cdf(x, x.mar$f1)
#  px.2 <- get.pdf.cdf(x, x.mar$f2)
#  py.1 <- get.pdf.cdf(y, y.mar$f1)
#  py.2 <- get.pdf.cdf(y, y.mar$f2)

  px.1 <- x.mar$f1.value
  px.2 <- x.mar$f2.value
  py.1 <- y.mar$f1.value
  py.2 <- y.mar$f2.value

  if(copula.txt=="gaussian"){
    rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  
    
    if(!fix.rho2)
      rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
    else
      rho2 <- 0
  } else {

    if(copula.txt=="clayton"){
      rho1 <- mle.clayton.copula(px.1$cdf, py.1$cdf, e.z)  
    
      if(!fix.rho2)
        rho2 <- mle.clayton.copula(px.2$cdf, py.2$cdf, 1-e.z) 
      else
        rho2 <- 0
    }    
  }
    
  p <- sum(e.z)/length(e.z) 

  pdf.cdf <- list(px.1=px.1, px.2=px.2, py.1=py.1, py.2=py.2)

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar,
              pdf.cdf=pdf.cdf))
}




# updated
# mixture likelihood of two gaussian copula
# nonparametric and ranked transformed
loglik.2gaussian.copula.value <- function(x, y, p, rho1, rho2, pdf.cdf){

  px.1 <- pdf.cdf$px.1
  px.2 <- pdf.cdf$px.2
  py.1 <- pdf.cdf$py.1
  py.2 <- pdf.cdf$py.2

  c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
  c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)

  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}



# updated
# mixture likelihood of two gaussian copula
# nonparametric and ranked transformed
loglik.2copula.value <- function(x, y, p, rho1, rho2, pdf.cdf, copula.txt){

  px.1 <- pdf.cdf$px.1
  px.2 <- pdf.cdf$px.2
  py.1 <- pdf.cdf$py.1
  py.2 <- pdf.cdf$py.2

  if(copula.txt=="gaussian"){
    c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
    c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  } else {
    if(copula.txt=="clayton"){
      c1 <- clayton.cop.den(px.1$cdf, py.1$cdf, rho1)
      c2 <- clayton.cop.den(px.2$cdf, py.2$cdf, rho2)
    }
  }

  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}



# EM for 2 Gaussian, speed up computation, unfinished

em.2gaussian.quick <- function(x, y, p0, rho1.0, rho2.0, eps, fix.p=F, stoc=T, fix.rho2=T){

  x <- rank(x, tie="average")
  y <- rank(y, tie="average")

  # nbin=20
  xy.min <- min(x, y)
  xy.max <- max(x, y)
  binwidth <- (xy.max-xy.min)/50
  breaks <- seq(xy.min-binwidth/100, xy.max+binwidth/100, by=(xy.max-xy.min+binwidth/50)/50)
#  breaks <- seq(xy.min, xy.max, by=binwidth)
  

  # initiate marginals 
  # initialization: first p0 data has 
#  e.z <- e.step.2gaussian(x, y, p0, rho1.0, rho2.0, x0.mar, y0.mar) # this starting point assumes two components are overlapped

  e.z <- c(rep(0.9, round(length(x)*p0)), rep(0.1, length(x)-round(length(x)*p0)))

  if(!stoc)
    para <- m.step.2gaussian.value(x, y, e.z, breaks, fix.rho2)
  else 
    para <- m.step.2gaussian.stoc.value(x, y, e.z, breaks, fix.rho2)


  if(fix.p){
    p <- p0
  } else {
    p <- para$p  
  }

  if(fix.rho2){
    rho2 <- rho2.0
  } else {
    rho2 <- para$rho2
  }

#  rho1 <- 0.8
  rho1 <- para$rho1

  l0 <- loglik.2gaussian.copula.value(x, y, p, rho1, rho2, para$pdf.cdf)

  loglik.trace <- c()
  loglik.trace[1] <- l0
#  loglik.trace[2] <- l1
  to.run <- T

  i <- 2

  # this two lines to remove
#  x.mar <- est.mar.hist(x, e.z, breaks)
#  y.mar <- est.mar.hist(y, e.z, breaks)  
  
  while(to.run){

    e.z <- e.step.2gaussian.value(x, y, p, rho1, rho2, para$pdf.cdf) 
    if(!stoc)
      para <- m.step.2gaussian.value(x, y, e.z, breaks, fix.rho2)
    else
      para <- m.step.2gaussian.stoc.value(x, y, e.z, breaks, fix.rho2)

    # fix x.mar and y.mar : to remove
#    if(!stoc)
#      para <- m.step.2gaussian.value2(x, y, e.z, breaks, fix.rho2, x.mar, y.mar)
#    else
#      para <- m.step.2gaussian.stoc.value(x, y, e.z, breaks, fix.rho2)

    
    if(fix.p){
      p <- p0
    } else {
      p <- para$p  
    }

    if(fix.rho2){
      rho2 <- rho2.0
    } else {
      rho2 <- para$rho2
    }

#    rho1 <- 0.8
    rho1 <- para$rho1

  #  l0 <- l1
    l1 <- loglik.2gaussian.copula.value(x, y, p, rho1, rho2, para$pdf.cdf)
    loglik.trace[i] <- l1

#cat("l1=", l1, "\n") 

    # Aitken acceleration criterion
    if(i > 2){
      l.inf <- loglik.trace[i-2] + (loglik.trace[i-1] - loglik.trace[i-2])/(1-(loglik.trace[i]-loglik.trace[i-1])/(loglik.trace[i-1]-loglik.trace[i-2])) 
      to.run <- abs(l.inf - loglik.trace[i]) > eps 
#cat("para=", "p=", para$p, " rho1=", rho1, " rho2=", rho2, "\n")
#cat("l.inf=", l.inf, "\n")
#cat(l.inf-loglik.trace[i], "\n")     
    }

    i <- i+1
  }

  bic <- -2*l1 + (2*(length(breaks)-1+1)+1-fix.p-fix.rho2)*log(length(x)) # parameters
  return(list(para=list(p=para$p, rho1=rho1, rho2=rho2), 
              loglik=l1, bic=bic, e.z=e.z, conf.z = para$conf.z, 
              loglik.trace=loglik.trace, x.mar=para$x.mar, y.mar=para$y.mar,
              breaks=breaks))
}



em.2copula.quick <- function(x, y, p0, rho1.0, rho2.0, eps, fix.p=F, stoc=T, fix.rho2=T, copula.txt, nbin=50){

  x <- rank(x, tie="first")
  y <- rank(y, tie="first")

  # nbin=50
  xy.min <- min(x, y)
  xy.max <- max(x, y)
  binwidth <- (xy.max-xy.min)/50
  breaks <- seq(xy.min-binwidth/100, xy.max+binwidth/100, by=(xy.max-xy.min+binwidth/50)/nbin)  
#  breaks <- seq(xy.min, xy.max, by=binwidth)
  
  # initiate marginals 
  # initialization: first p0 data has 
#  e.z <- e.step.2gaussian(x, y, p0, rho1.0, rho2.0, x0.mar, y0.mar) # this starting point assumes two components are overlapped

  e.z <- c(rep(0.9, round(length(x)*p0)), rep(0.1, length(x)-round(length(x)*p0)))


  if(!stoc)
    para <- m.step.2copula.value(x, y, e.z, breaks, fix.rho2, copula.txt)
  else 
    para <- m.step.2copula.stoc.value(x, y, e.z, breaks, fix.rho2, copula.txt)

  if(fix.p){
    p <- p0
  } else {
    p <- para$p  
  }

  if(fix.rho2){
    rho2 <- rho2.0
  } else {
    rho2 <- para$rho2
  }

  l0 <- loglik.2copula.value(x, y, p, para$rho1, rho2, para$pdf.cdf, copula.txt)

  loglik.trace <- c()
  loglik.trace[1] <- l0
#  loglik.trace[2] <- l1
  to.run <- T

  i <- 2

  while(to.run){

    e.z <- e.step.2copula.value(x, y, p, para$rho1, rho2, para$pdf.cdf, copula.txt) 
    if(!stoc)
      para <- m.step.2copula.value(x, y, e.z, breaks, fix.rho2, copula.txt)
    else
      para <- m.step.2copula.stoc.value(x, y, e.z, breaks, fix.rho2, copula.txt)

    if(fix.p){
      p <- p0
    } else {
      p <- para$p  
    }

    if(fix.rho2){
      rho2 <- rho2.0
    } else {
      rho2 <- para$rho2
    }


  #  l0 <- l1
    l1 <- loglik.2copula.value(x, y, p, para$rho1, rho2, para$pdf.cdf, copula.txt)
    loglik.trace[i] <- l1

cat("l1=", l1, "\n") 

    # Aitken acceleration criterion
    if(i > 2){
      l.inf <- loglik.trace[i-2] + (loglik.trace[i-1] - loglik.trace[i-2])/(1-(loglik.trace[i]-loglik.trace[i-1])/(loglik.trace[i-1]-loglik.trace[i-2])) 
      to.run <- abs(l.inf - loglik.trace[i]) > eps 
cat("para=", "p=", para$p, " rho1=", para$rho1, " rho2=", rho2, "\n")
#cat("l.inf=", l.inf, "\n")
#cat(l.inf-loglik.trace[i], "\n")     
    }

    i <- i+1
  }

  bic <- -2*l1 + (2*(length(breaks)-1+1)+1-fix.p-fix.rho2)*log(length(x)) # parameters
  return(list(para=list(p=para$p, rho1=para$rho1, rho2=rho2), 
              loglik=l1, bic=bic, e.z=e.z, conf.z = para$conf.z, 
              loglik.trace=loglik.trace, x.mar=para$x.mar, y.mar=para$y.mar,
              breaks=breaks))
}


#######################
####################### fit EM procedure for the matched peaks
#######################
rm.unmatch <- function(sample1, sample2, p.value.impute=0){

  sample1.prune <- sample1[sample1$sig.value > p.value.impute & sample2$sig.value > p.value.impute,]
  sample2.prune <- sample2[sample1$sig.value > p.value.impute & sample2$sig.value > p.value.impute,]
 
  invisible(list(sample1=sample1.prune, sample2=sample2.prune))
}


# fit 2-component model
fit.em <- function(sample12, fix.rho2=T){

  prune.sample <- rm.unmatch(sample12$merge1, sample12$merge2)

  em.fit <- em.2gaussian.quick(-prune.sample$sample1$sig.value, -prune.sample$sample2$sig.value,
 p0=0.5, rho1.0=0.7, rho2.0=0, eps=0.01, fix.p=F, stoc=F, fix.rho2)

  invisible(list(em.fit=em.fit, data.pruned=prune.sample))
}



fit.2copula.em <- function(sample12, fix.rho2=T, copula.txt){

  prune.sample <- rm.unmatch(sample12$merge1, sample12$merge2)

#  o <- order(prune.sample$sample1)
#  n <- length(prune.sample$sample1)
    
#  para <- init(prune.sample$sample1$sig.value, prune.sample$sample2$sig.value, c(rep(0, round(n/3)), rep(c(0,1), round(n/6)), rep(1, n-round(n/3)-round(n/6))))

#  temp <- init.dist(f0, f1)
  para <- list()
  para$rho <- 0.6
  para$p <- 0.3
  para$mu <- 2.5
  para$sigma <- 1
##  para$mu <- -temp$mu
##  para$sigma <- temp$sigma
#cat("mu=", para$mu, "sigma=", para$sigma, "\n")
  
#  em.fit <- em.transform.1loop(-prune.sample$sample1, -prune.sample$sample2,
  cat("EM is running")
  em.fit <- em.transform(prune.sample$sample1$sig.value, prune.sample$sample2$sig.value, para$mu, para$sigma, para$rho, para$p, eps=0.01)

  invisible(list(em.fit=em.fit, data.pruned=prune.sample))
}




# fit 1-component model
fit.1.component <- function(data.pruned, breaks){

#  gaussian.1 <- fit.gaussian.1(-data.pruned$sample1$sig.value, -data.pruned$sample2$sig.value, breaks)
#  clayton.1 <- fit.clayton.1(-data.pruned$sample1$sig.value, -data.pruned$sample2$sig.value, breaks)

  gaussian.1 <- fit.gaussian.1(-data.pruned$sample1, -data.pruned$sample2, breaks)
  clayton.1 <- fit.clayton.1(-data.pruned$sample1, -data.pruned$sample2, breaks)

  return(list(gaussian.1=gaussian.1, clayton.1=clayton.1))
}



#################
# Fit a single component  
#################

# a single gaussian copula
# if breaks=NULL, use empirical pdf, otherwise use histogram estimate
fit.gaussian.1 <- function(x, y, breaks=NULL){

  # rank transformed and compute the empirical cdf
  t <- emp.mar.cdf.rank(x)
  s <- emp.mar.cdf.rank(y)

  mle.rho <- mle.gaussian.copula(t, s, rep(1, length(t)))

  c1 <- gaussian.cop.den(t, s, mle.rho)
cat("c1", sum(log(c1)), "\n")

  if(is.null(breaks)){
    f1 <- emp.mar.pdf.rank(t)
    f2 <- emp.mar.pdf.rank(s)
  } else {
    x.mar <- est.mar.hist(rank(x), rep(1, length(x)), breaks)
    y.mar <- est.mar.hist(rank(y), rep(1, length(y)), breaks)

    f1 <- x.mar$f1.value$pdf  # only one component
    f2 <- y.mar$f1.value$pdf
  }


cat("f1", sum(log(f1)), "\n")
cat("f2", sum(log(f2)), "\n")

  loglik <- sum(log(c1)+log(f1)+log(f2))

  bic <- -2*loglik + log(length(t))*(1+length(breaks)-1)

  return(list(rho=mle.rho, loglik=loglik, bic=bic))
}


# a single Clayton copula
fit.clayton.1 <- function(x, y, breaks=NULL){

  # rank transformed and compute the empirical cdf
  t <- emp.mar.cdf.rank(x)
  s <- emp.mar.cdf.rank(y)

  mle.rho <- mle.clayton.copula(t, s, rep(1, length(t)))

  c1 <- clayton.cop.den(t, s, mle.rho)

  if(is.null(breaks)){
    f1 <- emp.mar.pdf.rank(t)
    f2 <- emp.mar.pdf.rank(s)
  } else {
    x.mar <- est.mar.hist(rank(x), rep(1, length(x)), breaks)
    y.mar <- est.mar.hist(rank(y), rep(1, length(y)), breaks)

    f1 <- x.mar$f1.value$pdf  # only one component
    f2 <- y.mar$f1.value$pdf
  }

  loglik <- sum(log(c1)+log(f1)+log(f2))

  bic <- -2*loglik + log(length(t))*(1+length(breaks)-1)

  return(list(rho=mle.rho, tau=rho/(rho+2), loglik=loglik, bic=bic)) 
}

## obsolete function (01-06-2010)
## compute the average posterior probability to belong to the random component
## for peaks selected at different cutoffs 
comp.uri.ez <- function(tt, u, v, e.z){

   u.t <- quantile(u, prob=(1-tt))
   v.t <- quantile(v, prob=(1-tt))

 #  ez <- mean(e.z[u >= u.t & v >=u.t]) Is this wrong?
   ez <- mean(e.z[u >= u.t & v >=v.t])

   return(ez)
}

## obsolete function (01-06-2010)
# compute the largest posterior error probability corresponding to
# the square centered at the origin and spanned top tt% on both coordinates
# so the consistent low rank ones are excluded
# boundary.txt: either "max" or "min", if it is error prob, use "max"
comp.ez.cutoff <- function(tt, u, v, e.z, boundary.txt){

   u.t <- quantile(u, prob=(1-tt))
   v.t <- quantile(v, prob=(1-tt))

   if(boundary.txt == "max"){
 #    ez.bound <- max(e.z[u >= u.t & v >=u.t])
     ez.bound <- max(e.z[u >= u.t & v >=v.t])
   } else {
 #    ez.bound <- min(e.z[u >= u.t & v >=u.t])
     ez.bound <- min(e.z[u >= u.t & v >=v.t])     
   }

   return(ez.bound)

}

# created: 01-06-2010
# Output IDR at various number of selected peaks
# Find cutoff (idr cutoff, sig.value cutoff on each replicate) for specified IDR level
# IDR definition is similar to FDR
get.ez.tt <- function(em.fit, idr.level=c(0.01, 0.05, 0.1)){

#  u <- em.fit$data.pruned$sample1$sig.value
#  v <- em.fit$data.pruned$sample2$sig.value
  u <- em.fit$data.pruned$sample1
  v <- em.fit$data.pruned$sample2
  
  e.z <-  1-em.fit$em.fit$e.z # this is the error prob
  
  o <- order(e.z)
  e.z.ordered <- e.z[o]
  n.select <- c(1:length(e.z))
  IDR <- cumsum(e.z.ordered)/n.select

  u.o <- u[o]
  v.o <- v[o]

  n.level <- length(idr.level)
#  sig.value1 <- rep(NA, n.level)
#  sig.value2 <- rep(NA, n.level)
  ez.cutoff <- rep(NA, n.level)
  n.selected <- rep(NA, n.level)
  
  for(i in 1:length(idr.level)){

    # find which uri.ez is closet to fdr.level
    index <- which.min(abs(IDR - idr.level[i]))
#    sig.value1[i] <- min(u.o[1:index])
#    sig.value2[i] <- min(v.o[1:index])
    ez.cutoff[i] <- e.z[index]      
    n.selected[i] <- sum(e.z<=ez.cutoff[i])    
  }   

  # output the cutoff of posterior probability, number of selected overlapped peaks 
#  map.uv <- cbind(ez.cutoff, sig.value1, sig.value2, n.selected)

  map.uv <- cbind(ez.cutoff, n.selected)

  return(list(n=n.select, IDR=IDR, idr.level=idr.level, map.uv=map.uv))
}   
  
#  return(list(n=tt*length(u), uri.ez=uri.ez,  fdr.level=fdr.level,  map.uv=map.uv, e.z=e.z, error.prob.cutoff=error.prob.cutoff))  
  




### compute the mean of the marginals
get.mar.mean <- function(em.out){

  x.f1 <- em.out$x.mar$f1
  x.f2 <- em.out$x.mar$f2

  y.f1 <- em.out$y.mar$f1
  y.f2 <- em.out$y.mar$f2

  x.stat1 <- get.hist.mean(x.f1)
  x.stat2 <- get.hist.mean(x.f2)
  y.stat1 <- get.hist.mean(y.f1)
  y.stat2 <- get.hist.mean(y.f2)

  return(list(x.mean1=x.stat1$mean, x.mean2=x.stat2$mean, 
              y.mean1=y.stat1$mean, y.mean2=y.stat2$mean,
              x.sd1=x.stat1$sd, x.sd2=x.stat2$sd, 
              y.sd1=y.stat1$sd, y.sd2=y.stat2$sd
              ))

}


# compute the mean of marginals
get.hist.mean  <- function(x.f){

  nbreaks <- length(x.f$breaks)
  x.bin <- x.f$breaks[-1]-x.f$breaks[-nbreaks]

  x.mid <- (x.f$breaks[-nbreaks]+x.f$breaks[-1])/2
  x.mean <- sum(x.mid*x.f$density*x.bin)
  x.sd <- sqrt(sum(x.mid*x.mid*x.f$density*x.bin)-x.mean^2)
  
  return(list(mean=x.mean, sd=x.sd))
}

get.hist.var <- function(x.f){

  nbreaks <- length(x.f$breaks)
  x.bin <- x.f$breaks[-1]-x.f$breaks[-nbreaks]

  x.mid <- (x.f$breaks[-nbreaks]+x.f$breaks[-1])/2
  x.mean <- sum(x.mid*x.f$density*x.bin)

  return(mean=x.mean)  
}

plot.ez.group <- function(ez.list, plot.dir, file.name=NULL, legend.txt, y.lim=NULL, xlab.txt="num of significant peaks",  ylab.txt="IDR", col.txt=NULL, title.txt=NULL){

  if(is.null(col.txt))
    col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")

  n.entry <- length(ez.list)
  x <- rep(NA, n.entry)
  y.max <- rep(NA, n.entry)

  for(i in 1:n.entry){
    x[i] <- max(ez.list[[i]]$n)
      
    y.max[i] <- max(ez.list[[i]]$IDR)
  
  }

  if(is.null(y.lim))
    y.lim <- c(0, max(y.max))

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "ez.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  }


  
  plot(c(0, max(x)), y.lim, ylim=y.lim, type="n", xlab=xlab.txt, ylab=ylab.txt, lwd=4, cex=5, cex.lab=2)

  q <- seq(0.01, 0.99, by=0.01)
  
  for(i in 1:length(ez.list)){

    n.plot <- round(quantile(ez.list[[i]]$n, prob=q))
    IDR.plot <- ez.list[[i]]$IDR[n.plot]
    lines(n.plot, IDR.plot, col=col.txt[i], cex=2, lwd=5)    
  }


#  legend(0, y.lim[2], legend=legend.txt, col=col.txt[1:length(col.txt)], lty=1, lwd=5, cex=2)

  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }
  
}



#############################################################################
#############################################################################
# statistics about peaks selected on the individual replicates
#
# idr.level: the consistency cutoff, say 0.05
# uri.output: a list of uri.output from consistency analysis generated by batch-consistency-analysis.r
# ez.list : a list of IDRs computed from get.ez.tt using the same idr.level
#
##################

get.ez.tt.all <- function(em.fit, all.data1, all.data2, idr.level=c(0.01, 0.05, 0.1)){

  u <- em.fit$data.pruned$sample1$sig.value
  v <- em.fit$data.pruned$sample2$sig.value
#  u <- em.fit$data.pruned$sample1
#  v <- em.fit$data.pruned$sample2
  
  e.z <-  1-em.fit$em.fit$e.z # this is the error prob
  
  o <- order(e.z)
  e.z.ordered <- e.z[o]
  n.select <- c(1:length(e.z))
  IDR <- cumsum(e.z.ordered)/n.select

  u.o <- u[o]
  v.o <- v[o]

  n.level <- length(idr.level)
#  sig.value1 <- rep(NA, n.level)
#  sig.value2 <- rep(NA, n.level)
  ez.cutoff <- rep(NA, n.level)
  n.selected <- rep(NA, n.level)
  npeak.rep1 <- rep(NA, n.level)
  npeak.rep2 <- rep(NA, n.level)
  
  for(i in 1:length(idr.level)){

    # find which uri.ez is closet to fdr.level
    index <- which.min(abs(IDR - idr.level[i]))
#    sig.value1[i] <- min(u.o[1:index])
#    sig.value2[i] <- min(v.o[1:index])
    ez.cutoff[i] <- e.z.ordered[index]      # fixed on 02/20/10
    n.selected[i] <- sum(e.z<=ez.cutoff[i])
#    npeak.rep1[i] <- sum(all.data1["sig.value"] >= sig.value1[i])
#    npeak.rep2[i] <- sum(all.data2["sig.value"] >= sig.value2[i])     
  }   

  # output the cutoff of posterior probability, number of selected overlapped peaks 
  map.uv <- cbind(ez.cutoff, n.selected)

  return(list(n=n.select, IDR=IDR, idr.level=idr.level, map.uv=map.uv))
}   
  
#  return(list(n=tt*length(u), uri.ez=uri.ez,  fdr.level=fdr.level,  map.uv=map.uv, e.z=e.z, error.prob.cutoff=error.prob.cutoff))  
  





####### the following is for determining thresholds for merged dataset


###### new fitting procedure




# 1. rank pairs

# 2. initialization
# 3. convert to pseudo-number

# 4. EM

# need plugin and test
# find the middle point between the bins
get.pseudo.mix <- function(x, mu, sigma, rho, p){

  
  # first compute cdf for points on the grid
  # generate 200 points between [-3, mu+3*sigma]
  nw <- 1000
  w <- seq(min(-3, mu-3*sigma), max(mu+3*sigma, 3), length=nw) 
  w.cdf <- p*pnorm(w, mean=mu, sd=sigma) + (1-p)*pnorm(w, mean=0, sd=1)

  i <- 1

  quan.x <- rep(NA, length(x))

  for(i in c(1:nw)){
    index <- which(x >= w.cdf[i] & x < w.cdf[i+1])
    quan.x[index] <- (x[index]-w.cdf[i])*(w[i+1]-w[i])/(w.cdf[i+1]-w.cdf[i]) +w[i]
  }

  index <- which(x < w.cdf[1])
  if(length(index)>0)
    quan.x[index] <- w[1]

  index <- which(x > w.cdf[nw])
  if(length(index)>0)
    quan.x[index] <- w[nw]  
  
#  linear.ext <- function(x, w, w.cdf){
  # linear interpolation
#    index.up <- which(w.cdf>= x)[1]
#    left.index <- which(w.cdf <=x)
#    index.down <- left.index[length(left.index)]
#    quan.x <- (w[index.up] + w[index.down])/2  
#  }
  
#  x.pseudo <- sapply(x, linear.ext, w=w, w.cdf=w.cdf)

#  invisible(x.pseudo)
  invisible(quan.x)
}


# EM to compute the latent structure
# steps:
# 1. raw values are first transformed into pseudovalues
# 2. EM is used to compute the underlining structure, which is a mixture
#    of two normals
em.transform <- function(x, y, mu, sigma, rho, p, eps){
  
  x.cdf.func <- ecdf(x)
  y.cdf.func <- ecdf(y)
  afactor <- length(x)/(length(x)+1)
  x.cdf <- x.cdf.func(x)*afactor
  y.cdf <- y.cdf.func(y)*afactor
  
  # initialization
  para <- list()
  para$mu <- mu
  para$sigma <- sigma
  para$rho <- rho
  para$p <- p  

  j <- 1
  to.run <- T
  loglik.trace <- c()
  loglik.inner.trace <- c()
  
  #to.run.inner <- T
  z.1 <- get.pseudo.mix(x.cdf, para$mu, para$sigma, para$rho, para$p)
  z.2 <- get.pseudo.mix(y.cdf, para$mu, para$sigma, para$rho, para$p)

#  cat("length(z1)", length(z.1), "\n")
  while(to.run){
    
    # get pseudo value in each cycle
#    z.1 <- get.pseudo.mix(x.cdf, para$mu, para$sigma, para$rho, para$p)
#    z.2 <- get.pseudo.mix(y.cdf, para$mu, para$sigma, para$rho, para$p)

    i <- 1
    while(to.run){
      
      # EM for latent structure
      e.z <- e.step.2normal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)
      para <- m.step.2normal(z.1, z.2, e.z)
#para$rho <- rho
#para$p <- p    
#para$mu <- mu
#para$sigma <- sigma    
      if(i > 1)
        l.old <- l.new
    
      # this is just the mixture likelihood of two-component Gaussian
      l.new <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)

      loglik.inner.trace[i] <- l.new 

      if(i > 1){
        to.run <- loglik.inner.trace[i]-loglik.inner.trace[i-1]>eps         
      }
        
    
#      if(i > 2){
#        l.inf <- loglik.inner.trace[i-2] + (loglik.inner.trace[i-1] - loglik.inner.trace[i-2])/(1-(loglik.inner.trace[i]-loglik.inner.trace[i-1])/(loglik.inner.trace[i-1]-loglik.inner.trace[i-2]))

#        if(loglik.inner.trace[i-1]!=loglik.inner.trace[i-2])
#          to.run <- abs(l.inf - loglik.inner.trace[i]) > eps
#        else
#          to.run <- F
          
#      }

      cat("loglik.inner.trace[", i, "]=", loglik.inner.trace[i], "\n")
    cat("mu=", para$mu, "sigma=", para$sigma, "p=", para$p, "rho=", para$rho, "\n\n")
      
      i <- i+1
    }
    

    # get pseudo value in each cycle
    z.1 <- get.pseudo.mix(x.cdf, para$mu, para$sigma, para$rho, para$p)
    z.2 <- get.pseudo.mix(y.cdf, para$mu, para$sigma, para$rho, para$p)

    if(j > 1)
      l.old.outer <- l.new.outer

    l.new.outer <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)

    loglik.trace[j] <- l.new.outer
    
    if(j == 1)
      to.run <- T
    else{ # stop when iteration>100
      if(j > 100)
        to.run <- F
      else
        to.run <- l.new.outer - l.old.outer > eps
    }

#    if(j %% 10==0)
      cat("loglik.trace[", j, "]=", loglik.trace[j], "\n")
    cat("mu=", para$mu, "sigma=", para$sigma, "p=", para$p, "rho=", para$rho, "\n")
    
    j <- j+1
  }

  bic <- -2*l.new + 4*log(length(z.1))
  
  return(list(para=list(p=para$p, rho=para$rho, mu=para$mu, sigma=para$sigma),
              loglik=l.new, bic=bic, e.z=e.z, loglik.trace=loglik.trace))
}  




# compute log-likelihood for mixture of two bivariate normals
loglik.2binormal <- function(z.1, z.2, mu, sigma, rho, p){

  l.m <- sum(d.binormal(z.1, z.2, 0, 1, 0)+log(p*exp(d.binormal(z.1, z.2, mu, sigma, rho)-d.binormal(z.1, z.2, 0, 1, 0))+(1-p)))
  
#  l.m <- sum((p*d.binormal(z.1, z.2, mu, sigma, rho) + (1-p)*d.binormal(z.1, z.2, 0, 1, 0)))
  return(l.m) 
}

# check this when rho=1

# density of binomial distribution with equal mean and sigma on both dimensions
d.binormal <- function(z.1, z.2, mu, sigma, rho){

  loglik <- (-log(2)-log(pi)-2*log(sigma) - log(1-rho^2)/2 - (0.5/(1-rho^2)/sigma^2)*((z.1-mu)^2 -2*rho*(z.1-mu)*(z.2-mu) + (z.2-mu)^2))

  return(loglik)
}

# E-step for computing the latent strucutre
# e.z is the prob to be in the consistent group
# e.step for estimating posterior prob
# z.1 and z.2 can be vectors or scalars
e.step.2normal <- function(z.1, z.2, mu, sigma, rho, p){

  e.z <- p/((1-p)*exp(d.binormal(z.1, z.2, 0, 1, 0)-d.binormal(z.1, z.2, mu, sigma, rho))+ p)
  
  invisible(e.z)
}

# M-step for computing the latent structure
# m.step for estimating proportion, mean, sd and correlation coefficient
m.step.2normal <- function(z.1, z.2, e.z){

  p <- mean(e.z)
  mu <- sum((z.1+z.2)*e.z)/2/sum(e.z) 
  sigma <- sqrt(sum(e.z*((z.1-mu)^2+(z.2-mu)^2))/2/sum(e.z))
  rho <- 2*sum(e.z*(z.1-mu)*(z.2-mu))/(sum(e.z*((z.1-mu)^2+(z.2-mu)^2)))

  return(list(p=p, mu=mu, sigma=sigma, rho=rho))
}


# assume top p percent of observations are true
# x and y are ranks, estimate
init <- function(x, y, x.label){
  
  x.o <- order(x)

  x.ordered <- x[x.o]
  y.ordered <- y[x.o]
  x.label.ordered <- x.label[x.o]
  
  n <- length(x)
  p <- sum(x.label)/n
  
  rho <- cor(x.ordered[1:ceiling(p*n)], y.ordered[1:ceiling(p*n)])

  temp <- find.mu.sigma(x.ordered, x.label.ordered)
  mu <- temp$mu
  sigma <- temp$sigma
  
  invisible(list(mu=mu, sigma=sigma, rho=rho, p=p))

}

# find mu and sigma if the distributions of marginal ranks are known
# take the medians of the two dist and map back to the original
init.dist <- function(f0, f1){

  # take the median in f0
  index.median.0 <- which(f0$cdf>0.5)[1]
  q.0.small <- f0$cdf[index.median.0] # because f0 and f1 have the same bins
  q.1.small <- f1$cdf[index.median.0]

  # take the median in f1
  index.median.1 <- which(f1$cdf>0.5)[1]
  q.0.big <- f0$cdf[index.median.1] # because f0 and f1 have the same bins
  q.1.big <- f1$cdf[index.median.1]

  # find pseudo value for x.middle[1] on normal(0,1) 
  pseudo.small.0 <- qnorm(q.0.small, mean=0, sd=1)
  pseudo.small.1 <- qnorm(q.1.small, mean=0, sd=1)

  # find pseudo value for x.middle[2] on normal(0,1) 
  pseudo.big.0 <- qnorm(q.0.big, mean=0, sd=1)
  pseudo.big.1 <- qnorm(q.1.big, mean=0, sd=1)

  mu <- (pseudo.small.0*pseudo.big.1 - pseudo.small.1*pseudo.big.0)/(pseudo.big.1-pseudo.small.1) 

  sigma <- (pseudo.small.0-mu)/pseudo.small.1

  return(list(mu=mu, sigma=sigma))  
}

# generate labels

# find the part of data with overlap

# find the percentile on noise and signal

# Suppose there are signal and noise components, with mean=0 and sd=1 for noise
# x and x.label are the rank of the observations and their labels,
# find the mean and sd of the other component
# x.label takes values of 0 and 1
find.mu.sigma <- function(x, x.label){

  x.0 <- x[x.label==0]
  x.1 <- x[x.label==1]

  n.x0 <- length(x.0)
  n.x1 <- length(x.1)

  x.end <- c(min(x.0), min(x.1), max(x.0), max(x.1))
  o <- order(x.end)
  x.middle <- x.end[o][c(2,3)]

  # the smaller end of the overlap
  q.1.small <- mean(x.1 <= x.middle[1])*n.x1/(n.x1+1)
  q.0.small <- mean(x.0 <= x.middle[1])*n.x0/(n.x0+1)

  # the bigger end of the overlap
  q.1.big <- mean(x.1 <= x.middle[2])*n.x1/(n.x1+1)
  q.0.big <- mean(x.0 <= x.middle[2])*n.x0/(n.x0+1)

  # find pseudo value for x.middle[1] on normal(0,1) 
  pseudo.small.0 <- qnorm(q.0.small, mean=0, sd=1)
  pseudo.small.1 <- qnorm(q.1.small, mean=0, sd=1)

  # find pseudo value for x.middle[2] on normal(0,1) 
  pseudo.big.0 <- qnorm(q.0.big, mean=0, sd=1)
  pseudo.big.1 <- qnorm(q.1.big, mean=0, sd=1)

  mu <- (pseudo.small.0*pseudo.big.1 - pseudo.small.1*pseudo.big.0)/(pseudo.big.1-pseudo.small.1) 

  sigma <- (pseudo.small.0-mu)/pseudo.small.1

  return(list(mu=mu, sigma=sigma))
}
