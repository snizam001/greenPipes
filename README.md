# greenPipe
     Flaxible, scalable and integrative analysis of epigenomics, proteomics and transcriptomic datasets. 
     See license before cloning/using this tool

# citation: 
     
     Unpublished work

INSTALLATION:
       
       git clone https://github.com/snizam001/greenPipe.git
       cd greenPipe
       pip install ./
       
USAGE:
     
     For help type:
     greenPipe --help
     
     greenPipe [options]

     optional arguments:
       -h, --help            show this help message and exit
       --effectiveGenomeSize EFFECTIVEGENOMESIZE
                             callPeaks, idr, coverageTracks, initHeatmap mode: effective genome size. human (2913022398: hg38, 2864785220: hg19), mouse (2652783500: GRCm38, 2620345972: GRCm37) and fruitfly (142573017:dm6, 162367812: dm3). (default is 2913022398 which is equivalent for the human
                             genome)
       --refgenome REFGENOME
                             alignment mode: bowtie2 index file of reference genome. For example: if bowtie2 indexes is present in folder /home/databases/genomes and name of indexes are: GRCh38.p13.1.bt2, GRCh38.p13.2.bt2 etc; then give path /home/databases/genomes/GRCh38.p13
       --spikein SPIKEIN     alignment mode: bowtie2 index file of spikeIn (Drosophilla). Give path of index as suggested for --referenceGenome
       --alignParam ALIGNPARAM
                             alignment mode: The alignment parameters of the bowtie2/bwa. The default parameter of bowtie2 in this pipeline is: --dovetail --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700 (based on Meers et al. (2019)). For single-end default
                             parameters of the program of bwa was used.If you want to change it, give parameter as comma seperated values e.g. for bowtie2 --no-unal,--no-mixed,--no-discordant,-I,0,-X,1500 and for bwa see the manual of bwa (for mem). If --gpu is True, then see the options of the nvBowtie:
                             https://nvlabs.github.io/nvbio/nvbowtie_page.html
       --gpu {True,False}    If you have good source of the GPU. Use this option to activate. Wherever neccessary tool automatically recognize and use the GPU source.
       --SelectReads SELECTREADS
                             equalRead mode: How many reads should be use randomly? If not given minimum number of read will be choosed using all experimental files
       --spikeNormPeak {True,False}
                             initpeakcalling mode: Should include SpikeIn normalizationin peak calling ? (default is True)
       --reverseName_equalRead {True,False}
                             equalRead mode: If the bamfile folder already contains *.original.bam files and you are selecting random reads, then mode equalRead will denote error. Therefore either manually rename *.original.bam to *.bam or use this function to automatically rename files.
       --blackListedRegions BLACKLISTEDREGIONS
                             qcExperiment, qcBamfiles, callPeaks, idr, coverageTracks, initHeatmap, heatmap mode: regions which are black listed by ENCODE in bed format.
       --pMethod {homer,macs2,seacr}
                             callPeak mode: peak calling algorithm. Default is homer
       --pStyle PSTYLE       callPeak mode: style of the peaks (narrow/broad/both). Default is the narrow. For each samples users can provide comma seperated stylese.g. both,narrow,narrow,broad
       --pFdrHomer PFDRHOMER
                             callPeak mode: FDR/poisson in the homer peak calling. Default 0.001
       --pPvalueHomer PPVALUEHOMER
                             callPeak mode: p value in the homer peak calling. Default is 0.0001
       --pFcHomer PFCHOMER   callPeak mode: Fold changes in the homer peak calling. Default is 4.0
       --pDistHomer {fdr,poisson}
                             callPeak mode: Distribution (FDR or poisson)
       --pControl {True,False}
                             callPeak mode: Include control in the homer peak calling? Default is True
       --pSpike {True,False}
                             callPeak mode: Normalize reads by spike-in in the homer peak calling. Default is True
       --pSeacrMode {stringent,relaxed}
                             callPeak mode: Mode when calling peaks using SEACR. Default is strigent
       --pSeacrThreshold PSEACRTHRESHOLD
                             callPeak mode: Threshold value if control is not used in peak calling using SEACR. Default is 0.01
       --pOpts POPTS         callPeak mode: Pipeline uses all default parameters. If you want to change something, in the peakcalling from homer and macs2 you can give here as comma seperated values
       --genomeFile GENOMEFILE
                             callPeak mode: With SEACR mode --genomeFile is require which contains the length of each chromosome. These files can be downloaded from https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes.Preferebly generate your own genome file from fasta (used in the
                             alignment) by using samtools faidx genome.fa && cut -f1,2 genome.fa.fai > genomeFile.txt. By default program uses hg38 genome file (hg38.genome)
       --idrExprs IDREXPRS   idr mode: List of the tagDirectories of the experiments (IDR mode). If --idrMethod is homer then, give as follows: /home/dir1/exp1_rep1,/home/dir1/exp1_rep2;/home/dir1/exp2_rep1,/home/dir1/exp2_rep2,/home/dir1/exp2_rep3. Give full path. If --idrMethod is macs2 then
                             give list of the bamfiles
       --idrCtrl IDRCTRL     idr mode: List of the tagDirectories of the control (IDR mode). If --idrMethod is homer then, give as follows: /home/dir1/ctrl1_rep1,/home/dir1/ctrl1_rep2;/home/dir1/ctrl2_rep1,/home/dir1/ctrl2_rep2,/home/dir1/ctrl2_rep3. Give full path. If --idrMethod is macs2 then
                             give list of the bamfiles
       --idrName IDRNAME     idr mode: List of the Name of the experiment (IDR mode). Give as follows: expr1;expr2
       --idrExprSpike IDREXPRSPIKE
                             idr mode: List of the experiment bamfiles of the spikeIn, if --idrSpike is True. Format is same. Give full path
       --idrCtrlSpike IDRCTRLSPIKE
                             idr mode: List of the control bamfiles of the spikeIn, if --idrSpike is True. Format is same. Give full path
       --idrControl {True,False}
                             idr mode: Should include the control in IDR peak calling (IDR mode). Default is True.
       --idrSpike {True,False}
                             idr mode: Should include the spikeIn in IDR peak calling (IDR mode). Default is True
       --idrStyle {factor,histone}
                             idr mode: narrow or broad peaks (IDR mode). Default is narrow (factor) or broad (histone)You can give a single style for all experiments e.g. factor or different style for each experiments e.g. factor;histone
       --idrOutput IDROUTPUT
                             idr mode: Prefix of the output file for each experiment (IDR mode) e.g. idr_expr1;idr_Expr2
       --idrMethod {homer,macs2}
                             idr mode: Peak calling method: homer or macs2 in IDR
       --overFiles OVERFILES
                             PeakComparison mode: Input peak file for the comparison. Give as CSV input e.g. dir1/dir2/peak1.bed,dir1/dir2/peak2.bed. If not give, pipeline automatically identify the IDR peaks or peaks from the output folder. Give the full path of file.
       --overDist OVERDIST   PeakComparison mode: The distance between peaks to be consider as overlapping. Default is 400 base-pairs
       --annpeakFiles ANNPEAKFILES
                             UserAnnotation/annotation mode: in this mode you can provide peaks which can be annotated by user provided annotations files (--annFiles). Multiple peak files can be given as comma seperated values.
       --annFiles ANNFILES   UserAnnotation mode: provide bed files of your region of interest e.g. h3k4me3, jun/fos motif location etc. Users can provide more than one annotation bed files as comma seperated values e.g. f1.bed,f2.bed
       --annSize ANNSIZE     UserAnnotation mode: maximum distance between the boundary of peaks (--annpeakFiles) and annotation bed files (annFiles). For each annotation bed files provide this maximum distance e.g. for h3k4me3 500, for jun/fos 100 as 500,100
       --annName ANNNAME     UserAnnotation mode: give name of your annotations e.g. give h3k4me3,jun/fos for file f1.bed,f2.bed
       --annPrefix ANNPREFIX
                             UserAnnotation/annotation mode: prefix of the output annotated files. Use this option if peak information is given by users by --annpeakFiles
       --covSpike {True,False}
                             coverageTracks mode: During calculation of the whole genome coverage in bins, reads should be normalize according to spikeIn. This normalization will be calculated within the samples given in the --inputfile. Default is False
       --covOtherOptions COVOTHEROPTIONS
                             coverageTracks mode: Pipeline uses all default parameters. If you want to change something, besides --bl, --effectiveGenomeSize and -p, the other options of bamCoverage (deepTools) can be provided here as comma seperated values.
       --cMotif {True,False}
                             cutfrequency mode: cutfrequency mode needs location of the motifs in the whole genome If you do not have this file then create it by using --cMotif True
       --cMotifFile CMOTIFFILE
                             cutfrequency mode: If the --cMotif is False, then give the file which contains the location of the motifs in the whole genome. It should be in HOMER bed file format
       --cMaN CMAN           cutfrequency mode: Provide the MA number of the JASPER motif Go to http://jaspar.genereg.net/ for finding this number e.g. for TP53 number is MA0106
       --cGVersion CGVERSION
                             cutfrequency/annotation: In the cutfrequency mode, if --cMotif is true, provide the version of the genome, so that location of the motifs can be generated. This one is also require in the annotation mode
       --cCenter CCENTER     cutfrequency mode: Give the center of the motif
       --initHeatmapOtherOptions INITHEATMAPOTHEROPTIONS
                             initHeatmap mode: Pipeline uses all default parameters. If you want to change something, besides --bl, --effectiveGenomeSize and -p, the other options of bamCompare (deepTools) can be provided here as comma seperated values.
       --hMOpt HMOPT         heatmap mode: Pipeline uses all default parameters. If you want to change something, besides --missingDataAsZero,-bl,--smartLabels,-p,--metagene,--samplesLabel, the other options of computeMatrix (deepTools) can be provided here as comma seperated values.(moreover,
                             by default -a 5000, -b 5000 are used. It can be changed here.)
       --hPOpt HPOPT         heatmap mode: Pipeline uses all default parameters. If you want to change something, besides --refPointLabel and --dpi, the other options of plotHeatmap (deepTools) can be provided here as comma seperated values.(moreover, by default--colorMap RdBu are used. It can
                             be changed here.)
       --hSpike {True,False}
                             heatmap mode: Should include spike in preparation of bamcompare files. Default is True
       --hCovComp {compare,coverage,NA}
                             heatmap mode: Heatmap should be plotted for bamcoverage (spikeIn normalized) or bamcompare files? Leave it, if you want to given your own files using --hInFiles
       --hInFiles HINFILES   heatmap mode: Input bamcoverage or bamcompare files in csv format. If file is not given, tool will indentify files from the --inputfile
       --hInNames HINNAMES   heatmap mode: Name of the samples in csv format. Equal to the number of files in --hInFiles.
       --hRegionMode {tss,metagene,bed,peaks}
                             heatmap mode: mode of the region? on tss, on given peak file or select automatically total peak file or differential peaks from Peak folder or metagene. Default is peaks mode.
       --hGtf HGTF           heatmap mode: Path of the GTF file, if --hRegionMode is tss or metagene
       --hBed HBED           heatmap mode: List of bed files on which you want to generate heatmap in csv format, if --hRegionMode bed
       --hDiffPeaks {True,False}
                             heatmap mode: In the --hRegionMode peaks, should program consider differential peaks also ? Default is the True
       --lProt LPROT         piggyBack mode: Files with name of the proteins per line. These are those proteins which are significantly enriched in the mass spectrometry experiment and you want to check if these proteins are providing any piggy-back binding to your protein of interest. You can
                             find the list of proteins on basis of the stochiometry (iBaq values) observed in your MassSpectrometry experiment. This mode has limitation because it can predict the piggy-back binding events only for those proteins which are binding to the DNA (not histone modifications
                             proteins) and their PWM are known. This list should not contain any special character e.g. ; or - or :
       --mPwm MPWM           piggyBack and doughnut mode: Provide the MA number of the JASPER motif Go to http://jaspar.genereg.net/ for finding this number e.g. for TP53 number is MA0106. You can provide more than 1 motifs as comma seperated values.
       --sProt SPROT         piggyBack and doughnut mode: provide the sequence of protein of interest for which you are searching piggy-back binding event e.g. for NFYA you can give CCAAT, for JUN you can give TGANTCA, for ATF7 you can give TGANNTCA. You can provide more than one sequence using
                             comma seperated values. In the sequence IUPAC codes are also applicable.
       --mPeak MPEAK         piggyBack: Path of the peak file of greenPipe experiment in bed format. It is optional. If not given, tool will automatically search peak file on basis of your --outputdir and --inputfile. If provided give the full path of peaks as comma seperated values e.g.
                             /home/xyz/exp1.Clean.bed,/home/xyz/exp2.Clean.bed.
       --mDist MDIST         piggyBack and doughnut mode: Search motifs from this distance from the centre of peaks. Default is 400
       --distPiggy DISTPIGGY
                             piggyBack and doughnut mode: Distance to find piggy back binding events. Default is 400
       --Species {human,mouse,rat,fruitfly,nematode,zebrafish,thale-cress,frog,pig}
                             piggyBack and doughnut mode: Taxon ID or name of Species according to NDBI e.g. for home sapiens taxon ID is 9606 and name is human
       --sFasta SFASTA       piggyBack/doughnut/annotation mode: Fasta file of the species genome
       --gVer GVER           piggyBack and doughnut mode: Genome version of the fasta file. Default is hg38
       --mPvalue MPVALUE     piggyBack mode: Motif finding pvalue for piggy back binding events. Default is 0.05
       --sDist SDIST         piggyBack and doughnut mode: Search sequence of the motif from the center of peak. Default is 200. best practice: --sDist should be half of the --mDist
       --mPrefix MPREFIX     piggyBack mode: If giving your own bed peak files using --mPeak then give prefix of outputfor each bedfiles
       --dPeakfiles DPEAKFILES
                             doughnut mode: Path of the peak file for doughnut mode in bed format as comma seperated values e.g. /home/xyz/exp1.Clean.bed,/home/xyz/exp2.Clean.bed.
       --dNames DNAMES       doughnut mode: name of the experiment for each bed file given with --dPeakfiles e.g. exp1,exp2.

     required arguments
        _____________________:
       --modes MODES         run mode. Multiple modes can be provided as comma seperated values e.g. --modes qc,alignment. Choices are: qc, contamination, alignment, equalRead, qcExperiment, initPeakCalling, qcTagDirectories, callPeaks, idr, doughnut, PeakComparison, annotation,
                             UserAnnotation, coverageTracks, cutfrequency, initHeatmap, heatmap, piggyBack, rnaIntegrate, stats
       --outputdir OUTPUTDIR
                             output directory (provide full path)

     common arguments
        _____________________:
       --inputdir INPUTDIR   input directory having fastq files (provide full path)
       --inputfile INPUTFILE
                             input file in .txt format
       --libraryType {single,pair}
                             type of the library
       --threads THREADS     number of threads (default is: total CPU - 2)

      Authors: (a) Sheikh Nizamuddin: snizam001@gmail.com, (b) H.Th. Marc Timmers: m.timmers@dkfz-heidelberg.de 


     #--------------------------------------------------------------------------

     Format of the SampleInfo.txt
     
     Control_R1     Control_R2     Experiment_R1  Experiment_R2  Prefix/Name     
     IgG-neg_R1.fastq.gz	IgG-neg_R2.fastq.gz	TBP2-neg_R1.fastq.gz	TBP2-neg_R2.fastq.gz	TBP_neg
     
     Example command to run:
     
     1. Prepare the SampleInfo.txt
     2. Create Bowtie2 index of mm10 mouse genome index
     3. Run this command: 
     
     greenPipe --modes qc,alignment,initPeakCalling,
     callPeaks,annotation,coverageTracks,initHeatmap,heatmap 
     --inputdir /media/sheikh/Swift0/dTAGSystemTBP/CUTRUN2/Fastq 
     --inputfile ./SampleInfo.txt --threads 30 --outputdir 
     /media/sheikh/Swift0/dTAGSystemTBP/CUTRUN2 
     --refgenome /media/txpn/nvme/Databases/mm10/mm10 
     --spikein /media/txpn/nvme/Databases/Drosophilla/Bowtie2/
     Drosophila_melanogaster --blackListedRegions mm10.blacklist.bed 
     --pStyle narrow --covSpike True --hCovComp coverage 
     --hDiffPeaks False 
     
