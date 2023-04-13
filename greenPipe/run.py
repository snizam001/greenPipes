#!/usr/bin/env python3

import pkg_resources
import sys
import subprocess
import argparse
import os
import pandas as pd
import re
import glob
from itertools import combinations
from termcolor import colored
from datetime import datetime
from itertools import product
import itertools
import mygene
mg = mygene.MyGeneInfo()
from Bio.Seq import Seq
import pyfaidx
import pkg_resources as psource
import filetype
from multiprocessing import Pool
from itertools import repeat
from greenPipe import universal
from greenPipe import linksFile as lf
from greenPipe import qcFastQ
from greenPipe import contamination
from greenPipe import align
from greenPipe import equalRead
from greenPipe import qcBamfiles
from greenPipe import initPeakCalling
from greenPipe import qcTagD
from greenPipe import ann
from greenPipe import callPeak
from greenPipe import covTracks
from greenPipe import greenCutFrq
from greenPipe import heatMap
from greenPipe import idr
from greenPipe import filterMotifSeq
from greenPipe import massSpectro
from greenPipe import comparePeak
from greenPipe import stats
from greenPipe import rna
from greenPipe import minReadPrediction
from greenPipe import hugo2Entrez

#from greenPipe import generateIndex
#_____________________________________________________________________________________________________

__description__ = """

greenPipe pipeline (version 3.0): 23rd January, 2023

#--note: for modes using homer, first customize homer if you are using species other than human. See http://homer.ucsd.edu/homer/introduction/update.html
#--All fastqfiles must be in zip Format


1. qc mode
    greenPipe --modes qc --outputdir $(pwd) --inputdir $(pwd)/Fastq --inputfile sampleInfo.txt --libraryType pair

2. alignment mode
    greenPipe --modes alignment --outputdir $(pwd) --inputdir $(pwd)/Fastq --inputfile sampleInfo.txt --libraryType pair
    --refgenome /home/sheikh/Databases/hg38/GenCode/GRCh38.p13
    -spikein /home/sheikh/Databases/Drosophila_genome/Drosophila_melanogaster

3. equalRead mode
    greenPipe --modes equalRead --outputdir $(pwd) --inputdir $(pwd)/Fastq --inputfile sampleInfo.txt --libraryType pair
    --reverseName_equalRead True --SelectReads 10000000

9. doughnut mode

    best practice: --sDist should be half of the --mDist

    greenPipe --modes doughnut --dPeakfiles ./Peaks/Jun_merged_narrow-homer.Clean.bed --dNames jun --sDist 200
    --gVer hg38 --sFasta /home/sheikh/Databases/hg38/GenCode/GRCh38.p13.fa --mDist 400 --sProt TGA.TCA,TGA..TCA --mPwm MA0491 --outputdir $(pwd)



#____________ > differential peaks should be in the outputdir + DifferentialPeaks
"""

epilogue = """
Authors:
        (a) Sheikh Nizamuddin:
        snizam001@gmail.com,
        (b) H.Th. Marc Timmers:
        m.timmers@dkfz-heidelberg.de
"""
#_____________________________________________________________________________________________________

parser = argparse.ArgumentParser('greenPipe.py',
                                 usage = __description__,
                                 epilog=colored(epilogue, 'yellow', attrs = ['bold']))

reqNamed = parser.add_argument_group('required arguments\n   _____________________')

reqNamed.add_argument("--modes",
                      help=colored("run mode.", 'green', attrs = ['bold']) + " Multiple modes can be provided as "+
                      "comma seperated values e.g. --modes qc,alignment. " +
                      "Choices are: " +
                      colored("qc, "+
                      "contamination, "+
                      "alignment, "+
                      "equalRead, "+
                      "qcExperiment, "+
                      'initPeakCalling, '+
                      'qcTagDirectories, '+
                      'callPeaks, '+
                      'idr, '+
                      'doughnut, '+
                      'PeakComparison, '+
                      'annotation, '
                      'UserAnnotation, '+
                      'coverageTracks, '+
                      'cutfrequency, '+
                      'initHeatmap, '+
                      'heatmap, '+
#                      'heatmapCoverage, ' +
                      'piggyBack, '+
                      'rnaIntegrate, '+
                      'stats',
                      'green', attrs = ['bold']),
                      type=str,
                      required=True)

reqNamed.add_argument("--outputdir",
                      help="output directory (provide full path)",
                      type=str,
                      required=True)

comNamed = parser.add_argument_group('common arguments\n   _____________________')

comNamed.add_argument("--inputdir",
                      help="input directory having fastq files (provide full path)",
                      default="NA")

comNamed.add_argument("--inputfile",
                      help="input file in .txt format",
                      default="NA")

comNamed.add_argument("--libraryType",
                      help='type of the library',
                      choices=['single','pair'],
                      default="NA")
"""
comNamed.add_argument("--experimentType",
                      help='type of the experiment',
                      choices=['greencutrun','cutrun','cuttag'],
                      default="NA")
"""
comNamed.add_argument("--threads",
                    help="number of threads (default is: total CPU - 2)",
                    type=int,
                    default=(os.cpu_count()-2))

parser.add_argument("--effectiveGenomeSize",
                    help=colored("callPeaks, idr, coverageTracks, initHeatmap mode: ", 'green', attrs = ['bold']) +
                    "effective genome size. "+
                    "human (2913022398: hg38, 2864785220: hg19), "+
                    "mouse (2652783500: GRCm38, 2620345972: GRCm37) and "+
                    "fruitfly (142573017:dm6, 162367812: dm3). "+
                    "(default is 2913022398 which is equivalent for the human genome)",
                    type=int,
                    default=2913022398)

parser.add_argument("--refgenome",
                    help=colored("alignment mode: ", 'green', attrs = ['bold']) +
                    "bowtie2 index file of reference genome. "+
                    "For example: if bowtie2 indexes is present in folder "+
                    "/home/databases/genomes and name of indexes are: "+
                    "GRCh38.p13.1.bt2, GRCh38.p13.2.bt2 etc; then give "+
                    "path /home/databases/genomes/GRCh38.p13",
                    default='None')

parser.add_argument("--spikein",
                    help=colored("alignment mode: ", 'green', attrs = ['bold']) +
                    "bowtie2 index file of spikeIn (Drosophilla). "+
                    "Give path of index as suggested for --referenceGenome",
                    default='None')

parser.add_argument("--alignParam",
                    help=colored("alignment mode: ", 'green', attrs = ['bold']) +
                    "The alignment parameters of the bowtie2/bwa"+
                    ". The default parameter of bowtie2 in this pipeline is: --dovetail --local "+
                    "--very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700 (based on Meers et al. (2019))."+
                    " For single-end default "+
                    "parameters of the program of bwa was used." +
                    "If you want to change it, give parameter as comma seperated values in large bracket e.g. for bowtie2 "+
                    "[--no-unal,--no-mixed,--no-discordant,-I,0,-X,1500] and for bwa see the manual of bwa (for mem)."+
                    " If --gpu is True, then see the options of the nvBowtie: https://nvlabs.github.io/nvbio/nvbowtie_page.html",
                    default='None')

parser.add_argument("--gpu",
                    help="If you have good source of the GPU. Use this option to activate. "+
                    "Wherever neccessary tool automatically recognize and use the GPU source.",
                    default='False',
                    choices=['True','False']
                    )

parser.add_argument("--SelectReads",
                    help=colored("equalRead mode: ", 'green', attrs = ['bold'])+
                    "How many reads should be use randomly? "+
                    "If not given minimum number of read will be "+
                    "choosed using all experimental files",
                    type = int,
                    default = 0
                    )

parser.add_argument("--spikeNormPeak",
                    help=colored("initpeakcalling mode: ", 'green', attrs = ['bold']) +
                    "Should include SpikeIn normalization"+
                    "in peak calling ? (default is True)",
                    type = str,
                    choices=['True','False'],
                    default = 'True'
                    )

parser.add_argument("--reverseName_equalRead",
                    help=colored("equalRead mode: ", 'green', attrs = ['bold']) +
                    'If the bamfile folder already contains *.original.bam files'+
                    ' and you are selecting random reads, then mode equalRead' +
                    ' will denote error. Therefore either manually rename *.original.bam ' +
                    ' to *.bam or use this function to automatically rename files. ',
                    type= str,
                    choices=['True','False'],
                    default = 'False'
                   )

parser.add_argument("--blackListedRegions",
                    help=colored("qcExperiment, qcBamfiles, callPeaks, idr, "+
                      "coverageTracks, initHeatmap, heatmap mode: ", 'green', attrs = ['bold']) +
                    "regions which are black listed by ENCODE"+
                    " in bed format.",
                    type = str,
                    default='None'
                    )

parser.add_argument("--pMethod",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "peak calling algorithm. Default is homer",
                    type=str,
                    choices=['homer',
                             'macs2',
                             'seacr'],
                    default = 'homer'
                   )

parser.add_argument("--pStyle",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "style of the peaks (narrow/broad/both). Default is the narrow. "+
                    "For each samples users can provide comma seperated styles"+
                    "e.g. both,narrow,narrow,broad",
                    type=str,
                    default = 'narrow'
                    )

parser.add_argument("--pFdrHomer",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "FDR/poisson in the homer peak calling. Default 0.001",
                    type=str,
                    default='0.001'
                    )

parser.add_argument("--pPvalueHomer",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "p value in the homer peak calling. Default is 0.0001",
                    type=str,
                    default='0.0001'
                    )

parser.add_argument("--pFcHomer",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Fold changes in the homer peak calling. Default is 4.0",
                    type=str,
                    default='4.0'
                    )

parser.add_argument("--pDistHomer",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Distribution (FDR or poisson)",
                    type=str,
                    choices = ['fdr','poisson'],
                    default='fdr'
                    )

parser.add_argument("--pControl",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Include control in the homer peak calling? Default is True",
                    type=str,
                    choices=['True',
                             'False'],
                    default = 'True'
                   )

parser.add_argument("--pSpike",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Normalize reads by spike-in in the homer peak calling. Default is True",
                    type=str,
                    choices=['True',
                             'False'],
                    default = 'True'
                   )

parser.add_argument("--pSeacrMode",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Mode when calling peaks using SEACR. Default is strigent",
                    type=str,
                    choices=['stringent',
                             'relaxed'],
                    default = 'stringent'
                   )
parser.add_argument("--pSeacrThreshold",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Threshold value if control is not used in peak calling using SEACR"+
                    ". Default is 0.01",
                    type=str,
                    default = '0.01'
                   )

parser.add_argument("--pOpts",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "Pipeline uses all default parameters. If you want to change something, "+
                    "in the peakcalling from homer and macs2 you can give here as comma seperated values in bracket.",
                    type = str,
                    default = 'None'
                   )

parser.add_argument("--genomeFile",
                    help=colored("callPeaks mode: ", 'green', attrs = ['bold']) +
                    "With SEACR mode --genomeFile is require which contains the length of each chromosome."+
                    " These files can be downloaded from  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes."+
                    "Preferebly generate your own genome file from fasta (used in the alignment) by using "+
                    "samtools faidx genome.fa && cut -f1,2 genome.fa.fai > genomeFile.txt."+
                    " By default program uses hg38 genome file (hg38.genome)",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrExprs",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "List of the tagDirectories of the experiments (IDR mode). "+
                    "If --idrMethod is homer then, give as follows: "+
                    " /home/dir1/exp1_rep1,/home/dir1/exp1_rep2;/home/dir1/exp2_rep1,/home/dir1/exp2_rep2,/home/dir1/exp2_rep3."+
                    " Give full path" +
                    ". If --idrMethod is macs2 then give list of the bamfiles",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrCtrl",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "List of the tagDirectories of the control (IDR mode). If --idrMethod is homer then, give as follows: "+
                    " /home/dir1/ctrl1_rep1,/home/dir1/ctrl1_rep2;/home/dir1/ctrl2_rep1,/home/dir1/ctrl2_rep2,/home/dir1/ctrl2_rep3. Give full path"+
                    ". If --idrMethod is macs2 then give list of the bamfiles",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrName",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "List of the Name of the experiment (IDR mode). Give as follows: "+
                    " expr1;expr2",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrExprSpike",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "List of the experiment bamfiles of the spikeIn, if --idrSpike is True. Format is same. Give full path",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrCtrlSpike",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "List of the control bamfiles of the spikeIn, if --idrSpike is True. Format is same. Give full path",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrControl",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "Should include the control in IDR peak calling (IDR mode). Default is True.",
                    type=str,
                    default = 'True',
                    choices = ["True",
                               "False"]
                   )

parser.add_argument("--idrSpike",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "Should include the spikeIn in IDR peak calling (IDR mode). Default is True",
                    type=str,
                    default = 'True',
                    choices = ["True",
                               "False"]
                   )

parser.add_argument("--idrStyle",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "narrow or broad peaks (IDR mode). Default is narrow (factor) or broad (histone)"+
                    "You can give a single style for all experiments e.g. factor or different style "+
                    "for each experiments e.g. factor;histone",
                    type=str,
                    default = 'factor',
                    choices = ["factor",
                               "histone"]
                   )

parser.add_argument("--idrOutput",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "Prefix of the output file for each experiment  (IDR mode) e.g. idr_expr1;idr_Expr2",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--idrMethod",
                    help=colored("idr mode: ", 'green', attrs = ['bold']) +
                    "Peak calling method: homer or macs2 in IDR",
                    type=str,
                    default = 'homer',
                    choices =['homer','macs2']
                   )

parser.add_argument("--overFiles",
                    help=colored("PeakComparison mode: ", 'green', attrs = ['bold']) +
                    "Input peak file for the comparison. Give as CSV input e.g. dir1/dir2/peak1.bed,dir1/dir2/peak2.bed"+
                    ". If not give, pipeline automatically identify the IDR peaks or peaks from the output folder"+
                    ". Give the full path of file.",
                    type=str,
                    default = 'NA'
                   )

parser.add_argument("--overDist",
                    help=colored("PeakComparison mode: ", 'green', attrs = ['bold']) +
                    "The distance between peaks to be consider as overlapping. "+
                    "Default is 400 base-pairs",
                    type=str,
                    default = '400'
                   )

parser.add_argument("--annpeakFiles",
                    help=colored("UserAnnotation/annotation mode: ", 'green', attrs = ['bold']) +
                    "in this mode you can provide peaks which can be annotated "+
                    "by user provided annotations files (--annFiles). "+
                    "Multiple peak files can be given as comma seperated values.",
                    type = str,
                    default='None')

parser.add_argument("--annFiles",
                    help=colored("UserAnnotation mode: ", 'green', attrs = ['bold']) +
                    "provide bed files of your region of interest "+
                    "e.g. h3k4me3, jun/fos motif location etc. Users can provide more than one" +
                    " annotation bed files as comma seperated values e.g. f1.bed,f2.bed",
                    type = str,
                    default='None')

parser.add_argument("--annSize",
                    help=colored("UserAnnotation mode: ", 'green', attrs = ['bold']) +
                    "maximum distance between the boundary of peaks (--annpeakFiles) "+
                    "and annotation bed files (annFiles). For each annotation bed files"+
                    " provide this maximum distance e.g. for h3k4me3 500, for jun/fos 100 "+
                    "as 500,100",
                    type = str,
                    default='None')

parser.add_argument("--annName",
                    help=colored("UserAnnotation mode: ", 'green', attrs = ['bold']) +
                    "give name of your annotations "+
                    "e.g. give h3k4me3,jun/fos for file f1.bed,f2.bed",
                    type = str,
                    default='None')

parser.add_argument("--annPrefix",
                    help=colored("UserAnnotation/annotation mode: ", 'green', attrs = ['bold']) +
                    "prefix of the output annotated files. Use this option if peak information "+
                    "is given by users by --annpeakFiles",
                    type = str,
                    default='None')

parser.add_argument("--covSpike",
                    help=colored("coverageTracks mode: ", 'green', attrs = ['bold']) +
                    "During calculation of the whole genome coverage in bins, reads should be "+
                    "normalize according to spikeIn. This normalization will be calculated within the "+
                    "samples given in the --inputfile. Default is False",
                    type = str,
                    choices=['True',
                             'False'],
                    default = 'False'
                   )

parser.add_argument("--covSpike_NormalizationFormula",
                    help=colored("coverageTracks mode: ", 'green', attrs = ['bold']) +
                    "To normalize coverage tracks with spikeIn, choose formula among these. "+
                    "Suppose, x1 ... xn is spikeIn reads per human reads in "+
                    "sample s1 ... sn. Then, per human read spikeIn reads will be: x1/s1"+
                    ", ... xn/sn. Suppose this ratio is r1 .. rn, then scaling factor is (1)"+
                    " min(r1...rn)/r1... min(r1..rn)/rn and (2) 0.05/r1 ... 0.05/rn."+
                    " Default is 1. Notice that if you use method 1, you can compare coverage within"+
                    " given set of samples only. If you want to compare new samples, then you have to "
                    "generate coverage file again.",
                    type = int,
                    choices=['1',
                             '2'],
                    default = '1'
                   )

parser.add_argument("--covOtherOptions",
                    help=colored("coverageTracks mode: ", 'green', attrs = ['bold']) +
                    "Pipeline uses all default parameters. If you want to change something, "+
                    "besides --bl, --effectiveGenomeSize and -p, the other options of "+
                    " bamCoverage (deepTools) can be provided here as comma seperated values but in bracket."+
                    " ",
                    type = str,
                    default = 'None'
                   )

parser.add_argument("--cMotif",
                    help=colored("cutfrequency mode: ", 'green', attrs = ['bold']) +
                    "cutfrequency mode needs location of the motifs in the whole genome "+
                    "If you do not have this file then create it by using --cMotif True ",
                    type = str,
                    choices=['True',
                             'False'],
                    default = 'False'
                   )

parser.add_argument("--cMotifFile",
                    help=colored("cutfrequency mode: ", 'green', attrs = ['bold']) +
                    "If the --cMotif is False, then give the file which contains the location of the "+
                    "motifs in the whole genome. It should be in HOMER bed file format",
                    type = str,
                    default = 'None'
                   )

parser.add_argument("--cMaN",
                    help=colored("cutfrequency mode: ", 'green', attrs = ['bold']) +
                    "Provide the MA number of the JASPER motif "+
                    "Go to http://jaspar.genereg.net/ for finding this number e.g. for TP53 "+
                    " number is MA0106",
                    type = str,
                    default="None"
                   )

parser.add_argument("--cGVersion",
                    help=colored("cutfrequency/annotation: ", 'green', attrs = ['bold']) +
                    "In the cutfrequency mode, if --cMotif is true, provide the version of the genome, "+
                    "so that location of the motifs can be generated. This one is also require in "+
                    "the annotation mode "+
                    " if you are using other species than human. Specify here. For mouse use mm10.",
                    type = str,
                    default = 'hg38'
                   )

parser.add_argument("--cCenter",
                    help=colored("cutfrequency mode: ", 'green', attrs = ['bold']) +
                    "Give the center of the motif",
                    type = str,
                    default="None"
                   )

parser.add_argument("--initHeatmapOtherOptions",
                    help=colored("initHeatmap mode: ", 'green', attrs = ['bold']) +
                    "Pipeline uses all default parameters. If you want to change something, "+
                    "besides --bl, --effectiveGenomeSize and -p, the other options of "+
                    " bamCompare (deepTools) can be provided here as comma seperated values."+
                    " ",
                    type = str,
                    default = 'None'
                   )

parser.add_argument("--hMOpt",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Pipeline uses all default parameters. If you want to change something, "+
                    "besides --missingDataAsZero,-bl,--smartLabels,-p,--metagene,--samplesLabel, the other options of "+
                    " computeMatrix (deepTools) can be provided here as comma seperated values."+
                    "(moreover, by default -a 5000, -b 5000 are used. It can be changed here.)",
                    type = str,
                    default = 'None'
                   )

parser.add_argument("--hPOpt",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Pipeline uses all default parameters. If you want to change something, "+
                    "besides --refPointLabel and --dpi, the other options of "+
                    " plotHeatmap (deepTools) can be provided here as comma seperated values."+
                    "(moreover, by default--colorMap GnBu are used. It can be changed here.)",
                    type = str,
                    default = 'None'
                   )

parser.add_argument("--hSpike",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Should include spike in preparation of bamcompare files. "+
                    "Default is True",
                    type = str,
                    choices=['True',
                             'False'],
                    default = 'True'
                   )

parser.add_argument("--hCovComp",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Heatmap should be plotted for bamcoverage (spikeIn normalized) or "+
                    "bamcompare files? Leave it, if you want to given your own files using --hInFiles",
                    type = str,
                    default = "NA",
                    choices=['compare',
                             'coverage',
                            "NA"]
                   )

parser.add_argument("--hCovMethod",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "If --hCovComp is bamcoverage, then either you (1) can use spikeIn normalize bamcoverage file or "+
                    " (2) simply a bamcoverage file which will be normalized during heatmap production using following"+
                    " formula: sum of reads within the bins / (spikeIn reads/10,000). Default is 1. But I recommend to use 2."+
                    " Leave it, if you want to given your own files using --hInFiles",
                    type = int,
                    default = 1,
                    choices=[1,
                             2,
                             "NA"]
                   )

parser.add_argument("--hInCounts",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "If you want to given your own files using --hInFiles and want to normalize it also, then provide "+
                    "normalization count for each samples in comma seperated value format. It is basically spikeIn-in-sample-i/10,000."+
                    " and use --hCovMethod 2.",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--hInFiles",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Input bamcoverage or bamcompare files in csv format" +
                    ". If file is not given, tool will indentify files from the "+
                    "--inputfile",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--hInNames",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Name of the samples in csv format. Equal to the number of files in --hInFiles.",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--hRegionMode",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "mode of the region? on tss, on given peak file or select automatically total peak file"+
                    " or differential peaks from Peak folder or metagene. Default is peaks mode.",
                    type = str,
                    default = "peaks",
                    choices=['tss',
                             'metagene',
                             'bed',
                             "peaks"]
                   )

parser.add_argument("--hGtf",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "Path of the GTF file, if --hRegionMode is tss or metagene",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--hBed",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "List of bed files on which you want to generate heatmap in csv format, if --hRegionMode bed",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--hDiffPeaks",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "In the --hRegionMode peaks, should program consider differential peaks also ?"+
                    " Default is the True",
                    type = str,
                    choices=['True',
                             'False'],
                    default = 'True'
                   )
"""
parser.add_argument("--hCoverage",
                    help=colored("heatmap mode: ", 'green', attrs = ['bold']) +
                    "If you want to generate heatmap coverage plot. Default is False.",
                    type = str,
                    default = 'False',
                    choices=['True',
                             'False'],
                   )
"""
parser.add_argument("--lProt",
                    help=colored("piggyBack mode: ", 'green', attrs = ['bold']) +
                    "Files with name of the proteins per line."+
                    " These are those proteins which are significantly enriched in the mass spectrometry experiment "+
                    " and you want to check if these proteins are providing any piggy-back binding to your "+
                    "protein of interest. You can find the list of proteins on basis of the stochiometry "+
                    "(iBaq values) observed in your MassSpectrometry experiment.  "+
                    "This mode has limitation because it can predict the piggy-back binding events only for those proteins "+
                    "which are binding to the DNA (not histone modifications proteins) and their PWM are known. "+
                    "This list should not contain any special character e.g. ; or - or :",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--mPwm",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "Provide the MA number of the JASPER motif "+
                    "Go to http://jaspar.genereg.net/ for finding this number e.g. for TP53 "+
                    " number is MA0106"+
                    ". You can provide more than 1 motifs as comma seperated values. ",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--sProt",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "provide the sequence of protein of interest for which you are searching "+
                    "piggy-back binding event e.g. for NFYA you can give CCAAT, for JUN you can give "+
                    "TGANTCA, for ATF7 you can give TGANNTCA. You can provide more than one sequence using "+
                    "comma seperated values. In the sequence IUPAC codes are also applicable. ",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--mPeak",
                    help=colored("piggyBack: ", 'green', attrs = ['bold']) +
                    "Path of the peak file of greenPipe experiment in bed format. "+
                    "It is optional. If not given, tool will automatically search peak file on basis of "+
                    "your --outputdir and --inputfile. If provided give the full path of peaks as comma "+
                    "seperated values e.g. /home/xyz/exp1.Clean.bed,/home/xyz/exp2.Clean.bed. ",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--mDist",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "Search motifs from this distance from the centre of peaks. Default is 400",
                    default=400,
                    type=int
                   )

parser.add_argument("--distPiggy",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "Distance to find piggy back binding events. Default is 400",
                    default=400,
                    type=int
                   )

parser.add_argument("--Species",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "Taxon ID or name of Species according to NDBI e.g. "+
                    "for home sapiens taxon ID is 9606 and name is human",
                    default="NA",
                    type=str,
                    choices=['human',
                             "mouse",
                             "rat",
                             "fruitfly",
                             "nematode",
                             "zebrafish",
                             "thale-cress",
                             "frog",
                             "pig"]
                   )

parser.add_argument("--sFasta",
                    help=colored("piggyBack/doughnut/annotation mode: ", 'green', attrs = ['bold']) +
                    "Fasta file of the species genome",
                    default="NA",
                    type=str
                   )

parser.add_argument("--gVer",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "Genome version of the fasta file. Default is hg38",
                    default="hg38",
                    type=str
                   )

parser.add_argument("--mPvalue",
                    help=colored("piggyBack mode: ", 'green', attrs = ['bold']) +
                    "Motif finding pvalue for piggy back binding events. Default is 0.05",
                    default=0.05,
                    type=int
                   )

parser.add_argument("--sDist",
                    help=colored("piggyBack and doughnut mode: ", 'green', attrs = ['bold']) +
                    "Search sequence of the motif from the center of peak. Default is 200."+
                    " best practice: --sDist should be half of the --mDist",
                    default=200,
                    type=int
                   )

parser.add_argument("--mPrefix",
                    help=colored("piggyBack mode: ", 'green', attrs = ['bold']) +
                    "If giving your own bed peak files using --mPeak then give prefix of output"+
                    "for each bedfiles",
                    default="NA",
                    type = str
                   )

parser.add_argument("--dPeakfiles",
                    help=colored("doughnut mode: ", 'green', attrs = ['bold']) +
                    "Path of the peak file for doughnut mode in bed format as comma "+
                    "seperated values e.g. /home/xyz/exp1.Clean.bed,/home/xyz/exp2.Clean.bed. ",
                    type = str,
                    default = "NA"
                   )

parser.add_argument("--dNames",
                    help=colored("doughnut mode: ", 'green', attrs = ['bold']) +
                    "name of the experiment for each bed file given with --dPeakfiles "+
                    "e.g. exp1,exp2. ",
                    type = str,
                    default = "NA"
                   )

args = parser.parse_args()
#------
outputdir=args.outputdir
inputdir=args.inputdir
inputfile=args.inputfile
modes=args.modes
libraryType=args.libraryType
#experimentType=args.experimentType
#------
threads=args.threads
effectiveGenomeSize=args.effectiveGenomeSize
gpu=args.gpu
#------
spikein=args.spikein
refgenome=args.refgenome
alignParam=args.alignParam
#------
SelectReads=args.SelectReads
#------
spikeNormPeak=args.spikeNormPeak
#------
reverseName_equalRead=args.reverseName_equalRead
#------
blackListedRegions=args.blackListedRegions
#------
annpeakFiles=args.annpeakFiles
annFiles=args.annFiles
annSize=args.annSize
annName=args.annName
annPrefix=args.annPrefix
#------
pMethod=args.pMethod
pStyle=args.pStyle
pFdrHomer=args.pFdrHomer
pPvalueHomer=args.pPvalueHomer
pFcHomer=args.pFcHomer
pDistHomer=args.pDistHomer
pControl=args.pControl
pSpike=args.pSpike
pSeacrMode=args.pSeacrMode
pSeacrThreshold=args.pSeacrThreshold
genomeFile=args.genomeFile
pOpts=args.pOpts
#------
idrExprs=args.idrExprs
idrCtrl=args.idrCtrl
idrName=args.idrName
idrControl=args.idrControl
idrStyle=args.idrStyle
idrOutput=args.idrOutput
idrMethod=args.idrMethod
idrExprSpike=args.idrExprSpike
idrCtrlSpike=args.idrCtrlSpike
idrSpike=args.idrSpike
#-------
overDist=args.overDist
overFiles=args.overFiles
#-------
covSpike=args.covSpike
covOtherOptions=args.covOtherOptions
covSpike_NormalizationFormula=args.covSpike_NormalizationFormula
#-------
cMotif=args.cMotif
cMaN=args.cMaN
cGVersion=args.cGVersion
cCenter=args.cCenter
cMotifFile=args.cMotifFile
#--------
hSpike=args.hSpike
initHeatmapOtherOptions=args.initHeatmapOtherOptions
#--------
hInFiles=args.hInFiles
hInNames=args.hInNames
hRegionMode=args.hRegionMode
hGtf=args.hGtf
hBed=args.hBed
hCovComp=args.hCovComp
hCovMethod=args.hCovMethod
hDiffPeaks=args.hDiffPeaks
hMOpt=args.hMOpt
hPOpt=args.hPOpt
hInCounts=args.hInCounts
#hCoverage=args.hCoverage
#----
dPeakfiles=args.dPeakfiles
dNames = args.dNames
#----
lProt=args.lProt
mPwm=args.mPwm
sProt=args.sProt
mPeak=args.mPeak
mDist=args.mDist
distPiggy=args.distPiggy
Species=args.Species
sFasta=args.sFasta
gVer=args.gVer
sDist=args.sDist
mPvalue=args.mPvalue
mPrefix=args.mPrefix
#-- logfile
#_____________________________________________________________________________________________________
with open(outputdir+'/'+'log.txt','w') as logfile:
    logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

#-- few definitions
#_____________________________________________________________________________________________________

def checkBlackList (blackListedRegions):
    if blackListedRegions == 'None':
        print(colored("Provide blacklisted regions of ENCODE using "+
                      "--blackListedRegions. ",
                      'red',
                      attrs = ['bold']
                     )
             )
        exit()
    elif not os.path.exists(blackListedRegions):
        print(colored("blackListedRegions file does not exist.",
                      'red',
                      attrs = ['bold']
                     )
             )
        exit()
    elif os.stat(blackListedRegions).st_size == 0:
        print(colored("blackListedRegions file is empty.",
                      'red',
                      attrs = ['bold']
                     )
             )
        exit()

#--- check python packages installed or not
#_____________________________________________________________________________________________________
mypackages=['pkg_resources',
            'sys',
            'subprocess',
            'argparse',
            'os',
            'termcolor',
            're',
            'pandas',
            'importlib',
            'glob',
            'itertools',
            'datetime',
            'mygene',
            'Bio',
            'pyfaidx',
            'upsetplot',
            'more_itertools',
            'pathlib',
            'gzip'
            ]

print(colored('-> Checking if python packages are installed or not',
              'green',
              attrs=['bold']))

print(colored('-> Install IDR always from conda, if it is not installed.',
              'green',
              attrs=['bold']))

for mypackage in mypackages:
    if mypackage in sys.modules:
        print(mypackage + "is installed and already present in sys.modules")
    else:
        print("can't find the " +
              mypackage +
              "module, trying to install it")
        c=['conda',
           'install',
           mypackage]

        universal.run_cmd(c,outputdir)

        c=['conda',
           'install',
           '-c',
           'bioconda',
           mypackage ]

        universal.run_cmd(c,outputdir)

        c=['conda',
        'install',
        '-c',
        'conda-forge',
        mypackage]

        universal.run_cmd(c,outputdir)

#-- checking if other neccessary packages installed or not
#_____________________________________________________________________________________________________
#  trimgalore, Fastqc, bowtie2, samtools, Homer, bedtools, deeptools, Real_differentialPeaks.R present or not?
myothercommands=['makeTagDirectory',
                 'bedtools',
                 'bamCoverage',
                 'samtools',
                 'trim_galore',
                 'fastqc',
                 'bowtie2',
                 'fastq_screen',
                 'pv',
                 'pigz'
                ]

print(colored('-> Checking if other packages are installed or not',
              'green',
              attrs=['bold']))

for myothercommand in myothercommands:
    if myothercommand == 'makeTagDirectory':
        mytool='homer'
    elif myothercommand == 'bamCoverage':
        mytool='deeptools'
    else:
        mytool=myothercommand
    if mytool!='fastqc':
        try:
            subprocess.call([myothercommand],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            print(mytool + ' packages in installed in the system')
        except FileNotFoundError:
            print(colored(myothercommand + ": packages is absent in your system. "+
                          "Please install this tools before running script.",
                          'green',
                          attrs=['bold']))
            exit()
    else:
        try:
            subprocess.call([myothercommand,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            print(mytool + ' packages in installed in the system')
        except FileNotFoundError:
            print(colored(myothercommand + ': packages is absent in your system. '+
                          "Please install this tools before running script" +
                          " If it is fastq_screen, see the documentation: "+
                          "https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_documentation.html." +
                          " Install the related database of the fastq_screen also "+
                          "for checking the contamination using following command: fastq_screen --get_genomes",
                          'green',
                          attrs=['bold']))
            exit()
    #-----
try:
    subprocess.call(['idr','--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    print('idr' + ' packages in installed in the system')
except FileNotFoundError:
    print(colored('idr : packages is absent in your system. '+
                  "Installing it by using following command:"+
                  "conda install -c bioconda idr",
                  'green',
                  attrs=['bold']))
    c=['conda',
    'install',
    '-c',
    'bioconda',
    'idr' ]

    universal.run_cmd(c,outputdir)

try:
    subprocess.call(['featureCounts'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    print(mytool + ' packages in installed in the system')
except FileNotFoundError:
    print(colored(myothercommand + ": packages is absent in your system. "+
                  "Installing it.",
                  'green',
                  attrs=['bold']))
    c=['conda',
    'install',
    '-c',
    'bioconda',
    'subread' ]

    universal.run_cmd(c,outputdir)

#-- all about required arguments
#_______________________________________________________________________________

def main ():
    for mode in modes.split(','):
        if mode != 'idr' and mode != 'doughnut' and (mode == 'PeakComparison' and overFiles == 'NA') and mode != "UserAnnotation":
            if inputdir=="NA" or inputfile == "NA" or libraryType=="NA": #or experimentType == "NA"
                print(colored("Options --inputdir, --inputfile, --libraryType"+ #or --experimentType
                              " is missing",
                              "green",
                              attrs = ["bold"]
                             )
                     )

        if not os.path.exists(outputdir):
            os.makedirs(outputdir)

        if mode!='idr' and mode != 'doughnut' and (mode == 'PeakComparison' and overFiles == 'NA'):

            if not os.path.exists(inputdir):
                print(colored('Error: input directory does not exist.',
                              'green', attrs=['bold']))
                exit()

            if not os.path.exists(inputfile):
                print(colored('inputfile does not exist.',
                             'green',attrs=['bold']))
                exit()
            else:
                myin=pd.read_csv(inputfile,header=None,sep='\t')
                columnF = 5
                if myin.shape[1] < columnF:
                    print(myin)
                    print(colored('Error: inputfile is not right. '+
                                  'Number of the columns are not equal to ' + columnF,
                                  'green', attrs=['bold']))
                    exit()
                elif myin.empty:
                    print(colored('Error: inputfile is empty',
                                  'green', attrs=['bold']))
                    exit()
                else:
                    if libraryType=='single':
                        I=[0,2]
                    elif libraryType=='pair':
                        I=[0,1,2,3]
                    notPresent=0;
                    for i in I:
                        for j in range(0,myin.shape[0]):
                            if not os.path.exists(inputdir + '/' + myin.iloc[j,i]):
                                notPresent=notPresent+1
                                print(colored('Error: file does not exist: '+ myin.iloc[j,i],
                                              'green', attrs=['bold']))
                    if notPresent > 0:
                        exit()

        #-- link the fastq files
        #_____________________________________________________________________________________________________
        if mode != 'idr' and mode != 'doughnut':
            if not os.path.exists(outputdir+'/Fastq/'):
                os.makedirs(outputdir+'/Fastq/')
#           if experimentType == "greencutrun":
            myin=pd.read_csv(inputfile,header=None,sep='\t')

            myin_Files=sum(myin.iloc[:,[0, 1, 2, 3]].values.tolist(), [])
            myin_Files=[item for item in myin_Files if str(item) != 'nan']

            print(colored('Checking format of input fastq files','green', attrs=['bold']))

            for myin_File in myin_Files:
                myin_FileType = filetype.guess(inputdir+'/'+myin_File)
                if myin_FileType is None:
                    print(colored('Error:'+inputdir+'/'+myin_File+' inputfile '+
                            '(FASTQ) should be in .gz format'+
                            '. Compress your file',
                            'green', attrs=['bold']))
                else:
                    print(colored(myin_File + ": " + str(myin_FileType), 'green', attrs = ['bold']))

            lf.linksFile(libraryType,
                         myin,
                         inputdir,
                         outputdir)

        #-- qc pass
        #______________________________________________________________________________________________________
        if mode == 'qc':
            total_cmd=[]
            total_fastqc=[]
#            if experimentType == "greencutrun":
            for i in range(0,myin.shape[0]):
                Name = myin.loc[i,4]
                x, y = qcFastQ.qc(libraryType, Name, inputdir, outputdir, threads)
                for j in x:
                    total_cmd.append(j)
                for j in y:
                    total_fastqc.append(j)
            print(total_cmd)
            with Pool (threads) as p:
                p.starmap(universal.run_cmd,
                          zip(total_cmd,
                              repeat(outputdir)
                             )
                         )
            print(total_fastqc)
            with Pool (threads) as p:
                p.starmap(universal.run_cmd,
                          zip(total_fastqc,
                              repeat(outputdir)
                             )
                         )

        #-- contamination
        #_____________________________________________________________________________________________________
        if mode == 'contamination':
          total_cmd=[]
#         if experimentType == "greencutrun":
          for i in range(0,myin.shape[0]):
              Name = myin.loc[i,4]
              x = contamination.contamination (libraryType, Name, inputdir, outputdir, threads)
              for j in x:
                  total_cmd.append(j)
          print(total_cmd)
          with Pool (threads) as p:
              p.starmap(universal.run_cmd,
                        zip(total_cmd,
                            repeat(outputdir)
                           )
                       )

        #-- alignment
        #_____________________________________________________________________________________________________
        if mode == 'alignment':

            if refgenome=='None' or spikein=='None':

                print(colored("--refgenome and --spikein is missing in the command. "+
                              "These are necessary files for the alignment mode",
                              'green', attrs=['bold']))
                exit()

#            if experimentType == "greencutrun":
            for i in range(0,myin.shape[0]):
                Name = myin.loc[i,4]
                align.alignment(libraryType, outputdir, Name, refgenome, spikein, threads, alignParam, gpu)

        #-- selecting equal number of the reads
        if mode == 'equalRead':
#            if experimentType == "greencutrun":
          Names=[]
          for i in range(0,myin.shape[0]):
              Names.append(myin.loc[i,4])

          if reverseName_equalRead == 'True':
              equalRead.reverseName (Names,
                                     outputdir,
                                     threads)
          equalRead.equalRead(libraryType,
                              outputdir,
                              Names,
                              SelectReads,
                              threads)

        #-- qc  of the bamfiles of samples
        #_____________________________________________________________________________________________________
        if mode == "qcExperiment":
            if not os.path.exists(outputdir+'/Bamfiles_QC/'):
              os.makedirs(outputdir+'/Bamfiles_QC/')

            checkBlackList(blackListedRegions)

            Names=[]
            FlagstatFiles=[]
            fastQs=[]
            qcBfiles=[]

            for i in range(0,myin.shape[0]):
                Name = myin.loc[i,4]
                Names.append(Name)
                FlagstatFiles.append(outputdir + "/Bamfiles/" + Name + "_expr.Flagstats.txt")
                FlagstatFiles.append(outputdir + "/SpikeIn/" + Name + "_expr.Flagstats.txt")
                FlagstatFiles.append(outputdir + "/Bamfiles/" + Name + "_control.Flagstats.txt")
                FlagstatFiles.append(outputdir + "/SpikeIn/" + Name + "_control.Flagstats.txt")

                qcBfiles.append(outputdir + "/Bamfiles/" + Name + "_expr.bam")

                if libraryType == "pair":
                    fastQs.append (outputdir + "/Fastq/" + Name + "_expr_R1.fastq.gz")
                    fastQs.append (outputdir + "/Fastq/" + Name + "_control_R1.fastq.gz")
                    fastQs.append (outputdir + "/Trim_galore/" + Name + "_expr_R1_val_1.fq.gz")
                    fastQs.append (outputdir + "/Trim_galore/" + Name + "_control_R1_val_1.fq.gz")

                elif libraryType == "single":
                    fastQs.append (outputdir + "/Fastq/" + Name + "_expr.fastq.gz")
                    fastQs.append (outputdir + "/Fastq/" + Name + "_control.fastq.gz")
                    fastQs.append (outputdir + "/Trim_galore/" + Name + "_expr_trimmed.fq.gz")
                    fastQs.append (outputdir + "/Trim_galore/" + Name + "_control_trimmed.fq.gz")

            qcBfiles = ','.join(str(e) for e in qcBfiles)
            qcBamfiles.countReads (Names, fastQs, FlagstatFiles, outputdir, libraryType)
            qcBamfiles.qc_bam (qcBfiles, outputdir, str(threads), blackListedRegions, Names)

        #-- initiate peak calling or generate tag directories
        #_____________________________________________________________________________________________________
        if mode == "initPeakCalling":
#            if experimentType == "greencutrun":
          Names = []
          for i in range(0,myin.shape[0]):
              Name = myin.loc[i,4]
              Names.append(Name)
              totalCmd=initPeakCalling.initPeakCalling(libraryType,
                                                       outputdir,
                                                       Names,
                                                       spikeNormPeak,
                                                       threads)

          with Pool (threads) as p:
              p.starmap(universal.run_cmd,
                    zip(totalCmd,
                        repeat(outputdir)
                       )
                       )
        #-- qc of the tagdirectories
        #_____________________________________________________________________________________________________
        if mode=="qcTagDirectories":

            qcTagR=psource.resource_filename(__name__, "rscripts/qcTag.R")

#            if experimentType == "greencutrun":
            Names = []
            for i in range(0,myin.shape[0]):
                Name = myin.loc[i,4]
                Names.append(Name)
                totalCmd=qcTagD.qcTagD(libraryType,
                              outputdir,
                              Names,
                              threads,
                              qcTagR)
            with Pool (threads) as p:
                p.starmap(universal.run_cmd,
                          zip(totalCmd,
                              repeat(outputdir)
                             )
                         )
        #-- calling peaks
        #_____________________________________________________________________________________________________
        if mode == "callPeaks":

            checkBlackList(blackListedRegions)

#            if experimentType == "greencutrun":
            Names = []
            for i in range(0,myin.shape[0]):
                Name = myin.loc[i,4]
                Names.append(Name)

            if pMethod == "homer":
                callPeak.callPeaksHomer (outputdir,
                                         Names,
                                         threads,
                                         pControl,
                                         pStyle,
                                         blackListedRegions,
                                         pFdrHomer,
                                         pPvalueHomer,
                                         pFcHomer,
                                         pDistHomer,
                                         pOpts
                                        )

            elif pMethod == "macs2":
                callPeak.callPeaksMacs2 (outputdir,
                                         Names,
                                         threads,
                                         pSpike,
                                         pControl,
                                         pStyle,
                                         effectiveGenomeSize,
                                         blackListedRegions,
                                         libraryType,
                                         pOpts
                                        )

            elif pMethod == "seacr":
                if genomeFile == "NA":
                    genome = psource.resource_filename(__name__, "data/hg38.genome")
                seacr=psource.resource_filename(__name__, "otherScripts/SEACR/SEACR_1.3.sh")
                callPeak.callPeaksSEACR (outputdir,
                                         Names,
                                         threads,
                                         pSpike,
                                         pControl,
                                         seacr,
                                         libraryType,
                                         pSeacrMode,
                                         pSeacrThreshold,
                                         blackListedRegions,
                                         genome
                                        )

        #-- IDR peak idr_macs2 ():
        #_____________________________________________________________________________________________________

        if mode == 'idr':

            checkBlackList(blackListedRegions)

            if idrControl == True and len(idrExprs.split(';')) != len(idrCtrls.split(';')):
                print(colored("number of the experiments and control are not same. Check --idrExprs and --idrCtrls",
                              "green",
                              attrs = ["bold"]
                             )
                     )

            if len(idrExprs.split(';')) != len(idrName.split(';')):
                print(colored("number of the experiment and names are not equal. Check --idrExprs and --idrName",
                              "green",
                              attrs = ["bold"]
                             )
                     )

            if len(idrStyle) == 1:
                idrStyles=[idrStyle]*len(idrExprs.split(';'))
            else:
                idrStyles=idrStyle.split(';')

            if idrMethod == 'homer':
                for i_idr in range(0,len(idrExprs.split(';'))):
                    idr.idr_homer (outputdir,
                             idrExprs.split(';')[i_idr],
                             idrCtrl.split(';')[i_idr],
                             idrName.split(';')[i_idr],
                             idrControl,
                             idrStyles[i_idr],
                             idrOutput.split(';')[i_idr],
                             blackListedRegions,
                             threads
                            )
            elif idrMethod == "macs2":
                if libraryType == "NA":
                    print(colored("Please define libraryType",
                                  "green",
                                  attrs = ["bold"]
                                 )
                         )
                    exit()
                print (colored("calling idr peaks using macs2 with libraryType" + libraryType,
                               "green",
                               attrs = ["bold"]
                              )
                      )
                for i_idr in range(0,len(idrExprs.split(';'))):
                    idr.idr_macs2(outputdir,
                                  idrSpike,
                                  idrExprs.split(';')[i_idr],
                                  idrCtrl.split(';')[i_idr],
                                  idrExprSpike.split(';')[i_idr],
                                  idrCtrlSpike.split(';')[i_idr],
                                  idrName.split(';')[i_idr],
                                  idrControl,
                                  idrStyles[i_idr],
                                  idrOutput.split(';')[i_idr],
                                  blackListedRegions,
                                  threads,
                                  effectiveGenomeSize,
                                  libraryType)


        #-- doughnut mode
        #______________________________________________________________________________________________________

        if mode == 'doughnut':
            dirs=[outputdir+'/'+'doughnut/']
            for d in dirs:
                if not os.path.exists(d):
                    os.makedirs(d)


            inFiles = dPeakfiles.split(',')

            if len(inFiles) != len(dNames.split(',')):
                print(colored("lenght of the names is not equal to the peak/bed files given in doughnut mode",
                              "green",
                              attrs = ["bold"]
                             )
                     )

            inFile_exist = 0
            for inFile in inFiles:
                if not os.path.exists(inFile):
                    print(colored(inFile + ": does not exist",
                                  'green',
                                  attrs = ["bold"]
                                 )
                         )
                    inFile_exist += 1

            if inFile_exist > 0:
                exit()

            x = -1
            for dName in dNames.split(','):
                o1 = outputdir+'/'+'doughnut/' + dName + '.motifs_pwm.txt'
                o2 = outputdir+'/'+'doughnut/' + dName + '.peaks_WithoutPwmSeq.txt'
                o3 = outputdir+'/'+'doughnut/' + dName + '.peaks_WithoutPwmSeq'
                tempDir = outputdir+'/'+'doughnut/'
                x+=1
                MassSpectro_Fasta=pyfaidx.Fasta(sFasta)
                filterMotifSeq.filterMotifSequence(inFiles[x],o1,o2,o3,mPwm,sProt,tempDir,gVer,mDist,threads,sDist,MassSpectro_Fasta,outputdir)


        #-- differential peaks or mode PeakComparison
        #_____________________________________________________________________________________________________

        if mode == 'PeakComparison':

            #--- finding the overlapping peaks
            if overFiles == "NA":

                if os.path.exists(outputdir+'/idr_homer/'):

                    if not os.path.exists(outputdir+'/OverlappingPeaks/idr_homer'):
                        os.makedirs(outputdir+'/OverlappingPeaks/idr_homer')

                    overFiless=glob.glob(
                        outputdir+'/idr_homer/*/TotalExperiment.peaks-top-set.Clean.bed'
                        )

                    overInFiles = []

                    if len(overFiless) > 1:
                        for overFile in overfiless:
                            overName = overFile.split('/')[-2]
                            os.symlink(overFile, outputdir + '/OverlappingPeaks/idr_homer/' + overName )
                            overInFiles.append(outputdir + '/OverlappingPeaks/idr_homer/' + overName )

                    comparePeak.overlapPeaks(overInFiles,
                        outputdir+'/OverlappingPeaks/idr_homer',
                        overDist,
                        outputdir,
                        outputdir+'/OverlappingPeaks/idr_homer/' )

                elif os.path.exists(outputdir+'/idr_macs2/'):

                    if not os.path.exists(outputdir+'/OverlappingPeaks/idr_macs2'):

                        os.makedirs(outputdir+'/OverlappingPeaks/idr_macs2')

                    overFiless=glob.glob(
                        outputdir+'/idr_macs2/*/TotalExperiment.peaks-top-set.Clean.bed'
                        )

                    overInFiles = []

                    if len(overFiless) > 1:
                        for overFile in overfiless:
                            overName = overFile.split('/')[-2]
                            os.symlink(overFile, outputdir + '/OverlappingPeaks/idr_macs2/' + overName )
                            overInFiles.append(outputdir + '/OverlappingPeaks/idr_macs2/' + overName )

                    comparePeak.overlapPeaks(overInFiles,
                        outputdir+'/OverlappingPeaks/idr_macs2',
                        overDist,
                        outputdir,
                        outputdir+'/OverlappingPeaks/idr_macs2/'
                        )

                elif os.path.exists(outputdir+'/Peaks/'):

                    if not os.path.exists(outputdir+'/OverlappingPeaks/Peaks'):
                        os.makedirs(outputdir+'/OverlappingPeaks/Peaks')

                    overFiless=glob.glob(
                        outputdir+'/Peaks/*.Clean.bed'
                        )

                    overInFiles = []

                    if len(overFiless) > 1:
                        for overFile in overFiless:
                            overName = overFile.split('/')[-1].replace("_all-homer.Clean.bed","").replace("_broad-homer.removed.bed","").replace("_narrow-homer.removed.bed","")
                            os.symlink(overFile, outputdir + '/OverlappingPeaks/Peaks/' + overName )
                            overInFiles.append(outputdir + '/OverlappingPeaks/Peaks/' + overName )

                    comparePeak.overlapPeaks(overInFiles,
                        outputdir+'/OverlappingPeaks/Peaks',
                        overDist,
                        outputdir,
                        outputdir+'/OverlappingPeaks/Peaks/'
                        )

            else:
                overInFiles = overFiles.split(',')
                if not os.path.exists(outputdir+'/OverlappingPeaks/UserGiven'):
                    os.makedirs(outputdir+'/OverlappingPeaks/UserGiven')
                    comparePeak.overlapPeaks(overInFiles,
                        outputdir+'/OverlappingPeaks/UserGiven',
                        overDist,
                        outputdir,
                        outputdir+'/OverlappingPeaks/UserGiven'
                        )

        #-- annotation of peaks: user provided peak file and/or annotation bed files
        #_____________________________________________________________________________________________________
        if mode == 'UserAnnotation':
            totalCmd = ann.UserAnn (outputdir,annpeakFiles,annFiles,annSize,annName,annPrefix,threads)
            print(totalCmd)
            with Pool (threads) as p:
                p.starmap(universal.run_cmd,
                          zip(totalCmd,
                              repeat(outputdir)
                             )
                         )
        if mode == 'annotation':
          print(colored("using genome version: "+cGVersion+" in annotation mode",
                        "green",
                        attrs= ["bold"]
                       )
               )
          ann.Ann (annpeakFiles, annPrefix, cGVersion, sFasta, outputdir, threads)


        #-- Coverage of the tracks with or without spikeIn normalization
        #_____________________________________________________________________________________________________
        if mode == 'coverageTracks':

            checkBlackList(blackListedRegions)

#            if experimentType == "greencutrun":
            Names = []
            for i in range(0,myin.shape[0]):
                Name = myin.loc[i,4]
                Names.append(Name)
            covTracks.covTracks (Names,
                                 threads,
                                 outputdir,
                                 covSpike,
                                 blackListedRegions,
                                 effectiveGenomeSize,
                                 libraryType,
                                 covOtherOptions,
                                 covSpike_NormalizationFormula)

        #-- Cut frequency mode
        #_____________________________________________________________________________________________________
        if mode == "cutfrequency":

            print(colored("using genome version: "+cGVersion+" in cutfrequency mode",
                          "green",
                          attrs= ["bold"]
                         )
                 )

            e=0
            if libraryType != "pair":
                print(colored("cutfrequency mode only work with pair-end datasets ",
                              'green', attrs=['bold']))
                e=e+1

            if cCenter == "None":
                print(colored("--cCenter is missing",
                              'green',
                              attrs=['bold']
                             )
                     )
                e=e+1

            if cMotif == "False" and cMotifFile == "None":
                print(colored("--cMotifFile is missing with --cMotif False",
                              'green',
                              attrs=['bold']
                             )
                     )
                e=e+1

            if cMaN == "None":
                print (colored("--cMaN is missing",
                               "green",
                               attrs = ["bold"]
                              )
                      )
                e=e+1

                if e > 0:
                    exit()

            if cMotif == "True":
                file=greenCutFrq.initgreenCutFrq (cMaN,
                                                  cGVersion,
                                                  outputdir,
                                                  str(threads))

#                if experimentType == "greencutrun":
                Names = []
                for i in range(0,myin.shape[0]):
                    Name = myin.loc[i,4]
                    Names.append(Name)
                greenCutFrq.cutFreq (Names,
                                     threads,
                                     outputdir,
                                     file,
                                     cCenter,
                                     cMaN)

            else:

                if e > 0:
                    exit()

#                if experimentType == "greencutrun":
                Names = []
                for i in range(0,myin.shape[0]):
                    Name = myin.loc[i,4]
                    Names.append(Name)
                greenCutFrq.cutFreq (Names,
                                     str(threads),
                                     outputdir,
                                     cMotifFile,
                                     cCenter,
                                     cMaN)

        #-- Heatmap mode
        #_____________________________________________________________________________________________________

        if mode == "initHeatmap":

            checkBlackList(blackListedRegions)

            for i in range(0,myin.shape[0]):
                Name=myin.loc[i,4]
                heatMap.initHeatMap (Name,
                                     threads,
                                     outputdir,
                                     blackListedRegions,
                                     hSpike,
                                     libraryType,
                                     str(effectiveGenomeSize),
                                     initHeatmapOtherOptions
                                    )
        if mode == "heatmap":

            if hInNames == "NA":
                hInName = []
                for i in range(0,myin.shape[0]):
                    hInName.append(myin.loc[i,4])
            else:
                hInName=hInNames

            if hInFiles == "NA":
                if hCovComp == "coverage":
                    hInFile = []
                    if hCovMethod == 1:
                        for i in range(0,myin.shape[0]):
                            hInFile.append(outputdir+"/bamcoverage/"+ myin.loc[i,4] + "_expr.ScaledSpikeIn.bw")
                    elif hCovMethod == 2:
                        for i in range(0,myin.shape[0]):
                            hInFile.append(outputdir+"/bamcoverage/"+ myin.loc[i,4] + "_expr.bw")
                elif hCovComp == "compare":
                    hInFile = []
                    for i in range(0,myin.shape[0]):
                        hInFile.append(outputdir+"/bamcompare/"+ myin.loc[i,4] + ".bw")
            else:
                hInFile=hInFiles

            print("================")
            print(hInFile)
            print(hInName)
            print("================")

            #-----
            if hInFiles == "NA" and hCovComp == "coverage" and hCovMethod == 2:
                h_samplesExpr=[]
                for i in range(0,myin.shape[0]):
                    Name = myin.loc[i,4]
                    h_samplesExpr.append(outputdir + '/SpikeIn/' + Name + '_expr.Flagstats.txt')
                vals = initPeakCalling.spike_normalization2(h_samplesExpr, libraryType)
                hInCounts = [x / 10000 for x  in vals]

            elif  hInFiles != "NA":
                if hInCounts == "NA":
                    print (colored("Using total reads for normalizing heatmap",
                                   "green",
                                   attrs = ["bold"]
                                  ))
                else:
                    if len(hInCounts.split(",")) != len(hInFiles.split(",")):
                        print (colored("Length of hInCounts is not equal to hInFiles",
                                       "green",
                                       attrs = ["bold"]
                                      ))
                    print (colored("Using normalization counts given in --hInCounts",
                                   "green",
                                   attrs = ["bold"]
                                  ))
            #------


            heatMap.heatmap (hInFile,
                             hInName,
                             threads,
                             outputdir,
                             hRegionMode,
                             hGtf,
                             hBed,
                             hCovComp,
                             hCovMethod,
                             blackListedRegions,
                             hDiffPeaks,
                             hMOpt,
                             hPOpt,
                             hInCounts
                             )

        #-- MassSpectrometry or piggy back binding event mode
        #_____________________________________________________________________________________________________
        if mode == "MassSpectrometry" or mode == "piggyBack":
            if mPeak == "NA":
                for i in len(0,myin.shape[0]):
                    mName = myin.loc[i,4]
                    peakFile =outputdir + '/Peaks/' + mName +'.Clean.bed'

                    massSpectro.piggBack (outputdir,
                                          mPwm,
                                          sProt,
                                          Species,
                                          sFasta,
                                          lProt,
                                          peakFile,
                                          mName,
                                          sDist,
                                          threads,
                                          gVer)
            else:
                if mPrefix == "NA":
                    print (colored("if giving your own bed files using --mPeak"+
                                   " then give prefix of the output using --mPrefix",
                                   "green",
                                   attrs = ["bold"]
                                  )
                          )
                    exit()
                peakFiles=mPeak.split(',')
                mName=mPrefix.split(",")

                for i in range(0,len(peakFiles)):
                    massSpectro.piggBack (outputdir,
                                          mPwm,
                                          sProt,
                                          Species,
                                          sFasta,
                                          lProt,
                                          peakFiles[i],
                                          mName[i],
                                          sDist,
                                          threads,
                                          gVer)


if __name__ == "__main__":
    main()
