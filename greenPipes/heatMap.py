#-- initiate_heatmap or bamcompare files preparation: here running for the single file each time
#____________________________________

import os
import pandas as pd
import subprocess
from datetime import datetime
from termcolor import colored
import glob
from greenPipes import initPeakCalling
from greenPipes import universal
import pkg_resources as psource

def initHeatMap (Name,threads,outputdir,blackListedRegions,heatmapSpikeIn,libraryType,effectiveGenomeSize, initHeatmapOtherOptions):
    cmd_rs=['bamCompare']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.'+
                          " This tools is the part of deepTools. Please install it.",
                          'red',
                          attrs=['bold']))
            exit()

    dirs=[outputdir+'/'+'bamcompare']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    experiments = ['_expr','_control']
    e_file=0
    for experiment in experiments:
        b_file=outputdir+ '/Bamfiles/' + Name + experiment + '.bam'
        if not os.path.exists(b_file):
            e_file=e_file+1
            print(colored(b_file+": does not exist",
                          'red',
                          attrs=['bold']
                         )
                 )
    if e_file > 0:
        exit("First align the reads before generating heatmaps")


    #-----
    eFile=outputdir + '/Bamfiles/' + Name + '_expr.bam'
    cFile=outputdir + '/Bamfiles/' + Name + '_control.bam'
    oFile=outputdir + '/bamcompare/' + Name + '.bw'

    #----
    if heatmapSpikeIn=='True':

        sCtrl =outputdir + '/SpikeIn/' + Name + '_control.Flagstats.txt'
        sExpr =outputdir + '/SpikeIn/' + Name + '_expr.Flagstats.txt'
        Ctrl  =outputdir + '/Bamfiles/' + Name + '_control.Flagstats.txt'
        Expr  =outputdir + '/Bamfiles/' + Name + '_expr.Flagstats.txt'

        oVal  =initPeakCalling.spike_normalization2([sCtrl,sExpr,Ctrl,Expr],libraryType)
        print(colored("sCtrl sExpr sCtrl sExpr","green",attrs=["bold"]))
        print(oVal)

        # This one is changed on the 12 May, 2023
        #oVal_norm = (oVal[0]/oVal[2])/(oVal[1]/oVal[3])
        oVal_norm = (oVal[0]/(oVal[0]+oVal[2]))/(oVal[1]/(oVal[1]+oVal[3]))
        oVal_norm = 1/oVal_norm
        if initHeatmapOtherOptions == "None":
            mycmd=['bamCompare',
                   '-b1', eFile,
                   '-b2', cFile,
                   '-o', oFile,
                   '-p', str(threads),
                   '-bl', blackListedRegions,
                   '--scaleFactors', '1:'+str(oVal_norm),
                   '--effectiveGenomeSize', effectiveGenomeSize ]

            universal.run_cmd(mycmd,outputdir)
        else:
            mycmd=['bamCompare',
                   '-b1', eFile,
                   '-b2', cFile,
                   '-o', oFile,
                   '-p', str(threads),
                   '-bl', blackListedRegions,
                   '--scaleFactors', '1:'+str(oVal_norm),
                   '--effectiveGenomeSize', effectiveGenomeSize ] + initHeatmapOtherOptions.replace("[","").replace("]","").split(',')

            universal.run_cmd(mycmd,outputdir)

    if heatmapSpikeIn=='False':
        if initHeatmapOtherOptions == "None":
            mycmd=['bamCompare',
                   '-b1', eFile,
                   '-b2', cFile,
                   '-o', oFile,
                   '-p', str(threads),
                   '-bl', blackListedRegions,
                   '--effectiveGenomeSize', effectiveGenomeSize ]

            universal.run_cmd(mycmd,outputdir)
        else:
            mycmd=['bamCompare',
                   '-b1', eFile,
                   '-b2', cFile,
                   '-o', oFile,
                   '-p', str(threads),
                   '-bl', blackListedRegions,
                   '--effectiveGenomeSize', effectiveGenomeSize ] + initHeatmapOtherOptions.replace("[","").replace("]","").split(',')

            universal.run_cmd(mycmd,outputdir)

def heatmap (inFiles,inNames,threads,outputdir,hRegionMode,gtf,hBed,hCovComp, hCovMethod, blackListedRegions,hDiffPeaks, hMOpt, hPOpt, inCounts, hPeakType):
    heatmapNormalize=psource.resource_filename(__name__, "rscripts/heatmapNormalize.R")

    if hMOpt == "None":
        hMOpt = [
        '-b', '5000',
        '-a', '5000']
    else:
        hMOpt = hMOpt.replace("[","").replace("]","").split(',')
        print("Using following options in heatmap production (computeMatrix)")
        print(hMOpt)

    if hPOpt == "None":
        hPOpt = ['--colorMap', 'GnBu',
        '--yAxisLabel', 'Coverage',
        '--heatmapWidth', '1.5',
        '--heatmapHeight', '5']

    else:
        hPOpt = hPOpt.replace("[","").replace("]","").split(',')

    if not isinstance(inFiles, list):
        inFiles = inFiles.split(",")
    if not isinstance(inNames, list):
        inNames = inNames.split(",")

    e_file=0
    if len(inFiles) != len(inNames):
        print(colored("length of the input file (--infiles) and names (--inNames) are not equal",
                      "red",
                      attrs = ["bold"]
                     )
             )
        e_file = e_file + 1
    print(inFiles)
    for inFile in inFiles:
        if not os.path.exists(inFile):
            print(colored(inFile + ":given input files (--inFiles) does not exist",
                          "red",
                          attrs = ["bold"]
                         )
                 )
            e_file = e_file + 1
    if e_file > 0:
        exit()


    cmd_rs=['computeMatrix','plotHeatmap']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.'+
                          " This tools is the part of deepTools. Please install it.",
                          'red',
                          attrs=['bold']))
            exit()

    #-----
    dirs=[outputdir+'/'+'HeatMaps']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    #-----
    if hRegionMode == "metagene":
        if inFiles == "NA":
            print (colored("--inFiles option is missing",
                           "red",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if gtf == "NA":
            print (colored("--gtf option is missing",
                           "red",
                           attrs = ["bold"]
                          )
                  )
            exit()

        if hCovMethod == 1:
            inCounts = [1]*len(inNames)
            inCounts = [str(inCounts) for inCounts in inCounts]
        elif hCovMethod == 2:
            inCounts = [str(inCounts) for inCounts in inCounts]

        if hCovMethod == 1:
            c=['computeMatrix',
               'scale-regions',
               '-R', gtf,
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--metagene',
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        if hCovMethod == 2:
            c=['computeMatrix',
               'scale-regions',
               '-R', gtf,
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--metagene',
               '--averageTypeBins', 'sum',
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.jpeg',
           '--dpi', '300'] + hPOpt

        universal.run_cmd(c,outputdir)

        c=['Rscript', heatmapNormalize,
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.gz',
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metageneCoverage',
        '-l', ",".join(inNames),
        '-c', ",".join(inCounts)
        ]

        universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene-normalized.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene-normalized.jpeg',
           '--dpi', '300',
           '--refPointLabel', "center"] + hPOpt

        universal.run_cmd(c,outputdir)

    #------------
    elif hRegionMode == "tss":
        if inFiles == "NA":
            print (colored("--inFiles option is missing",
                           "red",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if gtf == "NA":
            print (colored("--gtf option is missing",
                           "red",
                           attrs = ["bold"]
                          )
                  )
            exit()

        if hCovMethod == 1:
            inCounts = [1]*len(inNames)
            inCounts = [str(inCounts) for inCounts in inCounts]
        elif hCovMethod == 2:
            inCounts = [str(inCounts) for inCounts in inCounts]

        if hCovMethod == 1:
            c=['computeMatrix',
               'reference-point',
               '--referencePoint',
               'TSS',
               '-R', gtf,
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        elif hCovMethod == 2:
            c=['computeMatrix',
               'reference-point',
               '--referencePoint',
               'TSS',
               '-R', gtf,
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--averageTypeBins', 'sum',
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.jpeg',
           '--dpi', '300',
           '--refPointLabel', "TSS"] + hPOpt

        universal.run_cmd(c,outputdir)

        c=['Rscript', heatmapNormalize,
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.gz',
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss',
        '-l', ",".join(inNames),
        '-c', ",".join(inCounts)
        ]

        universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss-normalized.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss-normalized.jpeg',
           '--dpi', '300',
           '--refPointLabel', "center"] + hPOpt

        universal.run_cmd(c,outputdir)

    #------------
    elif hRegionMode == "bed":
        if inFiles == "NA":
            print (colored("--inFiles option is missing",
                           "red",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if hBed == "NA":
            print (colored("--hBed option is missing",
                           "red",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if hCovMethod == 1:
            inCounts = [1]*len(inNames)
            inCounts = [str(inCounts) for inCounts in inCounts]
        elif hCovMethod == 2:
            inCounts = [str(inCounts) for inCounts in inCounts]

        if hCovMethod == 1:
            c=['computeMatrix',
               'reference-point',
               '--referencePoint',
               'center',
               '-R'] + hBed.split(',') + [
                '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        elif hCovMethod == 2:
            c=['computeMatrix',
               'reference-point',
               '--referencePoint',
               'center',
               '-R'] + hBed.split(',') + [
                '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--averageTypeBins', 'sum',
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.jpeg',
           '--dpi', '300',
           '--refPointLabel', "center"] + hPOpt

        universal.run_cmd(c,outputdir)

        c=['Rscript', heatmapNormalize,
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz',
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed',
        '-l', ",".join(inNames),
        '-c', ",".join(inCounts)
        ]

        universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed-normalized.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed-normalized.jpeg',
           '--dpi', '300',
           '--refPointLabel', "center"] + hPOpt

        universal.run_cmd(c,outputdir)

    #------------
    elif hRegionMode == "peaks":
        print (colored("Using peak type: "+ hPeakType,
                       "green",
                       attrs = ["bold"]
                      )
              )
        hBed=[]
        if 'homer' in hPeakType:
            for inName in inNames:
                hBed.append(outputdir + '/Peaks/' + inName + '_' + hPeakType + '.Clean.bed')
        else:
            for inName in inNames:
                hBed.append(outputdir + '/Peaks/' + inName + '-' + hPeakType + '.bed')
        e_file=0
        for f in hBed:
            if not os.path.exists(f):
                print(colored(f+" : file does not exit. These files require in the --hRegionMode peaks." +
                              "Call peaks before this step.",
                              "red",
                              attrs = ["bold"]
                             )
                     )
                e_file = e_file + 1
        if e_file > 0:
            exit()
        print(hBed)
        if hCovMethod == 1:
            inCounts = [1]*len(inNames)
            inCounts = [str(inCounts) for inCounts in inCounts]
        elif hCovMethod == 2:
            inCounts = [str(inCounts) for inCounts in inCounts]

        if hCovMethod == 1:
            c=['computeMatrix',
               'reference-point',
               '--referencePoint',
               'center',
               '-R'] + hBed + [
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        if hCovMethod == 2:
            c=['computeMatrix',
               'reference-point',
               '--referencePoint',
               'center',
               '-R'] + hBed + [
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz',
               '--missingDataAsZero',
               '-bl', blackListedRegions,
               '--smartLabels',
               '-p', str(threads),
               '--averageTypeBins', 'sum',
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.jpeg',
           '--dpi', '300',
           '--refPointLabel', "center"] + hPOpt

        universal.run_cmd(c,outputdir)
        c=['Rscript', heatmapNormalize,
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz',
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks',
        '-l', ",".join(inNames),
        '-c', str(",".join(inCounts))
        ]

        universal.run_cmd(c,outputdir)

        c=['plotHeatmap',
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks-normalized.gz',
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks-normalized.jpeg',
           '--dpi', '300',
           '--refPointLabel', "center"] + hPOpt

        universal.run_cmd(c,outputdir)

#--
        if hDiffPeaks == "True":
            hBed=glob.glob(outputdir+ "/Peaks/DifferentialPeaks/" + "*.differentialPeaks.txt")
            if len(hBed) == 0:
                print (colored("You wanted to generate heatmap of the differential peak, but it"+
                               " seems that files does not exist. Can you check *differentialPeaks.txt files" +
                               " in the folder " + outputdir+ "/Peaks/DifferentialPeaks/"+
                               " or use --hDiffPeaks False",
                               "red",
                               attrs = ["bold"]
                              )
                      )
                exit()


            if hCovMethod == 1:
                inCounts = [1]*len(inNames)
                inCounts = [str(inCounts) for inCounts in inCounts]
            elif hCovMethod == 2:
                inCounts = [str(inCounts) for inCounts in inCounts]

            if hCovMethod == 1:
                c=['computeMatrix',
                   'reference-point',
                   '--referencePoint',
                   'center',
                   '-R'] + hBed + [
                   '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz',
                   '--missingDataAsZero',
                   '-bl', blackListedRegions,
                   '--smartLabels',
                   '-p', str(threads),
                   '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

                universal.run_cmd(c,outputdir)

            if hCovMethod == 2:
                c=['computeMatrix',
                   'reference-point',
                   '--referencePoint',
                   'center',
                   '-R'] + hBed + [
                   '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz',
                   '--missingDataAsZero',
                   '-bl', blackListedRegions,
                   '--smartLabels',
                   '-p', str(threads),
                   '--averageTypeBins', 'sum',
                   '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

                universal.run_cmd(c,outputdir)

            c=['plotHeatmap',
               '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz',
               '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.jpeg',
               '--dpi', '300',
               '--refPointLabel', "center"] + hPOpt

            universal.run_cmd(c,outputdir)

            c=['Rscript', heatmapNormalize,
            '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz',
            '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks',
            '-l', ",".join(inNames),
            '-c', ",".join(inCounts)
            ]

            universal.run_cmd(c,outputdir)


            c=['plotHeatmap',
               '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks-normalized.gz',
               '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks-normalized.jpeg',
               '--dpi', '300',
               '--refPointLabel', "center"] + hPOpt

            universal.run_cmd(c,outputdir)
