#-- initiate_heatmap or bamcompare files preparation: here running for the single file each time
#____________________________________

import os
import pandas as pd
import subprocess
from datetime import datetime
from termcolor import colored
import glob
from greenPipe import initPeakCalling
from greenPipe import universal
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
                          'green', 
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
                          'green',
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
        
        
        oVal_norm = (oVal[0]/oVal[2])/(oVal[1]/oVal[3])
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
                   '--effectiveGenomeSize', effectiveGenomeSize ] + initHeatmapOtherOptions.split(',')

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
                   '--effectiveGenomeSize', effectiveGenomeSize ] + initHeatmapOtherOptions.split(',')

            universal.run_cmd(mycmd,outputdir)            

def heatmap (inFiles,inNames,threads,outputdir,hRegionMode,gtf,hBed,hCovComp,blackListedRegions,hDiffPeaks, hMOpt, hPOpt):
    heatmapNormalize=psource.resource_filename(__name__, "rscripts/heatmapNormalize.R")
    if hMOpt == "None":
        hMOpt = [
        '-b', '5000', 
        '-a', '5000']
    else:
        hMOpt = hMOpt.split(',')

    if hPOpt == "None":
        hPOpt = ['--colorMap', 'RdBu']
    else:
        hPOpt = hPOpt.split(',')

    if not isinstance(inFiles, list):
        inFiles = inFiles.split(",")
    if not isinstance(inNames, list):
        inNames = inNames.split(",")

    e_file=0
    if len(inFiles) != len(inNames):
        print(colored("length of the input file (--infiles) and names (--inNames) are not equal",
                      "green",
                      attrs = ["bold"]
                     )
             )
        e_file = e_file + 1
    print(inFiles)
    for inFile in inFiles:
        if not os.path.exists(inFile):
            print(colored(inFile + ":given input files (--inFiles) does not exist",
                          "green",
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
                          'green', 
                          attrs=['bold']))
            exit()
            
    #-----       
    dirs=[outputdir+'/'+'HeatMaps']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
            
    #-----
    if regionMode == "metagene":
        if inFiles == "NA":
            print (colored("--inFiles option is missing",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if gtf == "NA":
            print (colored("--gtf option is missing",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()

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
        
        c=['plotHeatmap', 
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.gz', 
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.jpeg', 
           '--dpi', '300'] + hPOpt
        
        universal.run_cmd(c,outputdir)

        inCounts = [1]*len(inNames)

        c=[heatmapNormalize, 
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metagene.gz'
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-metageneCoverage'
        '-l', inNames,
        '-c', inCounts
        ]

        universal.run_cmd(c,outputdir)
    #------------
    elif regionMode == "tss":
        if inFiles == "NA":
            print (colored("--inFiles option is missing",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if gtf == "NA":
            print (colored("--gtf option is missing",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()

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
        
        c=['plotHeatmap', 
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.gz', 
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.jpeg',
           '--dpi', '300', 
           '--refPointLabel', "TSS"] + hPOpt
        
        universal.run_cmd(c,outputdir) 

        inCounts = [1]*len(inNames)

        c=[heatmapNormalize, 
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss.gz'
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-tss'
        '-l', inNames,
        '-c', inCounts
        ]
        
        universal.run_cmd(c,outputdir)
    #------------
    elif regionMode == "bed":
        if inFiles == "NA":
            print (colored("--inFiles option is missing",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()
        if hBed == "NA":
            print (colored("--hBed option is missing",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()

        c=['computeMatrix', 
           'reference-point',
           '--referencePoint',
           'center', 
           '-R', hBed, 
           '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz', 
           '--missingDataAsZero', 
           '-bl', blackListedRegions, 
           '--smartLabels', 
           '-p', str(threads),
           '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

        universal.run_cmd(c,outputdir)
        
        c=['plotHeatmap', 
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz', 
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.jpeg',
           '--dpi', '300', 
           '--refPointLabel', "center"] + hPOpt
        
        universal.run_cmd(c,outputdir) 

        inCounts = [1]*len(inNames)

        c=[heatmapNormalize, 
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed.gz'
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-bed'
        '-l', inNames,
        '-c', inCounts
        ]
        
        universal.run_cmd(c,outputdir)
    #------------
    elif regionMode == "peaks":
        hBed=[]
        for inName in inNames:
            hBed.append(outputdir + '/Peaks/' + inName + '.Clean.bed')
        e_file=0
        for f in hBed:
            if not os.path.exists(f):
                print(colored(f+" : file does not exit. These files require in the --regionMode peaks." + 
                              "Call peaks before this step.",
                              "green",
                              attrs = ["bold"]
                             )
                     )
                e_file = e_file + 1 
        if e_file > 0:
            exit()
            
        c=['computeMatrix', 
           'reference-point',
           '--referencePoint',
           'center', 
           '-R', hBed, 
           '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz', 
           '--missingDataAsZero', 
           '-bl', blackListedRegions, 
           '--smartLabels', 
           '-p', str(threads),
           '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

        universal.run_cmd(c,outputdir)
        
        c=['plotHeatmap', 
           '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz', 
           '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.jpeg',
           '--dpi', '300', 
           '--refPointLabel', "center"] + hPOpt
        
        universal.run_cmd(c,outputdir)   

        inCounts = [1]*len(inNames)

        c=[heatmapNormalize, 
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks.gz'
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-peaks'
        '-l', inNames,
        '-c', inCounts
        ]
        
        universal.run_cmd(c,outputdir)
        
#--
        if hDiffPeaks == "True":
            hBed=glob.glob(outputdir+ "/Peaks/DifferentialPeaks/" + "*.differentialPeaks.txt")
            if len(hBed) == 0:
                print (colored("You wanted to generate heatmap of the differential peak, but it"+
                               " seems that files does not exist. Can you check *differentialPeaks.txt files" +
                               " in the folder " + outputdir+ "/Peaks/DifferentialPeaks/"+
                               " or use --hDiffPeaks False",
                               "green",
                               attrs = ["bold"]
                              )
                      )
                exit()
                
            c=['computeMatrix', 
               'reference-point', 
               '--referencePoint',
               'center', 
               '-R', hBed, 
               '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz', 
               '--missingDataAsZero', 
               '-bl', blackListedRegions, 
               '--smartLabels', 
               '-p', str(threads),
               '--samplesLabel'] +  inNames + ['-S'] + inFiles + hMOpt

            universal.run_cmd(c,outputdir)

            c=['plotHeatmap', 
               '-m', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz', 
               '-out', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.jpeg',
               '--dpi', '300', 
               '--refPointLabel', "center"] + hPOpt 

            universal.run_cmd(c,outputdir)   

        inCounts = [1]*len(inNames)

        c=[heatmapNormalize, 
        '-i', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks.gz'
        '-o', outputdir+'/'+'HeatMaps/heatmap-'+hCovComp+'-differentialpeaks'
        '-l', inNames,
        '-c', inCounts
        ]
        
        universal.run_cmd(c,outputdir)
