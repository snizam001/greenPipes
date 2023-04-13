#-- Call peaks
#____________________________________

import os
import subprocess
from datetime import datetime
from termcolor import colored
from multiprocessing import Pool
from itertools import repeat
from greenPipe import qcTagD
import pandas as pd
from greenPipe import universal
from greenPipe import initPeakCalling
import glob

def run_cmd_file (mycmd,f,outputdir):
    print(colored(mycmd, 'green', attrs=['bold']))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open(f,'w') as oFile:
            process = subprocess.Popen(mycmd,
                                       stdout=oFile,
                                       stderr=logfile)
            stdout, stderr = process.communicate()
            stdout, stderr
#--- Analysis
def callPeaksHomer (outputdir,Names,threads,ControlPeak,styles,blackListedRegions,fdr_homer,pvalue_homer,foldChange_homer,pDistHomer,pOpts):
    if pOpts == "None":
        pOpt = []
    else:
        pOpt = pOpts.replace("[","").replace("]","").split(',')

    totalCmd=[]
    totalCmdF=[]
    totalCmd2=[]
    totalCmdFile=[]
    totalCmdFboth=[]
    totalCmdFileboth=[]

    dirs=[outputdir+'/'+'Peaks']
    dirs=dirs+[outputdir + '/Tagdirectories/']
    cmd_rs=['findPeaks','pos2bed.pl','intersectBed']

    #---
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is part of Homer or bedtools. It is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH',
                          'green', attrs=['bold']))
            exit()

    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    if len(styles.split(",")) > 1:
        if len(Names) != len(styles.split(",")):
            print(colored("number of the samples and style are not equal",
                          "green",
                          attrs = ["bold"]
                         )
                 )

    s = -1
    for Name in Names:

        #----
        s = s + 1
        if len(styles.split(',')) > 1:
            style = styles.split(",") [s]
        else:
            style = styles.split(",")[0]

        #----
        inFileExpr=outputdir + '/Tagdirectories/' + Name + '_expr'
        inFileCtrl=outputdir + '/Tagdirectories/' + Name + '_control'

        tagDirs=[inFileExpr,inFileCtrl]
        for tagDir in tagDirs:
            if not os.path.exists(tagDir):
                print(colored(tagDir+
                          ': Did you performed initialPeakCalling? The tagdirectory does not exist',
                          'green', attrs=['bold']))

        #----
        outPeak=outputdir + '/Peaks/' + Name

        #----
        narrow_on=0
        broad_on=0
        if style=='narrow' or style=='factor':
            narrow_on = 1
        if style=='broad' or style=='histone':
            broad_on = 1
        if style == "both":
            narrow_on = 1
            broad_on  = 1

        if narrow_on == 1:
            if ControlPeak == 'True':
                cmd=['findPeaks',
                     inFileExpr,
                     '-i', inFileCtrl,
                     '-style', 'factor',
                     '-C', '0',
                     '-o', outPeak + "_narrow-homer",
                     '-'+pDistHomer, str(fdr_homer),
                     '-P', str(pvalue_homer),
                     '-F', str(foldChange_homer)] + pOpt
            else:
                cmd=['findPeaks', inFileExpr,
                     '-style', 'factor',
                     '-C', '0',
                     '-o', outPeak + "_narrow-homer",
                     '-'+pDistHomer, str(fdr_homer),
                     '-P', str(pvalue_homer),
                     '-F', str(foldChange_homer)] + pOpt
            totalCmd.append(cmd)
            #-------
            cmd=['pos2bed.pl',
                 '-o', outPeak+'_narrow-homer.bed',
                 outPeak+"_narrow-homer"]
            totalCmd2.append(cmd)
            #-------
            cmd=['intersectBed',
                 '-a', outPeak+'_narrow-homer.bed',
                 '-b', blackListedRegions,
                 '-v']
            f=outPeak+'_narrow-homer.removed.bed'
            totalCmdF.append(cmd)
            totalCmdFile.append(f)
            #-------
            cmd=['grep', '^chr', outPeak+'_narrow-homer.removed.bed']
            f=outPeak+'_narrow-homer.Clean.bed'
            totalCmdF.append(cmd)
            totalCmdFile.append(f)

        if broad_on==1:
            if ControlPeak == 'True':
                cmd=['findPeaks', inFileExpr,
                     '-i', inFileCtrl,
                     '-style', 'histone',
                     '-C', '0',
                     '-o', outPeak + "_broad-homer",
                     '-'+pDistHomer, str(fdr_homer),
                     '-P', str(pvalue_homer),
                     '-F', str(foldChange_homer)] + pOpt
            else:
                cmd=['findPeaks', inFileExpr,
                     '-style', 'histone',
                     '-C', '0',
                     '-o', outPeak + "_broad-homer",
                     '-'+pDistHomer, str(fdr_homer),
                     '-P', str(pvalue_homer),
                     '-F', str(foldChange_homer)] + pOpt
            totalCmd.append(cmd)

            #-------
            cmd=['pos2bed.pl',
                 '-o', outPeak+'_broad-homer.bed',
                 outPeak+"_broad-homer"]
            totalCmd2.append(cmd)
            #-------
            cmd=['intersectBed',
                 '-a', outPeak+'_broad-homer.bed',
                 '-b', blackListedRegions,
                 '-v']
            f=outPeak+'_broad-homer.removed.bed'
            totalCmdF.append(cmd)
            totalCmdFile.append(f)
            #-------
            cmd=['grep', '^chr', outPeak+'_broad-homer.removed.bed']
            f=outPeak+'_broad-homer.Clean.bed'
            totalCmdF.append(cmd)
            totalCmdFile.append(f)

        if style == "both":
            totalCmdFboth.append(['cat',outPeak+'_broad-homer.Clean.bed',outPeak+'_narrow-homer.Clean.bed'])
            totalCmdFileboth.append(outPeak+".temp.txt")

            totalCmdFboth.append(['sort',"-k1,1","-k2,2n",outPeak+".temp.txt"])
            totalCmdFileboth.append(outPeak+".temp2.txt")

            totalCmdFboth.append(["bedtools","merge","-i",outPeak+".temp2.txt"])
            totalCmdFileboth.append(outPeak+'_all-homer.Clean.bed')

    #----
    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )

    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd2,
                  repeat(outputdir)
                 )
                 )

    with Pool (1) as p:
        p.starmap(run_cmd_file,
              zip(totalCmdF,
                  totalCmdFile,
                  repeat(outputdir)
                 )
                 )

    if style == "both":
        with Pool (1) as p:
            p.starmap(run_cmd_file,
                  zip(totalCmdFboth,
                      totalCmdFileboth,
                      repeat(outputdir)
                     )
                     )
        rmFiles = glob.glob(outputdir + '/Peaks/' + "*temp*.txt")
        for rmFile in rmFiles:
            cmd=['rm',rmFile]
            universal.run_cmd(cmd,outputdir)

def callPeaksMacs2 (outputdir,Names,threads,spikePeaks,ControlPeak,styles,effectiveGenomeSize,blackListedRegions,libraryType,pOpts):

    if pOpts == "None":
        pOpt = []
    else:
        pOpt = pOpts.replace("[","").replace("]","").split(',')

    totalCmd=[]
    totalCmdF=[]
    totalCmdFboth=[]
    totalCmdFileboth=[]
    totalCmdFile=[]
    totalCmd2=[]

    cmd_rs=['macs2','pos2bed.pl','intersectBed']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is part of macs or bedtools. It is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH',
                          'green', attrs=['bold']))
            exit()

    dirs=[outputdir+'/'+'Peaks']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    totalFile=[]
    for Name in Names:
        totalFile.append(outputdir + '/Bamfiles/' + Name + '_expr.bam')
        totalFile.append(outputdir + '/Bamfiles/' + Name + '_control.bam')
    totalFile.append(blackListedRegions)
    exist_file=0
    for checkFile in totalFile:
        if not os.path.exists(checkFile):
            exist_file=exist_file+1
            print(colored(blackListedRegions+": file does not exist",
                          'green',
                          attrs=['bold']
                         )
                 )
    if exist_file > 0:
        exit()

    if len(styles.split(",")) > 1:
        if len(Names) != len(styles.split(",")):
            print(colored("number of the samples and style are not equal",
                          "green",
                          attrs = ["bold"]
                         )
                 )

    s = -1
    for Name in Names:
        #----
        s = s + 1
        if len(styles.split(',')) > 1:
            style = styles.split(",") [s]
        else:
            style = styles.split(",")[0]


        #----
        narrow_on=0
        broad_on=0
        if style=='narrow' or style=='factor':
            narrow_on = 1
        if style=='broad' or style=='histone':
            broad_on = 1
        if style == "both":
            narrow_on = 1
            broad_on  = 1

        inExpr=outputdir + '/Bamfiles/' + Name + '_expr.bam'
        inCtrl=outputdir + '/Bamfiles/' + Name + '_control.bam'
        if spikePeaks=='False':
            if broad_on==1:
                c=['macs2',
                   'callpeak',
                   '-t', inExpr,
                   '-c', inCtrl,
                   '-f', 'BAM',
                   '-g', str(effectiveGenomeSize),
                   '-n', Name + "_broad-macs2",
                   '--outdir', dirs[0],
                   '--keep-dup','all',
                   '--seed','786',
                   '--broad'] + pOpt
                totalCmd.append(c)

            elif narrow_on == 1:
                c=['macs2',
                   'callpeak',
                   '-t', inExpr,
                   '-c', inCtrl,
                   '-f', 'BAM',
                   '-g', str(effectiveGenomeSize),
                   '-n', Name + "_narrow-macs2",
                   '--outdir', dirs[0],
                   '--keep-dup','all',
                   '--seed','786'] + pOpt
                totalCmd.append(c)

        elif ControlPeak=='False':
            if broad_on==1:
                c=['macs2',
                   'callpeak',
                   '-t', inExpr,
                   '-f', 'BAM',
                   '-g', str(effectiveGenomeSize),
                   '-n', Name + "_broad-macs2",
                   '--outdir', dirs[0],
                   '--keep-dup', 'all',
                   '--seed','786',
                   '--broad']
                totalCmd.append(c)

            elif narrow_on == 1:
                c=['macs2',
                   'callpeak',
                   '-t', inExpr,
                   '-f', 'BAM',
                   '-g', str(effectiveGenomeSize),
                   '-n', Name + "_narrow-macs2",
                   '--outdir', dirs[0],
                   '--keep-dup', 'all',
                   '--seed','786'] + pOpt
                totalCmd.append(c)

        else:
            sCtrl =outputdir + '/SpikeIn/' + Name + '_control.Flagstats.txt'
            sExpr =outputdir + '/SpikeIn/' + Name + '_expr.Flagstats.txt'
            Ctrl  =outputdir + '/Bamfiles/' + Name + '_control.Flagstats.txt'
            Expr  =outputdir + '/Bamfiles/' + Name + '_expr.Flagstats.txt'
            allFs=[sCtrl,sExpr,Ctrl,Expr]
            allF_exit = 0
            for allF in allFs:
                if not os.path.exists(allF):
                    print(colored(allF+" : file does not exist. This is required to calculate normalization ratio."+
                                  "green",
                                  attrs = ['bold']
                                 )
                         )

            oVal  =initPeakCalling.spike_normalization2([sCtrl,sExpr,Ctrl,Expr],libraryType)
            print(oVal)
            oVal_norm = (oVal[0]/oVal[2])/(oVal[1]/oVal[3])
            myratio=1/oVal_norm

            if broad_on==1:
                c=['macs2',
                   'callpeak',
                   '-t', inExpr,
                   '-c', inCtrl,
                   '-f', 'BAM',
                   '-g', str(effectiveGenomeSize),
                   '-n', Name + "_broad-macs2",
                   '--outdir', dirs[0],
                   '--keep-dup', 'all',
                   '--seed', '786',
                   '--broad',
                   '--ratio', str(myratio)] + pOpt
                totalCmd.append(c)

            elif narrow_on == 1:
                c=['macs2',
                   'callpeak',
                   '-t', inExpr,
                   '-c', inCtrl,
                   '-f', 'BAM',
                   '-g', str(effectiveGenomeSize),
                   '-n', Name + "_narrow-macs2",
                   '--outdir', dirs[0],
                   '--keep-dup', 'all',
                   '--seed', '786',
                   '--ratio', str(myratio)] + pOpt
                totalCmd.append(c)

        if broad_on == 1:
            c=['intersectBed',
               '-a', dirs[0] + '/' + Name + '_broad-macs2_peaks.broadPeak',
               '-b', blackListedRegions,
               '-v']
            totalCmdFile.append(dirs[0] + '/' + Name + '_broad-macs2.removed.bed')
            totalCmdF.append(c)

            c=['grep', '^chr',
               dirs[0] + '/' + Name + '_broad-macs2.removed.bed']
            totalCmdFile.append(dirs[0] + '/' + Name + '_broad-macs2.Clean.bed')
            totalCmdF.append(c)

        if narrow_on == 1:
            c=['intersectBed',
               '-a', dirs[0] + '/' + Name + '_narrow-macs2_summits.bed',
               '-b', blackListedRegions,
               '-v']
            totalCmdFile.append(dirs[0] + '/' + Name + '_narrow-macs2.removed.bed')
            totalCmdF.append(c)

            c=['grep', '^chr',
               dirs[0] + '/' + Name + '_narrow-macs2.removed.bed']
            totalCmdFile.append(dirs[0] + '/' + Name + '_narrow-macs2.Clean.bed')
            totalCmdF.append(c)

        if style == "both":
            totalCmdFboth.append(['cat',
              dirs[0] + '/' + Name + '_broad-macs2.Clean.bed',
              dirs[0] + '/' + Name +'_narrow-macs2.Clean.bed'])
            totalCmdFileboth.append(dirs[0] + '/' + Name+".temp.txt")

            totalCmdFboth.append(['sort',"-k1,1","-k2,2n",dirs[0] + '/' + Name+".temp.txt"])
            totalCmdFileboth.append(dirs[0] + '/' + Name+".temp2.txt")

            totalCmdFboth.append(["bedtools","merge","-i",dirs[0] + '/' + Name+".temp2.txt"])
            totalCmdFileboth.append(dirs[0] + '/' + Name+'_all-macs2.Clean.bed')

    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )

    print(totalCmdF)
    with Pool (1) as p:
        p.starmap(run_cmd_file,
              zip(totalCmdF,
                  totalCmdFile,
                  repeat(outputdir)
                 )
                 )

    if style == "both":
        with Pool (1) as p:
            p.starmap(run_cmd_file,
                  zip(totalCmdFboth,
                      totalCmdFileboth,
                      repeat(outputdir)
                     )
                     )
        rmFiles = glob.glob(outputdir + '/Peaks/' + "*temp*.txt")
        for rmFile in rmFiles:
            cmd=['rm',rmFile]
            universal.run_cmd(cmd,outputdir)
#----------------------------------------------------------
def callPeaksSEACR (outputdir,Names,threads,spikePeaks,ControlPeak,seacr,libraryType,seacrMode,seacrThreshold,blackListedRegions,genomeFile):
    if libraryType!='pair':
        print(colored('SEACR for single end reads is not implemented',
                      'green',
                      attrs=['bold']
                     )
             )
        exit()

    totalCmd=[]
    totalCmdF=[]
    totalCmdFile=[]
    totalGnomCov=[]
    totalGnomCovF=[]

    cmd_rs=[seacr,'samtools','bedtools']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH',
                          'green', attrs=['bold']))
            exit()

    dirs=[outputdir+'/'+'Peaks']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    totalFile=[]
    for Name in Names:
        totalFile.append(outputdir + '/Bamfiles/' + Name + '_expr.bam')
        totalFile.append(outputdir + '/Bamfiles/' + Name + '_control.bam')
    totalFile.append(blackListedRegions)
    exist_file=0
    for checkFile in totalFile:
        if not os.path.exists(checkFile):
            exist_file=exist_file+1
            print(colored(blackListedRegions+": file does not exist",
                          'green',
                          attrs=['bold']
                         )
                 )
    if exist_file > 0:
        exit()

    fragDir=outputdir+'/'+'Tagdirectories_qualities/'
    if not os.path.exists(fragDir):
      os.makedirs(fragDir)

    for Name in Names:

        exprFrag=fragDir + '/' + Name + '_expr.sorted.fragments.bed'
        ctrlFrag=fragDir + '/' + Name + '_control.sorted.fragments.bed'
        outPeak=outputdir + '/Peaks/' + Name + '-seacr'
        #-------
        if not os.path.exists(exprFrag) or not os.path.exists(ctrlFrag):
            qcTagD.qcTagD_fragmentLength(outputdir,Name,threads,"False")

        totalGnomCov.append(["bedtools",
                             "genomecov",
                             "-bg",
                             "-i", exprFrag,
                             "-g", genomeFile]
                           )
        totalGnomCovF.append(fragDir + '/' + Name + '_expr.gCov.bed')

        totalGnomCov.append(["bedtools",
                             "genomecov",
                             "-bg",
                             "-i", ctrlFrag,
                             "-g", genomeFile]
                           )
        totalGnomCovF.append(fragDir + '/' + Name + '_control.gCov.bed')

        #-------
        if spikePeaks=='False':
            c=['bash',
               seacr,
               fragDir + '/' + Name + '_expr.gCov.bed',
               fragDir + '/' + Name + '_control.gCov.bed',
               'norm',
               seacrMode,
               outPeak]

        elif ControlPeak=='False':
            c=['bash',
               seacr,
               fragDir + '/' + Name + '_expr.gCov.bed',
               seacrThreshold,
               'norm',
               seacrMode,
               outPeak]
        else:
            #-------
            c=['bash',
               seacr,
               fragDir + '/' + Name + '_expr.gCov.bed',
               fragDir + '/' + Name + '_control.gCov-normalized.bed',
               'non',
               seacrMode,
               outPeak]

        totalCmd.append(c)

        #-------
        cmd=['intersectBed',
             '-a', outPeak+'.'+ seacrMode + '.bed',
             '-b', blackListedRegions,
             '-v']
        totalCmdFile.append(outPeak+'.removed.bed')
        totalCmdF.append(cmd)

        #-------
        cmd=['grep', '^chr', outPeak+'.removed.bed']
        totalCmdFile.append(outPeak+'.Clean.bed')
        totalCmdF.append(cmd)


    with Pool (threads) as p:
        p.starmap(run_cmd_file,
              zip(totalGnomCov,
                  totalGnomCovF,
                  repeat(outputdir)
                 )
                 )

    if spikePeaks!='False' or ControlPeak!='False':
        for Name in Names:
            sCtrl =outputdir + '/SpikeIn/' + Name + '_control.Flagstats.txt'
            sExpr =outputdir + '/SpikeIn/' + Name + '_expr.Flagstats.txt'
            Ctrl  =outputdir + '/Bamfiles/' + Name + '_control.Flagstats.txt'
            Expr  =outputdir + '/Bamfiles/' + Name + '_expr.Flagstats.txt'
            allFs=[sCtrl,sExpr,Ctrl,Expr]
            allF_exit = 0
            for allF in allFs:
                if not os.path.exists(allF):
                    print(colored(allF+" : file does not exist. This is required to calculate normalization ratio."+
                                  "green",
                                  attrs = ['bold']
                                 )
                         )

            oVal  =initPeakCalling.spike_normalization2([sCtrl,sExpr,Ctrl,Expr],libraryType)
            print(oVal)
            oVal_norm = (oVal[0]/oVal[2])/(oVal[1]/oVal[3])
            myratio=1/oVal_norm
            print("------------------")
            print(myratio)

            seacrBed=pd.read_csv(fragDir + '/' + Name + '_control.gCov.bed',
                sep='\t',
                header=None)
            seacrBed.head()
            seacrBed.iloc[:,3] = seacrBed.iloc[:,3]*myratio

            seacrBed.to_csv(fragDir + '/' + Name + '_control.gCov-normalized.bed',
                sep='\t',
                header=None,
                index=None)

    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )

    with Pool (1) as p:
        p.starmap(run_cmd_file,
              zip(totalCmdF,
                  totalCmdFile,
                  repeat(outputdir)
                 )
                 )
