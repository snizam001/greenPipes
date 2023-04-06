#-- IDRs

from greenPipe import universal
from greenPipe import initPeakCalling
import pkg_resources as psource
import glob
import os
import subprocess
from multiprocessing import Pool
from itertools import repeat
from termcolor import colored
from datetime import datetime
import pandas as pd
from itertools import combinations
import math

#------------ IDRs for homer and macss
#-- expr: /home/dir1/exp1_rep1,/home/dir1/exp1_rep2;/home/dir1/exp2_rep1,/home/dir1/exp2_rep2,/home/dir1/exp2_rep3
#-- ctrl: /home/dir1/ctrl1_rep1,/home/dir1/ctrl1_rep2;/home/dir1/ctrl2_rep1,/home/dir1/ctrl2_rep2,/home/dir1/ctrl2_rep3
#-- names: exp1;expr2
#--- idrControl in run main
#-- Here I am using the script of Karmer: https://github.com/karmel/homer-idr. Mention it.

def run_cmd_file (mycmd,f,outputdir):
    print(colored(mycmd, 'green', attrs=['bold']))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        with open(f,'w') as file:
            logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
            process = subprocess.Popen(mycmd,
                                       stdout=file, 
                                       stderr=logfile)
            stdout, stderr = process.communicate()
            stdout, stderr 
    if process.returncode != 0:
        err_msg = ["error in the code, check log file:"]
        raise Exception(err_msg)

def macs2idrCommand (mysamples,peakList,inputFiletype,outputFile,IdrThreshold):
    c=['idr',
       '--samples'] + mysamples + ['--peak-list', peakList,
                                   '--input-file-type', inputFiletype,
                                   '--output-file', outputFile,
                                   '--rank', 'signal.value',
                                   '--soft-idr-threshold', str(IdrThreshold),
                                   '--plot']
    return(c)

        
def macs2PeakCalling (Style,CtrlInclude,SpikeInclude,Expr,Ctrl,OutDir,Prefix,Pvalue,effectiveGenomeSize,myratio):
    #---  
    if Style == "broad" or Style=="histone":
        if CtrlInclude == "True" and SpikeInclude == "False":
            c=['macs2',
               'callpeak',
               '-t', Expr,
               '-c', Ctrl,
               '-f', 'BAM',
               '-g', str(effectiveGenomeSize),
               '-n', Prefix,
               '--outdir', OutDir,
               '--keep-dup','all',
               '--seed','786',
               '--pvalue', Pvalue,
               '--broad']
        elif CtrlInclude == "False":
            c=['macs2',
               'callpeak',
               '-t', Expr,
               '-f', 'BAM',
               '-g', str(effectiveGenomeSize),
               '-n', Prefix,
               '--outdir', OutDir,
               '--keep-dup','all',
               '--seed','786',
               '--pvalue', Pvalue,
               '--broad']
        elif CtrlInclude == "True" and SpikeInclude == "True":
            c=['macs2',
               'callpeak',
               '-t', Expr,
               '-c', Ctrl,
               '-f', 'BAM',
               '-g', str(effectiveGenomeSize),
               '-n', Prefix,
               '--outdir', OutDir,
               '--keep-dup','all',
               '--seed','786',
               '--pvalue', Pvalue,
               '--broad',
               '--ratio', str(myratio)
              ]
            
    if Style == "narrow" or Style=="factor":
        if CtrlInclude == "True" and SpikeInclude == "False":
            c=['macs2',
               'callpeak',
               '-t', Expr,
               '-c', Ctrl,
               '-f', 'BAM',
               '-g', str(effectiveGenomeSize),
               '-n', Prefix,
               '--outdir', OutDir,
               '--keep-dup','all',
               '--seed','786',
               '--pvalue', Pvalue]
        elif CtrlInclude == "False":
            c=['macs2',
               'callpeak',
               '-t', Expr,
               '-f', 'BAM',
               '-g', str(effectiveGenomeSize),
               '-n', Prefix,
               '--outdir', OutDir,
               '--keep-dup','all',
               '--seed','786',
               '--pvalue', Pvalue]
        elif CtrlInclude == "True" and SpikeInclude == "True":
            c=['macs2',
               'callpeak',
               '-t', Expr,
               '-c', Ctrl,
               '-f', 'BAM',
               '-g', str(effectiveGenomeSize),
               '-n', Prefix,
               '--outdir', OutDir,
               '--keep-dup','all',
               '--seed','786',
               '--pvalue', Pvalue,
               '--ratio', str(myratio)
              ]
    return(c)
            
def macs2SpikeInratio (idrExprSpike, idrCtrlSpike, n_expr, n_ctrl, outputdir, idrName, threads, libraryType):
    idrExprSpikes=idrExprSpike.split(',')
    idrCtrlSpikes=idrCtrlSpike.split(',')
    #------
    if len(idrExprSpikes) != len(n_expr) or len(idrCtrlSpikes) != len(n_ctrl):
        print ( colored ("Number of the SpikeIn files are not equal",
                         "green",
                         attrs = ["bold"]
                        )
              )
        exit()
    #------
    replicate = 0
    n_ctrlSpike = []
    for ctr in idrCtrlSpikes:
        replicate += 1
        source = ctr

        if not os.path.exists(source):
            print(colored(source + ": file does not exist",
                          'green',
                          attrs = ["bold"]
                         )
                 )
            bam_present+=1

        dst = outputdir+'/idr_macs2/bamfiles/' + idrName + "_ctrlSpike_replicate" + str(replicate)
        try:
            if not os.path.exists(dst):
                os.symlink(source, dst)
        except OSError:
            pass
        n_ctrlSpike.append(dst)
    #------
    replicate = 0
    n_exprSpike = []
    for expr in idrExprSpikes:
        replicate += 1
        source = expr

        if not os.path.exists(source):
            print(colored(source + ": file does not exist",
                          'green',
                          attrs = ["bold"]
                         )
                 )
            bam_present+=1

        dst = outputdir+'/idr_macs2/bamfiles/' + idrName + "_exprSpike_replicate" + str(replicate)
        try:
            if not os.path.exists(dst):
                os.symlink(source, dst)
        except OSError:
            pass
        n_exprSpike.append(dst)

    #------
    totalCmd     = []
    totalCmdfile = []
    totalFiles = n_exprSpike + n_ctrlSpike + n_expr + n_ctrl
    for t in totalFiles:
        f=t+'.Flagstats.txt'
        c=['samtools','flagstat',t]

        totalCmdfile.append(f)
        totalCmd.append(c)
            
    with Pool (threads) as p:
        p.starmap(run_cmd_file,
              zip(totalCmd,
                  totalCmdfile,
                  repeat(outputdir)
                 )
                 )
    #------
    oValFiles   = []
    for f in n_exprSpike:
        oValFiles.append(f+".Flagstats.txt")
    c_exprSpike = initPeakCalling.spike_normalization2(oValFiles, libraryType)

    oValFiles   = []
    for f in n_ctrlSpike:
        oValFiles.append(f+".Flagstats.txt")
    c_ctrlSpike = initPeakCalling.spike_normalization2(oValFiles, libraryType)
    
    oValFiles   = []
    for f in n_expr:
        oValFiles.append(f+".Flagstats.txt")
    c_expr = initPeakCalling.spike_normalization2(oValFiles, libraryType)
    
    oValFiles   = []
    for f in n_ctrl:
        oValFiles.append(f+".Flagstats.txt")
    c_ctrl = initPeakCalling.spike_normalization2(oValFiles, libraryType)
    
    #------
    c_totalCtrl = sum(c_ctrl)
    c_totalCtrlSpike = sum(c_ctrlSpike)
    
    #------
    myratios=[]
    
    for i in range(0,len(n_expr)):
        myratios.append((c_exprSpike[i]/c_expr[i])/(c_totalCtrlSpike/c_totalCtrl))

    #------
    TotalExprRatio=(
        (sum(c_exprSpike)/sum(c_expr))/(c_totalCtrlSpike/c_totalCtrl)
    )
    
    #------
    return(myratios, TotalExprRatio)



def idr_homer (outputdir,idrExpr,idrCtrl,idrName,idrControl,idrStyle,idrOutput,blackListedRegions,threads):
    
    idrScript=psource.resource_filename(__name__, "otherScripts/idr/run_idr.py")
    totalCmd=[]
    totalCmd2=[]
    
    dirs=[outputdir+'/idr_homer/',
          outputdir+'/idr_homer/TagDir',
          outputdir+'/idr_homer/peaks/replicates/',
          outputdir+'/idr_homer/peaks/pooled',
          outputdir+'/idr_homer/peaks/pseudoreps',
          outputdir+'/idr_homer/peaks/pooled-pseudoreps',
          outputdir+'/idr_homer/pseudoreps/individual',
          outputdir+'/idr_homer/pseudoreps/pooled',
          outputdir+'/idr_macs2/'+idrName
         ]
          
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    cmd_rs=['findPeaks']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.',
                          'green', attrs=['bold']))
            exit()
    
    
    #--- link the tagdirectories to the idr/TagDir folder
    #_____________________________________________________
    replicate = 0
    n_expr = []
    tag_present=0
    for exp in idrExpr.split(','):
        replicate += 1
        source = exp
        if not os.path.exists(source):
            print(colored(source + ": file does not exist",
                          'green',
                          attrs = ["bold"]
                         )
                 )
            tag_present+=1
        dst = outputdir+'/idr_homer/TagDir/' + idrName + "_expr_replicate" + str(replicate)
        try:
            if not os.path.exists(dst):
                os.symlink(source, dst)
        except OSError:
            pass
        n_expr.append(dst)

    if idrControl=="True":
        replicate = 0
        n_ctrl = []
        for ctr in idrCtrl.split(','):
            replicate += 1
            source = ctr
            if not os.path.exists(source):
                print(colored(source + ": file does not exist",
                              'green',
                              attrs = ["bold"]
                             )
                     )
                tag_present+=1
            dst = outputdir+'/idr_homer/TagDir/' + idrName + "_ctrl_replicate" + str(replicate)
            try:
                if not os.path.exists(dst):
                    os.symlink(source, dst)
            except OSError:
                pass
            n_ctrl.append(dst)
    
    if tag_present > 0:
        exit()
        
    #-- Creating pooled directory
    c=['makeTagDirectory',
       outputdir+'/idr_homer/TagDir/TotalExperiment',
       "-d"]+n_expr
    
    totalCmd.append(c)
    
    if idrControl == "True":
        c=['makeTagDirectory',
           outputdir+'/idr_homer/TagDir/TotalControl',
           "-d"]+n_ctrl

        totalCmd.append(c)
        
    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )


    #-- Calling peaks
    totalCmd=[]
    replicate = 0
    for i in range(0,len(n_expr)):
        replicate+=1
        if idrControl == "True":
            c=['findPeaks',
               n_expr[i],
               '-i', outputdir+'/idr_homer/TagDir/TotalControl',
               '-o', outputdir + '/idr_homer/peaks/replicates/' + idrName + '_replicate' + str(replicate) + '.peaks.txt',
               '-P','.1',
               '-LP','.1',
               '-C', '0', # version 3 
               '-poisson', '.1',
               '-style', idrStyle]
            
        else:
             c=['findPeaks',
               n_expr[i],
               '-o', outputdir + '/idr_homer/peaks/replicates/' + idrName + '_replicate' + str(replicate) + '.peaks.txt',
               '-P','.1',
               '-LP','.1',
               '-C', '0', # version 3 
               '-poisson', '.1',
               '-style', idrStyle]
        
        totalCmd.append(c)
            
    if idrControl == "True":
        c=['findPeaks',
           outputdir+'/idr_homer/TagDir/TotalExperiment',
           '-i', outputdir+'/idr_homer/TagDir/TotalControl',
           '-o', outputdir + '/idr_homer/peaks/pooled/TotalExperiment.peaks.txt',
           '-P','.1',
           '-LP','.1',
           '-C', '0', # version 3 
           '-poisson', '.1',
           '-style', idrStyle]
        totalCmd.append(c)
    else:
        c=['findPeaks',
           outputdir+'/idr_homer/TagDir/TotalExperiment',
           '-o', outputdir + '/idr_homer/peaks/pooled/TotalExperiment.peaks.txt',
           '-P','.1',
           '-LP','.1',
           '-C', '0', # version 3 
           '-poisson', '.1',
           '-style', idrStyle]
        totalCmd.append(c)        
    
    #-- creating pseudoreplicates

    c=['python',
       idrScript,
       'pseudoreplicate',
       '-d'] + n_expr + ['-o', outputdir+'/idr_homer/pseudoreps/individual']
    totalCmd.append(c)
    
    c=['python',
       idrScript,
       'pseudoreplicate',
       '-d',
       outputdir+'/idr_homer/TagDir/TotalExperiment',
       '-o', outputdir+'/idr_homer/pseudoreps/pooled']
    totalCmd.append(c)
    
    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )

    #--- calling peaks in psuedo-replicates
    totalCmd=[]
    print ("_______________________________  nizam")
    files=glob.glob(outputdir+'/idr_homer/pseudoreps/individual/*')
    if idrControl=="True":
        for f in files:
            c=['findPeaks', f, 
               '-P', '.1',
               '-LP', '.1', 
               '-C', '0', # version 3 
               '-poisson', '.1',
               '-style', idrStyle,
               '-i', outputdir+'/idr_homer/TagDir/TotalControl',
               '-o', f+'_peaks.txt'
              ]
            totalCmd.append(c)
            
            c=['mv',
               f+'_peaks.txt',
               outputdir+'/idr_homer/peaks/pseudoreps'
              ]
            totalCmd2.append(c)
        
    else:
         for f in files:
            c=['findPeaks', f,
               '-P', '.1',
               '-LP', '.1',
               '-C', '0', # version 3 
               '-poisson', '.1',
               '-style', idrStyle,
               '-o', f+'_peaks.txt'
              ]

            totalCmd.append(c)       
        
            c=['mv',
               f+'_peaks.txt',
               outputdir+'/idr_homer/peaks/pseudoreps'
              ]
            totalCmd2.append(c)
        
    files=glob.glob(outputdir+'/idr_homer/pseudoreps/pooled/*')
    if idrControl == "True":
        for f in files:
            c=['findPeaks', f, 
               '-P', '.1',
               '-LP', '.1', 
               '-C', '0', # version 3 
               '-poisson', '.1',
               '-style', idrStyle,
               '-i', outputdir+'/idr_homer/TagDir/TotalControl',
               '-o', f+'_peaks.txt'
              ]
            totalCmd.append(c)

            c=['mv',
               f+'_peaks.txt',
               outputdir+'/idr_homer/peaks/pooled-pseudoreps'
              ]    
            totalCmd2.append(c)
    else:
        for f in files:
            c=['findPeaks', f,
               '-P', '.1',
               '-LP', '.1', 
               '-C', '0', # version 3 
               '-poisson', '.1',
               '-style', idrStyle,
               '-o', f+'_peaks.txt'
              ]
            totalCmd.append(c)

            c=['mv',
               f+'_peaks.txt',
               outputdir+'/idr_homer/peaks/pooled-pseudoreps'
              ]    
            totalCmd2.append(c)

        
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

    
    #--- calling IDRs
    p = glob.glob(outputdir+'/idr_homer/peaks/replicates/'+idrName+'*')
    if len(p) == 0:
        print (colored("Peaks are not present in idr_homer/peaks/replicates/ folder for "+idrName,
                       "green",
                       attrs = ["bold"]
                      )
              )
        
    pr = glob.glob(outputdir+'/idr_homer/peaks/pseudoreps/'+idrName+'*')
    if len(pr) == 0:
        print (colored("Peaks are not present in idr_homer/peaks/pseudoreps/ folder for "+idrName,
                       "green",
                       attrs = ["bold"]
                      )
              )
        
    ppr = glob.glob(outputdir+'/idr_homer/peaks/pooled-pseudoreps/'+"TotalExperiment"+'*')
    if len(ppr) == 0:
        print (colored("Peaks are not present in idr_homer/peaks/pooled-pseudoreps/ folder for "+idrName,
                       "green",
                       attrs = ["bold"]
                      )
              )
    pooled = glob.glob(outputdir+'/idr_homer/peaks/pooled/'+"TotalExperiment"+'*')
    if len(pooled) == 0:
        print (colored("Peaks are not present in idr_homer/peaks/pooled/ folder for "+idrName,
                       "green",
                       attrs = ["bold"]
                      )
              )
        
    c=['python', idrScript, 'idr', '-p'] + p + ['-pr'] + pr + ['-ppr'] + ppr + ['--pooled_peaks'] + pooled + ['-o', outputdir + "/idr_homer/" + idrOutput]

    universal.run_cmd(c,outputdir)

    #------- TotalExperiment.peaks-top-set.txt
   
    pd.read_csv(outputdir + "/idr_homer/" + idrOutput + '/TotalExperiment.peaks-top-set.txt',
                sep="\t").iloc[:,[1,2,3,0,7,4]].to_csv(outputdir + "/idr_homer/" + idrOutput + '/TotalExperiment.peaks-top-set.bed',
                                                       sep="\t",header=None,index=None)

    cmd=['intersectBed', 
         '-a', outputdir + "/idr_homer/" + idrOutput + '/TotalExperiment.peaks-top-set.bed', 
         '-b', blackListedRegions, 
         '-v']

    file=outputdir + "/idr_homer/" + idrOutput + '/TotalExperiment.peaks-top-set.removed.bed'
    with open(file,'w') as f:
        universal.run_cmd_file(cmd,f,outputdir)

    #-------
    cmd=['grep', '^chr', outputdir + "/idr_homer/" + idrOutput + '/TotalExperiment.peaks-top-set.removed.bed'] 
    file=outputdir + "/idr_homer/" + idrOutput + '/TotalExperiment.peaks-top-set.Clean.bed'
    with open(file,'w') as f:
        universal.run_cmd_file(cmd,f,outputdir)
    
#---------------------------------------------
#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

def idr_macs2 (outputdir,idrSpike,idrExpr,idrCtrl,idrExprSpike,idrCtrlSpike,idrName,idrControl,idrStyle,idrOutput,blackListedRegions,threads,effectiveGenomeSize,libraryType):
    
    idrScript=psource.resource_filename(__name__, "otherScripts/idr_macs2/bin/idr")
    totalCmd=[]
    totalCmd2=[]

    dirs=[outputdir+'/idr_macs2/',
          outputdir+'/idr_macs2/bamfiles/',
          outputdir+'/idr_macs2/peaks/replicates/',
          outputdir+'/idr_macs2/peaks/pooled',
          outputdir+'/idr_macs2/peaks/pseudoreps',
          outputdir+'/idr_macs2/peaks/pooled-pseudoreps',
          outputdir+'/idr_macs2/pseudoreps/individual',
          outputdir+'/idr_macs2/pseudoreps/pooled',
          outputdir+'/idr_macs2/'+idrName 
         ]

    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    cmd_rs=['macs2','samtools']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.',
                          'green', attrs=['bold']))
            exit() 

            
    #--- link the bamfiles
    #)))))))))))))))))))))))))))))))))))))))))))))))))
    replicate = 0
    n_expr = []
    bam_present = 0
    for exp in idrExpr.split(','):
        replicate += 1
        source = exp
        
        if not os.path.exists(source):
            print(colored(source + ": file does not exist",
                          'green',
                          attrs = ["bold"]
                         )
                 )
            bam_present+=1
        
        dst = outputdir+'/idr_macs2/bamfiles/' + idrName + "_expr_replicate" + str(replicate)
        try:
            if not os.path.exists(dst):
                os.symlink(source, dst)
        except OSError:
            pass
        n_expr.append(dst)

    if idrControl=="True":
        replicate = 0
        n_ctrl = []
        for ctr in idrCtrl.split(','):
            replicate += 1
            source = ctr
            
            if not os.path.exists(source):
                print(colored(source + ": file does not exist",
                              'green',
                              attrs = ["bold"]
                             )
                     )
                bam_present+=1

            dst = outputdir+'/idr_macs2/bamfiles/' + idrName + "_ctrl_replicate" + str(replicate)
            try:
                if not os.path.exists(dst):
                    os.symlink(source, dst)
            except OSError:
                pass
            n_ctrl.append(dst)
            
    if bam_present > 0:
        exit()
        
    #-- Creating pooled bamfiles
    #))))))))))))))))))))))))))))
    c=['samtools',
       'merge',
       '-f',
       '-O','BAM',
       '-@',str(threads),
       outputdir+'/idr_macs2/bamfiles/TotalExperiment']+n_expr
    
    totalCmd.append(c)
    
    if idrControl == "True":
        c=['samtools',
           'merge',
           '-f',
           '-O','BAM',
           '-@',str(threads),
           outputdir+'/idr_macs2/bamfiles/TotalControl']+n_ctrl

        totalCmd.append(c)
        

    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )

    #-- identifying the spikeIn ratio for calling peaks in the IDRs
    if idrControl == "True":
        if idrSpike == "True":
            myratios, TotalExprRatio = macs2SpikeInratio (idrExprSpike, idrCtrlSpike, n_expr, n_ctrl, outputdir, idrName, threads, libraryType)
        else:
            TotalExprRatio = "NA"
            myratios = ["NA"] * len(n_expr)
            
    #-- Calling peaks
    #)))))))))))))))))))))))))))))))))))))))))))))))))
    totalCmd=[]
    inExpr = outputdir+'/idr_macs2/bamfiles/TotalExperiment'
    inCtrl  = outputdir+'/idr_macs2/bamfiles/TotalControl'
    oDir   = outputdir + '/idr_macs2/peaks/pooled'
    Prefix = 'TotalExperiment'
    
    c=macs2PeakCalling (idrStyle, idrControl, idrSpike, inExpr, inCtrl, oDir, Prefix, '0.2', effectiveGenomeSize, TotalExprRatio)
    totalCmd.append(c)
    #--
    replicate = 0
    for i in range(0,len(n_expr)):
        replicate+=1
        inExpr = n_expr[i]
        inCtrl = outputdir+'/idr_macs2/bamfiles/TotalControl'
        oDir   = outputdir + '/idr_macs2/peaks/replicates/'
        Prefix = idrName + '_replicate' + str(replicate)
        c=macs2PeakCalling (idrStyle, idrControl, idrSpike, inExpr, inCtrl, oDir, Prefix, '0.2', effectiveGenomeSize, myratios[i])
        totalCmd.append(c)
                
    #-- creating pseudoreplicates
    #)))))))))))))))))))))))))))))))))))))))))))))))))
    cmd=['sambamba', 
         'view', 
         '-h', 
         '-t', str(threads), 
         '-s', str(0.5),
         '-f', 'bam', 
         '--subsampling-seed' + '=' + '786',
         outputdir+'/idr_macs2/bamfiles/TotalExperiment',
         '-o', outputdir+'/idr_macs2/pseudoreps/pooled/TotalExperiment-Pseudorep1'
        ]

    totalCmd.append(cmd)
    replicate = 0
    for i in range(0,len(n_expr)):
        replicate+=1    
        inFile = outputdir+'/idr_macs2/bamfiles/' + idrName + "_expr_replicate" + str(replicate)
        outFile= outputdir+'/idr_macs2/pseudoreps/individual/' + idrName + "_expr_replicate" + str(replicate) + '-Pseudorep1'
        cmd=['sambamba', 
             'view', 
             '-h', 
             '-t', str(threads), 
             '-s', str(0.5),
             '-f', 'bam', 
             '--subsampling-seed' + '=' + '786',
             inFile,
             '-o', outFile
            ]
        
        totalCmd.append(cmd)

    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )
    #-- 
    c1 = ['samtools',
          'view',
          outputdir+'/idr_macs2/pseudoreps/pooled/TotalExperiment-Pseudorep1'
         ]
    c2 = ['cut','-f1']
    with open(outputdir+'/idr_macs2/total.temp','w') as oFile:
        p1=subprocess.Popen(c1, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        p2=subprocess.Popen(c2, 
                            stdin=p1.stdout, 
                            stdout=oFile, 
                            stderr = subprocess.PIPE)
        stdout, stderr = p2.communicate()
        
    c1=['samtools',
        'view',
        '-h',
        outputdir+'/idr_macs2/bamfiles/TotalExperiment']
    
    c2=['fgrep',
        '-v', 
        '-f', 
        outputdir+'/idr_macs2/total.temp']
    
    c3=['samtools',
        'view',
        '-h', 
        '-',
        '-O','BAM',
        '-o',outputdir+'/idr_macs2/pseudoreps/pooled/TotalExperiment-Pseudorep2']
    
    p1=subprocess.Popen(c1, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    p2=subprocess.Popen(c2,
                        stdin  = p1.stdout,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE)
    p3=subprocess.Popen(c3, 
                        stdin  = p2.stdout,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE)
    stdout, stderr = p3.communicate()
        
    replicate = 0
    for i in range(0,len(n_expr)):
        replicate+=1    
        inFile = outputdir+'/idr_macs2/bamfiles/' + idrName + "_expr_replicate" + str(replicate)
        outFile= outputdir+'/idr_macs2/pseudoreps/individual/' + idrName + "_expr_replicate" + str(replicate) + '-Pseudorep2'
        c1 = ['samtools',
              'view',
              outputdir+'/idr_macs2/pseudoreps/individual/' + idrName + "_expr_replicate" + str(replicate) + '-Pseudorep1'
             ]
        c2 = ['cut','-f1']
        print(c1)
        print(c2)
        with open(outputdir+'/idr_macs2/total.temp','w') as oFile:
            p1=subprocess.Popen(c1, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            p2=subprocess.Popen(c2, 
                                stdin=p1.stdout, 
                                stdout=oFile, 
                                stderr = subprocess.PIPE)
            stdout, stderr = p2.communicate()

        c1=['samtools',
            'view',
            '-h',
            inFile]

        c2=['fgrep',
            '-v', 
            '-f', 
            outputdir+'/idr_macs2/total.temp']

        c3=['samtools',
            'view',
            '-h', 
            '-',
            '-O','BAM',
            '-o',outFile]

        p1=subprocess.Popen(c1, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        p2=subprocess.Popen(c2,
                            stdin  = p1.stdout,
                            stdout = subprocess.PIPE)
        p3=subprocess.Popen(c3, 
                            stdin  = p2.stdout,
                            stdout = subprocess.PIPE)
        stdout, stderr = p3.communicate()
        
        c=['rm',outputdir+'/idr_macs2/total.temp']
        universal.run_cmd(c,outputdir)

    #--- calling peaks in psuedo-replicates
    #)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
    totalCmd = []
    #--
    pseudos = ['-Pseudorep1','-Pseudorep2']
    for pseudo in pseudos:
        inExpr = outputdir+'/idr_macs2/pseudoreps/pooled/TotalExperiment' + pseudo
        inCtrl = outputdir+'/idr_macs2/bamfiles/TotalControl'
        oDir   = outputdir+'/idr_macs2/peaks/pooled-pseudoreps/'
        Prefix = 'TotalExperiment' + pseudo
    
        c=macs2PeakCalling (idrStyle, idrControl, idrSpike, inExpr, inCtrl, oDir, Prefix, '0.2', effectiveGenomeSize, TotalExprRatio)
        totalCmd.append(c)

    #----
    for pseudo in pseudos:
        replicate = 0
        for i in range(0,len(n_expr)):
            replicate+=1
            inExpr = outputdir+'/idr_macs2/pseudoreps/individual/' + idrName + "_expr_replicate" + str(replicate) + pseudo
            inCtrl = outputdir+'/idr_macs2/bamfiles/TotalControl'
            oDir   = outputdir + '/idr_macs2/peaks/pseudoreps'
            Prefix = idrName + '_replicate' + str(replicate) + pseudo
            c=macs2PeakCalling (idrStyle, idrControl, idrSpike, inExpr, inCtrl, oDir, Prefix, '0.2', effectiveGenomeSize, myratios[i])
            totalCmd.append(c)
        
    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
              zip(totalCmd,
                  repeat(outputdir)
                 )
                 )
    #--- calling IDRs
    #)))))))))))))))))
    IdrThreshold = 0.05
    IdrThresholdLog = -math.log(IdrThreshold,10)
    
    totalCmd = []
    if idrStyle == "narrow" or idrStyle=="factor":
        extensionFile = "_peaks.narrowPeak"
        inputFiletype = "narrowPeak"
    elif idrStyle == "broad" or idrStyle == "histone":
        extensionFile = "_peaks.broadPeak"
        inputFiletype = "broadPeak"
        
    mysamples = []
    
    if len(n_expr) == 2:
        #-- true replicates
        #))))))))))))))))))
        replicate = 0
        for i in range(0,len(n_expr)):
            replicate += 1
            mysamples.append(outputdir + '/idr_macs2/peaks/replicates/' + idrName + '_replicate' + str(replicate) + extensionFile)
            
        peakList = outputdir + '/idr_macs2/peaks/pooled/' + 'TotalExperiment' + extensionFile
        outputFile = outputdir + '/idr_macs2/' + idrName + '/' + 'trueRep' + extensionFile
        c = macs2idrCommand (mysamples,peakList,inputFiletype,outputFile,IdrThreshold)
        totalCmd.append(c)
  
        #-- pooled replicates
        #)))))))))))))))))))))
        mysamples=[
            outputdir + '/idr_macs2/peaks/pooled-pseudoreps/TotalExperiment-Pseudorep1' + extensionFile,
            outputdir + '/idr_macs2/peaks/pooled-pseudoreps/TotalExperiment-Pseudorep2' + extensionFile
        ]
        peakList = outputdir + '/idr_macs2/peaks/pooled/' + 'TotalExperiment' + extensionFile
        outputFile = outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoPooledRep' + extensionFile
        c= macs2idrCommand (mysamples,peakList,inputFiletype,outputFile,IdrThreshold)
        totalCmd.append(c)

        #-- pseudoreplicates of replicate 1
        #)))))))))))))))))))))
        mysamples=[
            outputdir + '/idr_macs2/peaks/pseudoreps/'+idrName+'_replicate1-Pseudorep1' + extensionFile,
            outputdir + '/idr_macs2/peaks/pseudoreps/'+idrName+'_replicate1-Pseudorep2' + extensionFile
        ]
        peakList = outputdir + '/idr_macs2/peaks/replicates/' +idrName+'_replicate1' + extensionFile
        
        outputFile = outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoRep1' + extensionFile
        c=macs2idrCommand (mysamples,peakList,inputFiletype,outputFile,IdrThreshold)
        totalCmd.append(c)
        
        #--- pseudoreplicates of replicate 2
        #)))))))))))))))))))))
        mysamples=[
            outputdir + '/idr_macs2/peaks/pseudoreps/'+idrName+'_replicate2-Pseudorep1' + extensionFile,
            outputdir + '/idr_macs2/peaks/pseudoreps/'+idrName+'_replicate2-Pseudorep2' + extensionFile
        ]
        peakList = outputdir + '/idr_macs2/peaks/replicates/' +idrName+'_replicate2' + extensionFile
        
        outputFile = outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoRep2' + extensionFile
        c=macs2idrCommand (mysamples,peakList,inputFiletype,outputFile,IdrThreshold)
        totalCmd.append(c)

        
        #--------
        with Pool (threads) as p:
            p.starmap(universal.run_cmd,
                  zip(totalCmd,
                      repeat(outputdir)
                     )
                     )
        #--------
        x = pd.read_csv (outputdir + '/idr_macs2/' + idrName + '/' + 'trueRep' + extensionFile,sep="\t", header=None)
        x = x[x.iloc[:,11]>IdrThresholdLog]
        x.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'trueRep' + extensionFile + "-IDRfilter",
                 header = None,
                 index  = None
                )
        Nt = x.shape[0]
        
        x.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'trueRep' + extensionFile + "-conservativeSets.txt",
                 header = None,
                 index  = None
                )      
        
        y = pd.read_csv (outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoPooledRep' + extensionFile, sep="\t", header=None)
        y = y[y.iloc[:,11]>IdrThresholdLog]
        y.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoPooledRep' + extensionFile + "-IDRfilter",
                 header = None,
                 index  = None
                )
        Np = y.shape[0]
        
        #-- rescue ratio
        #))))))))))))))))
        if min(Nt,Np) != 0:
            Rr = max(Nt,Np)/min(Nt,Np)
        else:
            print (colored("The number of the peaks is zero: (Nt = " + str(Nt) + ", Np = " + str(Np) + ")."+
                           " Exeriment quality is bad.",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()
        
        #-----
        
        z = pd.read_csv (outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoRep1' + extensionFile, sep="\t", header=None)
        z = z[z.iloc[:,11]>IdrThresholdLog]
        z.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoRep1' + extensionFile + "-IDRfilter",
                 header = None,
                 index  = None
                )
        N1 = z.shape[0]     
        
        
        
        k = pd.read_csv (outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoRep2' + extensionFile, sep="\t", header=None)
        k = k[k.iloc[:,11]>IdrThresholdLog]
        k.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'pseudoRep2' + extensionFile + "-IDRfilter",
                 header = None,
                 index  = None
                )
        N2 = k.shape[0]   
        
        #-- Self-consistency score
        #)))))))))))))))))))))))))
        if min(N1,N2) != 0:
            Sc = max(N1,N2)/min(N1,N2)
        else:
            print (colored("The number of the peaks is zero: (N1 = " + str(N1) + ", N2 = " + str(N2) + ")." +
                           " Experiment quality is bad.",
                           "green",
                           attrs = ["bold"]
                          )
                  )
            exit()
        #-- reporting about IDRs
        #))))))))))))))))))))))))
        if Sc > 2:
            print (colored("It seems that replicates for the "+
                           idrName + "are not true replicates because self consistency score > 2. Score is: " +  str(Sc),
                           "green",
                           attrs =["bold"]
                          )
                  )
        if Rr > 2 and Sc > 2:
            print (colored ("Both rescue ratio and self consistency are >2. Recommended to not use these IDR peaks.",
                            "green",
                            attrs = ["bold"]
                           )
                  )
        elif Rr > 2 or Sc > 2:
            print (colored ( "Reproducibility is in borderline. Its up to you if you want to use IDR peaks.",
                            "green",
                            attrs = ["bold"]
                           )
                  )
        #-- Optimal sets
        #)))))))))))))))
        
        if Nt > Np:
            x.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'trueRep' + extensionFile + "-optimalSets.txt",
                     header = None,
                     index  = None
                    )
        elif Nt < Np:
            y.to_csv(outputdir + '/idr_macs2/' + idrName + '/' + 'trueRep' + extensionFile + "-optimalSets.txt",
                     header = None,
                     index  = None
                    )
    else:
        print (colored ("The number of the replicates are greater than 2. Contact developer",
                        "red",
                        attrs = ["bold"]
                       )
              )
        

