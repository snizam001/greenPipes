#-- Quality of the homer TagDirectories 
#____________________________________
#if mode == 'quality_tagDirectories' or mode == 'all':
from multiprocessing import Pool
from itertools import repeat
import os
import subprocess
from datetime import datetime
from termcolor import colored
import pkg_resources as psource
import pandas as pd 
from greenPipe import universal

       
def qcTagD_fragmentLength (outputdir,Name,threads,qcSpike):
    dirs=outputdir+'/'+'Tagdirectories_qualities/'
    folders=[dirs]
    for folder in folders:
        if not os.path.exists(outputdir+'/'+folder):
            os.makedirs(outputdir+'/'+folder)
    #--------------
    if qcSpike != "False":
        if not (os.path.exists(outputdir + '/SpikeIn/' + Name + '_expr' + '.bam') or
                os.path.exists(outputdir + '/SpikeIn/' + Name + '_control' + '.bam')
               ):
            print(colored(outputdir + '/SpikeIn/' + Name + '_expr' + '.bam'+' or '+
                          outputdir + '/SpikeIn/' + Name + '_control' + '.bam' +
                          ':does not exit',
                          'green', attrs=['bold']))
            exit()
    if not (os.path.exists(outputdir + '/Bamfiles/' + Name + '_expr' + '.bam') or
            os.path.exists(outputdir + '/Bamfiles/' + Name + '_control' + '.bam')
           ):
        print(colored(outputdir + '/Bamfiles/' + Name + '_expr' + '.bam'+' or '+
                      outputdir + '/Bamfiles/' + Name + '_control' + '.bam' +
                      ':does not exit',
                      'green', attrs=['bold']))
        exit()
    #--------------
    for ctrl_expr in ['_expr','_control']:
        if qcSpike != "False":
            c=['samtools',
               'sort',
               '-n',
               '-O', 'BAM',
               '-@', str(threads),
               '-o', outputdir + '/SpikeIn/' + Name + ctrl_expr + '-sorted.bam',
               outputdir + '/SpikeIn/' + Name + ctrl_expr + '.bam']

            universal.run_cmd(c,outputdir)
        
        c=['samtools',
           'sort',
           '-n',
           '-O', 'BAM',
           '-@', str(threads),
           '-o', outputdir + '/Bamfiles/' + Name + ctrl_expr + '-sorted.bam',
           outputdir + '/Bamfiles/' + Name + ctrl_expr + '.bam']

        universal.run_cmd(c,outputdir)
        #--------------
        
        if qcSpike != "False":
            c=['bedtools',
               'bamtobed',
               '-bedpe',
               '-i', outputdir + '/SpikeIn/' + Name + ctrl_expr + '-sorted.bam']
            with open(dirs + Name + ctrl_expr  + '.spike.sorted.bed', "w+") as f:
                universal.run_cmd_file(c,f,outputdir)
                
        c=['bedtools',
           'bamtobed',
           '-bedpe',
           '-i', outputdir + '/Bamfiles/' + Name + ctrl_expr + '-sorted.bam']

        with open(dirs + Name + ctrl_expr  + '.sorted.bed', "w+") as f:
            universal.run_cmd_file(c,f,outputdir)
        #--------------
        if qcSpike != "False":
            seacrBed=pd.read_csv(dirs + Name + ctrl_expr  + '.spike.sorted.bed',
                                 sep='\t',header=None)
            seacrBed=seacrBed.iloc[:,[0,1,5]]
            seacrBed.columns=['Chr','Start','End']
            seacrBed=seacrBed[seacrBed['Chr'].str.contains("chr|Chr")]
            seacrBed=seacrBed[~seacrBed['Chr'].str.contains("_")]
            seacrBed.to_csv(dirs + '/' + Name + ctrl_expr  + '.spike.sorted.fragments.bed',
                            sep='\t',index=None,header=None)

            cmd=['sort','-k1,1','-k2,2n',dirs + '/' + Name + ctrl_expr  + '.spike.sorted.fragments.bed']
            with open(dirs + '/' + Name + ctrl_expr  + '.spike.sorted.fragments.sorted.bed', "w+") as f:
                universal.run_cmd_file(cmd,f,outputdir)

            cmd=['mv',dirs + '/' + Name + ctrl_expr  + '.spike.sorted.fragments.sorted.bed',
                 dirs + '/' + Name + ctrl_expr  + '.spike.sorted.fragments.bed']
            universal.run_cmd(cmd,outputdir)
        
        
        
        seacrBed=pd.read_csv(dirs + Name + ctrl_expr  + '.sorted.bed',
                             sep='\t',header=None)
        seacrBed=seacrBed.iloc[:,[0,1,5]]
        seacrBed.columns=['Chr','Start','End']
        seacrBed=seacrBed[seacrBed['Chr'].str.contains("chr|Chr")]
        seacrBed=seacrBed[~seacrBed['Chr'].str.contains("_")]
        seacrBed.to_csv(dirs + '/' + Name + ctrl_expr  + '.sorted.fragments.bed',
                        sep='\t',index=None,header=None)

        cmd=['sort','-k1,1','-k2,2n',dirs + '/' + Name + ctrl_expr  + '.sorted.fragments.bed']
        with open(dirs + '/' + Name + ctrl_expr  + '.sorted.fragments.sorted.bed', "w+") as f:
            universal.run_cmd_file(cmd,f,outputdir)

        cmd=['mv',dirs + '/' + Name + ctrl_expr  + '.sorted.fragments.sorted.bed',
             dirs + '/' + Name + ctrl_expr  + '.sorted.fragments.bed']
        universal.run_cmd(cmd,outputdir)

    if qcSpike == "False":    
        expr_fragmentLength=dirs + '/' + Name + '_expr.sorted.fragments.bed'
        ctrl_fragmentLength=dirs + '/' + Name + '_control.sorted.fragments.bed'
        return(expr_fragmentLength,ctrl_fragmentLength)
    else:
        expr_fragmentLength=dirs + '/' + Name + '_expr.sorted.fragments.bed'
        ctrl_fragmentLength=dirs + '/' + Name + '_control.sorted.fragments.bed'
        exprSpk_fragmentLength=dirs + '/' + Name + '_expr.spike.sorted.fragments.bed'
        ctrlSpk_fragmentLength=dirs + '/' + Name + '_control.spike.sorted.fragments.bed'
        return(expr_fragmentLength,ctrl_fragmentLength,exprSpk_fragmentLength,ctrlSpk_fragmentLength)
        
def qcTagD (libraryType,outputdir,Names,threads,qcTagR):
    totalCmd=[]
    bedpeCmd=[]
    cmd_rs=['samtools','bedtools',qcTagR] 
    #---  
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH', 
                          'green', attrs=['bold']))
            exit()

    dirs=[outputdir+'/'+'Tagdirectories_qualities/', outputdir + '/Tagdirectories/']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    for Name in Names:

        output_jpeg=dirs[0] + Name + '_expr' + '.jpeg'

        infile_expr=dirs[1] + Name + '_expr'
        infile_ctrl=dirs[1] + Name + '_control'

        expr_tagCounts=infile_expr+'/'+'tagCountDistribution.txt'
        ctrl_tagCounts=infile_ctrl+'/'+'tagCountDistribution.txt'
        expr_tagAutocorrelation=infile_expr+'/'+'tagAutocorrelation.txt'
        ctrl_tagAutocorrelation=infile_ctrl+'/'+'tagAutocorrelation.txt'

        tagFiles=[expr_tagCounts,ctrl_tagCounts,expr_tagAutocorrelation,ctrl_tagAutocorrelation]
        
        p_tagFile=0
        for tagFile in tagFiles:
            if not os.path.exists(tagFile):
                p_tagFile=p_tagFile+1
                print(colored(tagFile+' :does not exit',
                              'green', attrs=['bold']))
        if p_tagFile > 0:
            exit()            
        #--- length of the fragments in expreriments
        #____
        if libraryType == 'pair':
            expr_fragmentLength,ctrl_fragmentLength=qcTagD_fragmentLength(outputdir,Name,threads,"False")
            cmd=[qcTagR,
                 '-a',expr_tagCounts,
                 '-b',ctrl_tagCounts,
                 '-c',expr_tagAutocorrelation,
                 '-d',ctrl_tagAutocorrelation,
                 '-e',expr_fragmentLength,
                 '-f',ctrl_fragmentLength,
                 '-o',output_jpeg]
            totalCmd.append(cmd)
#            universal.run_cmd(cmd,outputdir)
        else:
            cmd=[qcTagR,
                 '-a',expr_tagCounts,
                 '-b',ctrl_tagCounts,
                 '-c',expr_tagAutocorrelation,
                 '-d',ctrl_tagAutocorrelation,
                 '-o',output_jpeg]
            totalCmd.append(cmd)
#            universal.run_cmd(cmd,outputdir)

        return(totalCmd)
