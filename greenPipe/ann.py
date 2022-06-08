#-- annotations
#_________________
import os
import subprocess
from datetime import datetime
from termcolor import colored
import glob
import pkg_resources as psource
from greenPipe import filterMotifSeq
from greenPipe import universal
from multiprocessing import Pool
from itertools import repeat

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

def UserAnn (outputdir,annpeakFiles,annFiles,annSize,annName,annPrefix,threads):
    annR=psource.resource_filename(__name__, "rscripts/ann.R")
    cmd_rs=[annR]
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.',
                          'green', 
                          attrs=['bold']
                         )
                 )
            exit()
    
    dirs=[outputdir+'/'+'UserAnnotation']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)


    myPeakFolders = []
    zz = ['Peaks', 'idr_homer', 'idr_macs2']
    if annpeakFiles=='None':
        for z in zz:
            if os.path.exists(outputdir + '/' + z + '/'):
                myPeakFolders = myPeakFolders + [outputdir + '/' + z + '/']
        annpeakFile = []
        for myPeakFolder in myPeakFolders:
            annpeakFile = annpeakFile + glob.glob(outputdir + '/' + myPeakFolder + '/' + '*.Clean.bed')
            print(colored("UserAnnotation mode is active. Using peaks from "+
                          outputdir + '/' + myPeakFolder + '/ folder and using following files: ',
                          'green',
                          attrs = ['bold']
                         )
                 )
        for kk in annpeakFile:
            print(kk)
        annPrefixes=[jj.split('/')[-1].replace('.Clean.bed','') for jj in annpeakFile]
    else:
        annpeakFile=annpeakFiles.split(',')
        print(colored("UserAnnotation mode is active. User provided peak files: ",
                      'green',
                      attrs = ['bold']
                     )
             )
        for kk in annpeakFile:
            print(kk)

        annPrefixes=annPrefix.split(',')
    if len(annPrefixes) != len(annpeakFile):
      print(colored("The lenght of the --annPrefix is not equal to the --annpeakFiles.",
                    'green',
                    attrs = ['bold']
                   )
           )
      exit()
    #-- checking if given file exists
    c_files=annpeakFiles.split(',')+annFiles.split(',')
    print(c_files)
    e_file=0
    for c_file in c_files:
        if not os.path.exists(c_file) or os.path.getsize(c_file)==0:
            e_file=e_file+1
            print(colored(c_file+": does not exist or empty",
                          'green',
                          attrs=['bold']
                         )
                 )
    if e_file > 0:
        exit()
    
    if annSize == 'None' or annFiles == "None" or annName == "None":
        print(colored("Either --annSize or --annFiles or --annName is missing",
                      'green',
                      attrs=['bold']
                     )
             )
        exit()

    totalCmd=[]

    for i in range(0,len(annpeakFile)):
        outFile=outputdir+'/UserAnnotation/'+annPrefixes[i]+'.txt'
        cmd=[annR,
             '-m', annSize,
             '-i', annpeakFile[i],
             '-d', annFiles,
             '-c', str(threads),
             '-n', annName,
             '-o', outFile]
        totalCmd.append(cmd)

    return(totalCmd)

def Ann (annpeakFiles, annPrefix, cGVersion, sFasta, outputdir, threads):
    cmd_rs=['annotatePeaks.pl', 'meme', 'homerTools', 'tomtom']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': Part of homer. It is not installed in your computer or not in the PATH.',
                          'green', 
                          attrs=['bold']
                         )
                 )
            exit()
    
    dirs=[
        
        outputdir+'/'+'Annotation', 
        outputdir+'/'+'GeneOntology',
        outputdir+'/'+'Motifs', 
        outputdir+'/'+'Motifs/Meme',
        outputdir+'/'+'Motifs/Homer'
        
    ]
    
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    if sFasta == 'NA':
        print(colored('--sFasta is missing',
                      'green',
                      attrs = ['bold']
                     )
             )
        exit()
        
    myPeakFolders = []
    zz = ['Peaks', 'idr_homer', 'idr_macs2']
    if annpeakFiles=='None':
        for z in zz:
            if os.path.exists(outputdir + '/' + z + '/'):
                myPeakFolders = myPeakFolders + [outputdir + '/' + z + '/']
            else:
                print ("Folder " + z + 'does not exist.')
        annpeakFile = []
        for myPeakFolder in myPeakFolders:
            annpeakFile = annpeakFile + glob.glob(outputdir + '/' + myPeakFolder + '/' + '*.Clean.bed')
            print(colored("annotation mode is active. Using peaks from "+
                          outputdir + '/' + myPeakFolder + '/ folder and using following files: ',
                          'green',
                          attrs = ['bold']
                         )
                 )
        for kk in annpeakFile:
            print(kk)
        annPrefixes=[jj.split('/')[-1].replace('.Clean.bed','') for jj in annpeakFile]
    else:
        annpeakFile=annpeakFiles.split(',')
        print(colored("annotation mode is active. User provided peak files: ",
                      'green',
                      attrs = ['bold']
                     )
             )
        for kk in annpeakFile:
            print(kk)

        annPrefixes=annPrefix.split(',')

    if len(annPrefixes) != len(annpeakFile):
      print(colored("The lenght of the --annPrefix is not equal to the --annpeakFiles.",
                    'green',
                    attrs = ['bold']
                   )
           )
      exit()
    #-- checking if given file exists
    c_files=annpeakFiles.split(',')
    print(c_files)
    e_file=0
    for c_file in c_files:
        if not os.path.exists(c_file) or os.path.getsize(c_file)==0:
            e_file=e_file+1
            print(colored(c_file+": does not exist or empty",
                          'green',
                          attrs=['bold']
                         )
                 )
    if e_file > 0:
        exit()
    #-- annotation and the gene ontologies
    #______________________________________________
    for i in range(0,len(annpeakFile)):
        outFile=outputdir+'/UserAnnotation/'+annPrefixes[i]+'.txt'
        cmd = [
            'annotatePeaks.pl',
            annpeakFile[i],
            cGVersion,
            '-go', outputdir + '/GeneOntology/' + annPrefixes[i] + '_GO',
            '-annStats', outputdir + '/Annotation/' + annPrefixes[i] + '-annStats.txt',
            '-cpu', str(threads)
        ]
        
        with open(outputdir + '/Annotation/' + annPrefixes[i] + '-annotation.txt', 'w') as anntotalCmdFile:
            universal.run_cmd_file(cmd,
                                   anntotalCmdFile,
                                   outputdir
                                  )
    #-- motifs through homer 
    #______________________________________________
    for i in range(0,len(annpeakFile)):
        filterMotifSeq.motifs_peaks (
            cGVersion,
            annpeakFile[i],
            outputdir + '/Motifs/Homer/' + annPrefixes[i] ,
            str(400),
            str(threads),
            outputdir
        )

    #-- motifs through meme 
    #______________________________________________
    totalCmd = []
    totalCmdFile = []
    for i in range(0,len(annpeakFile)):
        cmd = [
            'homerTools',
            'extract',
            annpeakFile[i],
            sFasta,
            '-fa']
        totalCmd.append(cmd)
        totalCmdFile.append(outputdir + '/Motifs/Meme/' + annPrefixes[i] + '.fa')
    with Pool (threads) as p:
        p.starmap(run_cmd_file,
              zip(totalCmd,
                  totalCmdFile,
                  repeat(outputdir)
                 )
                 )

    for i in range(0,len(annpeakFile)):
        cmd = [
            'meme',
            outputdir + '/Motifs/Meme/' + annPrefixes[i] + '.fa',
            '-oc',
            outputdir + '/Motifs/Meme/' + annPrefixes[i] + '_meme',
            '-dna',
            '-p', str(threads),
            '-nmotifs', '10'
        ]
        universal.run_cmd(cmd, outputdir)

    totalCmd = []
    target=psource.resource_filename(__name__, "data/JASPAR2018_CORE_non-redundant.meme")
    for i in range(0,len(annpeakFile)):
        cmd = [
        'tomtom',
        '-oc',outputdir + '/Motifs/Meme/' + annPrefixes[i] + '_tomtom',
        outputdir + '/Motifs/Meme/' + annPrefixes[i] + '_meme/meme.txt',
        target
        ]
        totalCmd.append(cmd)

    with Pool (threads) as p:
        p.starmap(universal.run_cmd,
                  zip(totalCmd,
                      repeat(outputdir)
                     )
                 )
