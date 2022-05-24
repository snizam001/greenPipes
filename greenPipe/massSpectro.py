#- ------------------- mass spectrometry or piggy back binding events 
import pkg_resources
import pkg_resources as psource
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

def run_cmd (c,outputdir):
    print(colored(c, 'green', attrs=['bold']))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        process = subprocess.Popen(c,
                                   stdout=subprocess.PIPE, 
                                   stderr=logfile)
        stdout, stderr = process.communicate()
        stdout, stderr
    if process.returncode != 0:
        err_msg = ["error in the code, check log file:"]
        raise Exception(err_msg)

def run_cmd_file (mycmd,f,outputdir):
    print(colored(mycmd, 'green', attrs=['bold']))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        process = subprocess.Popen(mycmd,
                                   stdout=f, 
                                   stderr=logfile)
        stdout, stderr = process.communicate()
        stdout, stderr
    if process.returncode != 0:
        err_msg = ["error in the code, check log file:"]
        raise Exception(err_msg)

def motifs_peaks (genomeversion,mypeak,out_dir,size,threads,outputdir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    jFile=psource.resource_filename(__name__, "data/Jasper2018_homer.txt")
    c=['findMotifsGenome.pl',
       mypeak,
       genomeversion,
       out_dir,
       '-size', str(size),
       '-nomotif',
       '-p', str(threads),
       '-mknown', jFile]
    run_cmd(c,outputdir)

def massHugo2Entrez (infile,Species):
    try:
        d=pd.read_csv(infile,sep='\t',header=None)
    except IOError:
        print(colored("massHugo2Entrez: " + infile+' does not exist',
                      'green', 
                      attrs=['bold']
                     )
             )
        exit()

    Species_mg=["human",
                "mouse",
                "rat",
                "fruitfly",
                "nematode",
                "zebrafish",
                "thale-cress",
                "frog",
                "pig"]
    
    if Species not in Species_mg:
        print(colored("massHugo2Entrez: " + 
                      Species+
                      ' is not correct. Possible options are (Case sensitive):',
                      'green', 
                      attrs=['bold']
                     ) 
             )
        print(Species_mg)
        exit()
        
    d.columns=['Gene']
    entrezIds=[]
    for j in range(0,d.shape[0]):
        if ";" in d.iloc[j,0]:
            print(d.iloc[j,0])
            print(colored('The file '+ infile +' contains ";" or any other special character.'+
                          ' Did you read the manual? it should not be present in the file', 
                          'green', 
                          attrs=['bold']
                         ) 
                 )
            exit()
        else:
            mgout=mg.query(d.iloc[j,0],
                           fields='entrezgene,taxid,symbol',
                           as_dataframe=True,
                           species=Species)
            if mgout.empty:
                print('Can not find entrez gene id for '+ 
                      d.iloc[j,0] +" for " + 
                      Species)
            else:
                entrezIds.append(mgout['entrezgene'].tolist()[0])
    print ('entrez ids for genes present in MassSpectrometry datasets: ')
    print (entrezIds)
    return(entrezIds)


def ambiguous_dnasequence (seq):
    ambig = {"R": ["A", "G"], 
             "V":["A", "C", "G"],
             "Y":["C","T"],
             "S":["G","C"],
             "W":["A","T"],
             "K":["G","T"],
             "M":["A","C"],
             "B":["C","G","T"],
             "D":["A","G","T"],
             "H":["A","C","T"],
             "N":["A","C","T","G"]}
    groups = itertools.groupby(seq, lambda char:char not in ambig)
    splits = []
    for b,group in groups:
        if b:
            splits.extend([[g] for g in group])
        else:
            for nuc in group:
                splits.append(ambig[nuc])
    return([''.join(p) for p in itertools.product(*splits)])


def filterMotifSequence (infile,
                         outfile1,
                         outfile2,
                         outfile3,
                         pwm,
                         sequence,
                         tempDir,
                         genomeversion,
                         size,
                         threads,
                         dist,
                         MassSpectro_Fasta,
                         outputdir):

    if not os.path.exists(infile):
        print(colored("filterMotifSequence: " + infile+' does not exist',
                      'green', 
                      attrs=['bold']
                     ) 
             )
    
    mycmd=['findMotifsGenome.pl', 
           infile, 
           genomeversion, 
           tempDir, 
           '-size', str(size), 
           '-nomotif', 
           '-p', str(threads), 
           '-find', pwm] 
    
    with open(outfile1, "w+") as f:
        run_cmd_file(mycmd,f,outputdir)
        
    pwmPeaks=pd.read_csv(outfile1,sep='\t')
    peakFile=pd.read_csv(infile,
                         sep='\t',
                         header=None)
    print ('Total number of the peaks in the ' + 
           infile + ' is :' + 
           str(peakFile.shape[0]))
    
    peakFile=peakFile[~peakFile.iloc[:,3].isin(pwmPeaks['PositionID'])]
    print ('Total number of the peaks after filtering peaks'+
           ' containing pwmProteinOfInterest motifs :' + 
           str(peakFile.shape[0]))
    
    seqs=ambiguous_dnasequence(sequence)
    peakWithSeq=[]
    
    for i in range(0,peakFile.shape[0]):
        ID=peakFile.iloc[i,0]
        Start=int((peakFile.iloc[i,1]+peakFile.iloc[i,2])/2)-dist
        End=int((peakFile.iloc[i,1]+peakFile.iloc[i,2])/2)+dist
        Name=peakFile.iloc[i,3]
        SeqPos=MassSpectro_Fasta[ID][Start:End].seq
        SeqNeg=SeqPos[::-1].translate(SeqPos[::-1].maketrans("ATGCatgc", "TACGtacg"))
        for seq in seqs:
            if seq in SeqPos or seq in SeqNeg:
                peakWithSeq.append(Name)

    peakFile=peakFile[~peakFile.iloc[:,3].isin(peakWithSeq)]
    print ('Total number of the peaks after filtering peaks containing'+
           ' pwmProteinOfInterest motifs + sequenceProteinOfInterest :' + 
           str(peakFile.shape[0]))

    peakFile.to_csv(outfile2,
                    header=None,
                    index=None,
                    sep='\t')
    
    motifs_peaks(genomeversion,
                 outfile2,
                 outfile3,
                 size,
                 threads,
                 outputdir)

    
def vennDiagram (outDir,infile,distance,Names):
    if len(infile) != len(Names):
        print(colored("Length of name and infile are matching",
                      'green', 
                      attrs=['bold']
                     ) 
             )
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        
    mycmd=['mergePeaks', 
           '-d', str(distance)] + infile + ['-venn', out_dir+'/'+'Comparison_venn', 
                                            '-matrix', out_dir+'/'+'Comparison_matrix']
    
    with open(out_dir+'/'+'Comparison.txt', "w+") as f:
        run_cmd_file(mycmd,f,outputdir)
        
    #---- cleaning of the files:1
    extensions=['Comparison_venn','Comparison.txt']
    
    for extension in extensions:
        venn = pd.read_csv(out_dir+extension,sep='\t',header=None)
        i=-1
        for Name in Names:
            i=i+1
            venn = venn.astype(str).replace(infile[i],Name,regex=True)
        venn.to_csv(outDir+"/"+Name+"."+extension,index = None, header=None,sep='\t')
    venn = pd.read_csv(outDir+"/"+Name+"."+'Comparison_venn',sep='\t',header=None)
    
    if (venn.shape[1]-2) < 5:
        vennR=psource.resource_filename(__name__, "rscripts/venn.R")
        mycmd=[vennR,
               'Comparison_venn',
               outDir+"/"+Name+"."+'Comparison_venn']
        run_cmd(mycmd)
    else:
        print(colored("Unable to plot Venn diagram as categories are >5",
                      'green', 
                      attrs=['bold']
                     ) 
             )

def jasper_motif_wrapper (Name,outFile):
    jFile=psource.resource_filename(__name__, "data/Jasper2018_homer.txt")
    on = 0;
    d = pd.read_csv(jFile,sep='\t',header=None)
    AT_freq=[]
    for i in range(0,d.shape[0]):
        if ">" in d.iloc[i,0]:
            if Name==d.iloc[i,1].split(':')[-1]:
                on = 1
            else:
                on = 0
        if on == 1:
            AT_freq.append(list(d.iloc[i,:]))
    pd.DataFrame(AT_freq).to_csv(outFile,sep="\t",header=None,index=None)
        
def intersect_mass_green (Jasper2IDs,mIn,mOut,pvalue):
    
    Common_geneIDs = []
    mInfo = pd.read_csv(mIn,sep='\t')
    mInfo = mInfo.loc[mInfo.iloc[:,4]<pvalue,]
    mEntrezIDs = []
    
    for mInfo_i in range(0,mInfo.shape[0]):
        mName = mInfo.iloc[mInfo_i,0]
        for Jasper2ID_i in range(0,Jasper2IDs.shape[0]):
            Jasper2IDs_name = Jasper2IDs.iloc[Jasper2ID_i,0]
            if mName == Jasper2IDs_name:
                mEntrezIDs.append(Jasper2IDs.iloc[Jasper2ID_i,2])
                
    mInfo['entrezIDs'] = mEntrezIDs
    
    Common_geneIDs = set(mInfo.loc[:,'entrezIDs']).intersection(set(entrezOfProteins))
    
    mInfo = mInfo[mInfo.entrezIDs.isin(Common_geneIDs)]
    
    mInfo.to_csv(mOut,
                 sep = '\t', 
                 index = None)
    
    return(mInfo)


def proof_of_piggyBack (inFile,
                        metaFile,
                        EntrezIDs,
                        HugoSymbol,
                        outDir,
                        outFile,
                        distanceOverlapPiggyBinding,
                        Name,
                        genomeversion,
                        tempDir,
                        size,
                        threads,
                        pwm,
                        outputdir):
    
    if not os.path.exists(metaFile):
        print(colored(metaFile+" does not exist",'green', attrs=['bold']) )
    else:
        metaEncode=pd.read_csv(metaFile,sep='\t')
    
    if metaEncode.shape == 0:
        print(colored(metaFile+" is empty",'green', attrs=['bold']) )
    
    #---
    jasper_motif_wrapper(HugoSymbol,pwm) 

    mycmd=['findMotifsGenome.pl', 
           inFile, 
           genomeversion, 
           tempDir, 
           '-size', str(size), 
           '-nomotif', 
           '-p', str(threads), 
           '-find', pwm] 
    
    with open(outFile, "w+") as f:
        run_cmd_file(mycmd,f,outputdir)
        
    pwmPeaks=pd.read_csv(outFile,sep='\t')
    
    peakFile=pd.read_csv(inFile,sep='\t',header=None)
    
    peakFile=peakFile[peakFile.iloc[:,3].isin(pwmPeaks['PositionID'])]
    
    peakFile.to_csv(outFile,sep="\t",header=None,index=None)
    
    #---
    
    selFile=metaEncode[metaEncode['EntrezIDs']==int(EntrezIDs)]
    print(colored("Go throught this file and judge if the output reliable or not: "+
                  inFile.replace('.txt','')+'_SelectedDataSetsEncode.txt',
                  'green', 
                  attrs=['bold']) )
    
    selFile.to_csv(inFile.replace('.txt','')+'_SelectedDataSetsEncode.txt',index=None,sep="\t")
    
    #---
    if selFile.shape[0] != 0:
        for i in range(0,selFile.shape[0]):
            ENCName=selFile.iloc[i,22].replace('-human','')
            ENCCellLines=selFile.iloc[i,10].replace(" ","_")
            selectedFile="/media/linux/HD/greeCUTRUN/Data/Encode_tf/" + selFile.iloc[i,0]+'.bed'    #------------ change here the package in a way that metadata automatically uploaded
            
            mycmd=['mergePeaks', 
                   '-d', str(distanceOverlapPiggyBinding)] + [selectedFile,outFile] + ['-venn', outDir+'/'+'Comparison']
            
            with open(outDir+'/'+Name+'vs'+ENCName+'-'+ENCCellLines+selFile.iloc[i,0]+'.Comparison.txt','w') as f:
                run_cmd_file(mycmd,f,outputdir)

            venn = pd.read_csv(outDir+'/'+'Comparison',sep='\t',header=None)
            venn = venn.astype(str).replace("/media/linux/HD/greeCUTRUN/Data/Encode_tf/" + selFile.iloc[i,0]+'.bed',ENCName,regex=True)
            venn = venn.astype(str).replace(outFile,Name,regex=True)
            venn.to_csv(outDir+'/'+'Comparison',index = None, header=None,sep='\t')
            
            
            vennR=psource.resource_filename(__name__, "rscripts/venn.R")
            mycmd=[vennR,
                   Name+'vs'+ENCName+'-'+ENCCellLines+selFile.iloc[i,0]+'.Comparison',
                   outDir+'/'+'Comparison']
            
            run_cmd(mycmd,outputdir)
            
            mycmd=['mv',
                   outDir+'/'+Name+'vs'+ENCName+'-'+ENCCellLines+selFile.iloc[i,0]+'.Comparison_venn.pdf',
                   outputdir+'/'+'MassSpectrometry/'+'VennDiagrams'+'/'+Name+'vs'+ENCName+'-'+ENCCellLines+selFile.iloc[i,0]+'.Comparison_venn.pdf']
            
            run_cmd(mycmd,outputdir)
            
            mycmd=['mv',
                   outDir+'/'+Name+'vs'+ENCName+'-'+ENCCellLines+selFile.iloc[i,0]+'.Comparison.txt',
                   outputdir+'/'+'MassSpectrometry/'+'VennDiagrams'+'/'+Name+'vs'+ENCName+'-'+ENCCellLines+selFile.iloc[i,0]+'.Comparison.txt']
            
            run_cmd(mycmd,outputdir)
    else:
        print('Suitable data is not available for the proof of piggyback binding event.'+
              ' Search in the SRA/GEO database if the narrowPeaks bed files are available '+
              'for the piggy back transcription factor in the cell lines of your interest.')
        
###############################################

def piggBack (outputdir,mPwm,sProt,Species,sFasta,lProt,peak,Name,size,threads,gVer):
    cmd_rs=['findMotifsGenome.pl','mergePeaks']
    #---  
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables of HOMER to the default PATH', 
                          'green', attrs=['bold']))
            exit()
            
    #---  
    out_motifs=outputdir+'/Peaks/Peaks_DNAmotifs/'
    tempDir=outputdir+"/"+"temporary/"
    dirs=[outputdir+'/'+'MassSpectrometry/',tempDir,out_motifs]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    if mPwm=='NA':
        print(colored('--mPwm is missing in MassSpectromety mode. Exiting!', 
                      'green', 
                      attrs=['bold']) )
        exit()

    if sProt=='NA':
        print(colored('--sProt is missing in MassSpectromety mode. Exiting!', 
                      'green', 
                      attrs=['bold']) )
        exit()

    if Species=='NA':
        print(colored('--Species is missing in MassSpectromety mode. Exiting!', 
                      'green', 
                      attrs=['bold']) )
        exit()

    if lProt == 'NA':
        print(colored('--lProt is missing in MassSpectromety mode. Exiting!', 
                      'green', 
                      attrs=['bold']) )
        exit()
    else:
        entrezOfProteins=massHugo2Entrez(lProt,Species)

    #--- 
    motifOfgreenPipelist=[]
    motifOfgreenPipe_withoutpwm_withoutsesquencelist=[]
    MassSpectro_Fasta=pyfaidx.Fasta(sFasta)

    #-- Reading peak files, finding significant Jasper motifs, find user provided motifs and unexplained peaks with motifs
    #____________________
    motifFile=out_motifs + Name + '/' + 'knownResults.txt'

    peakFile=peak
        
    #---- Motifs in the total peaks
    if not os.path.exists(motifFile):

        motifs_peaks(gVer,peakFile,out_motifs + Name,size,threads,outputdir)
        motifOfgreenPipelist.append(motifFile)

        if not os.path.exists(motifFile):

            print(colored('Unable to calculate/find significant motifs in peak file: '+ 
                          peakFile +
                          '. Did you called peaks in your samples? You can give peak file (homer format) manually using --peak', 
                          'green', 
                          attrs=['bold']) )
            exit()

    outfile1=outputdir+'/'+'MassSpectrometry/' + Name + '.motifs_pwm.txt'
    outfile2=outputdir+'/'+'MassSpectrometry/' + Name + '.peaks_WithoutPwmSeq.txt'
    outfile3=outputdir+'/'+'MassSpectrometry/' + Name + '.motifs_WithoutPwmSeq'

    filterMotifSequence(peakFile,
                        outfile1,
                        outfile2,
                        outfile2[:-4],
                        mPwm,
                        sProt,
                        tempDir,
                        gVer,
                        size,
                        threads,
                        mDist,
                        MassSpectro_Fasta,
                        outputdir)

    motifOfgreenPipe_withoutpwm_withoutsesquencelist.append(
        outputdir+'/'+'MassSpectrometry/' + Name + '.motifs_WithoutmPwmWithoutSequence/knownResults.txt'
    )

    #-- (Intersecting the proteins of the mass spectrometry and greenCUT&RUN)

    metadata_file   = psource.resource_filename(__name__, "data/Encode_tf/metadata.tsv")
    Jasper2IDs_file = psource.resource_filename(__name__, "data/Jasper_geneIDconversion.txt")
    mVenn           = outputdir+'/'+'MassSpectrometry/'+'VennDiagrams'
    
    if not os.path.exists(Jasper2IDs_file):
        print(colored("Jasper2IDs_file files does not exist. Installation of the current tools was not proper.",
                      'green', 
                      attrs=['bold']) )
    else:
        Jasper2IDs=pd.read_csv(Jasper2IDs_file,sep='\t',header=None)  

    if not os.path.exists(mVenn):
        os.makedirs(mVenn)

    #--- (total peaks)
    motifInfo=intersect_mass_green(
        Jasper2IDs,
        motifOfgreenPipelist[i],
        outputdir+'/'+'MassSpectrometry/'+ Name + '.totalPeaks.txt',
        mPvalue)      
    #--- (peaks without PWM and Sequence of motifs)
    motifInfo=intersect_mass_green(
        Jasper2IDs,
        motifOfgreenPipe_withoutpwm_withoutsesquencelist[i],
        outputdir+'/'+'MassSpectrometry/'+ Name + '.WithoutmPwmWithoutSequence.txt',
        mPvalue)

    for j in range(0,motifInfo.shape[0]):
        inFile=outputdir+'/'+'MassSpectrometry/' + Name + '.peaks_WithoutPwmSeq.txt'
        EntrezIDs=int(motifInfo.iloc[j,9])
        HugoSymbol=motifInfo.iloc[j,0].split(':')[-1]
        outDir=outputdir+'/'+'MassSpectrometry'
        outfile1=outputdir+'/'+'MassSpectrometry/' + Name + '.peaks_WithoutPwmSeqOnly'+HugoSymbol+'.txt'
        pwm=outputdir+'/'+'MassSpectrometry/'+'temporary.pwm'
        proof_of_piggyBack (inFile,metadata_file,EntrezIDs,HugoSymbol,outDir,outfile1,distPiggy,Name,gVer,tempDir,size,threads,pwm,outputdir)

