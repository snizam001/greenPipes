# differential peaks, overlapping peaks
import os
import subprocess
from datetime import datetime
from termcolor import colored
from multiprocessing import Pool
from itertools import repeat
from greenPipe import universal
from greenPipe import filterMotifSeq
import pandas as pd
import upsetplot
import glob

def realDiffPeaksHomer (rdTag1, rdTag2, rdPeak, rdPvalue, rdFoldChange, rdSize, rdOther):

	totalCmd = []
	totalCmdFile = []

	cmd_rs=['getDifferentialPeaks']
	for cmd_r in cmd_rs:
		try:

			subprocess.call([cmd_r,'--help'],
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)

		except FileNotFoundError:
			print(colored(cmd_r+
				': It is part of Homer or bedtools. '+
				'It is not installed in your computer or not in the PATH.'+
				' Install or copy the executables to the default PATH',
				'green', attrs=['bold']))
			exit()

	rdFile_present=0
	for rdFile in [rdTag1, rdTag2, rdPeak]:
		if not os.path.exists(rdFile):
			print(colored(rdFile + " does not exist. Did forget to call the peaks?"+
				" or forget to use initPeakCalling mode before running this",
				"green",
				attrs=['bold']
				)
			)
		rdFile_present += 1

	if rdFile_present > 0:
		exit()

	cmd1=['getDifferentialPeaks', 
		rdPeak, 
		rdTag1, rdTag2, 
		' -size', str(rdSize), 
		'-P', str(rdPvalue), 
		'-F', str(rdFoldChange)] + rdOther

	totalCmd.append(cmd1)

	cmd2=['getDifferentialPeaks', 
		rdPeak, 
		rdTag1, rdTag2, 
		' -size', str(rdSize), 
		'-P', str(1), 
		'-F', str(0)] + rdOther

	totalCmd.append(cmd2)





	
"""
























def differential_peaks (genomeversion,out_dir,out_dir_annotation,mypeaks,tag1,tag2,firstExpr,secondExpr,pvalueD,foldchangeD):
     if not os.path.exists(out_dir):
          os.makedirs(out_dir)
     if not os.path.exists(mypeak):
          print(mypeak+' does not exist')
     if not os.path.exists(tag1):
          print(tag1+' does not exist')
     if not os.path.exists(tag2):
          print(tag2+' does not exist')
     i=[firstExpr,secondExpr]
     mycmd=['getDifferentialPeaks', mypeaks, tag1, tag2, ' -size', '200', '-P', str(pvalueD), '-F', str(foldchangeD)]
     with open(out_dir+i[0]+'.vs.'+i[1]+'.differentialPeaks.txt', "w+") as f:
          run_cmd_file(mycmd,f)              
     mycmd=['getDifferentialPeaks', mypeaks, tag1, tag2, ' -size', '200', '-F', '0', '-P', '1']
     with open(out_dir+i[0]+'.vs.'+i[1]+'.TotalPeaks.txt', "w+") as f:
          run_cmd_file(mycmd,f)



#-- Differential peaks   
for i in combinations(samples,2):
     tag1 = outputdir + '/Tagdirectories/' + i[0] + '_expr'
     tag2 = outputdir + '/Tagdirectories/' + i[1] + '_expr'
     mypeak = outputdir + '/Peaks/' + i[0] + '.Clean.bed'
     differential_peaks(genomeversion,out_dir,out_dir_annotation,mypeak,tag1,tag2,i[0],i[1],pvalueD,foldchangeD)
     motifs_peaks (genomeversion,mypeak,out_motifs+i[0]+'.vs.'+i[1],size,threads)
     inpeak = out_dir+i[0]+'.vs.'+i[1]+'.differentialPeaks.txt'
     outfile = out_dir_annotation+i[0]+'.vs.'+i[1]
     annotation_peaks (inpeak,genomeversion,outfile)
     #---
     tag1 = outputdir + '/Tagdirectories/' + i[1] + '_expr'
     tag2 = outputdir + '/Tagdirectories/' + i[0] + '_expr'
     mypeak = outputdir + '/Peaks/' + i[1] + '.Clean.bed'
     differential_peaks(genomeversion,out_dir,out_dir_annotation,mypeak,tag1,tag2,i[1],i[0],pvalueD,foldchangeD)
     motifs_peaks (genomeversion,mypeak,out_motifs+i[1]+'.vs.'+i[0],size,threads)
     inpeak = out_dir+i[1]+'.vs.'+i[0]+'.differentialPeaks.txt'
     outfile = out_dir_annotation+i[1]+'.vs.'+i[0]
     annotation_peaks (inpeak,genomeversion,outfile)






def peakAnn (mypeaks,genomeversion,outdir):

    cmd=['annotatePeaks.pl', 
        mypeaks, 
        genomeversion, 
        '-go', outdir+'_GO', 
        '-genomeOntology', 
        outdir+'_GenomeOntology']
    
    retun (cmd)

def peakMotifs (genomeversion,mypeak,out_dir,size,threads):

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    jFile=psource.resource_filename(__name__, "data/Jasper2018_homer.txt")

    cmd=['findMotifsGenome.pl',
        mypeak, genomeversion, 
        out_dir,
        '-size', str(size),
        '-nomotif',
        '-p', str(threads),
        '-mknown', jFile]

    return(cmd)


































"""