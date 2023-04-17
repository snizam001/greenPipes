# differential peaks, overlapping peaks
import os
import subprocess
from datetime import datetime
from termcolor import colored
from multiprocessing import Pool
from itertools import repeat
from greenPipe import universal
from greenPipe import filterMotifSeq
from greenPipe import initPeakCalling
from greenPipe import ann
import pandas as pd
import upsetplot
import glob

def overlapPeaks (overFiles, overDir, overDist, outputdir, overReplace):

	cmd_rs=['mergePeaks']
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

	dirs=[overDir]
	for d in dirs:
		if not os.path.exists(d):
			os.makedirs(d)

	overFile_present=0
	for overFile in overFiles:
		if not os.path.exists(overFile):
			print(colored(overFile + " does not exist. Did forgot to call the peaks?",
				"green",
				attrs=['bold']
				)
			)
		overFile_present += 1

	if overFile_present > 0:
		exit()

	c = ['mergePeaks', '-d', str(overDist)] + overFiles + ['-venn', overDir+'/Comparison_venn',
	'-matrix', overDir+'/Comparison_matrix',
	'-prefix', overDir+'/'+'Comparison']

	universal.run_cmd(c, outputdir)

	#---- cleaning of the files:1
	extensions=['Comparison_venn','Comparison.txt']
	for extension in extensions:
		venn = pd.read_csv(overDir+'/'+extension,sep='\t',header=None)
		venn = venn.astype(str).replace(overReplace+['nan'],"",regex=True)
		venn.to_csv(overDir+'/'+extension,index = None, header=None,sep='\t')
	#----
	venn = pd.read_csv(overDir+'/'+'Comparison_venn',sep='\t',header=None)
	if (venn.shape[1]-2) < 5:

		mycmd=['homer_venn.plot.R',
		'Comparison_venn',
		overDir+'/'+'Comparison_venn']

		universal.run_cmd(mycmd, outputdir)
	else:
		print(colored('Unable to plot Venn diagram as categories are >5. '+
			'Instead of this try UpSet plot.',
			'green',
			attrs = ['bold']
			)
		)


def realDiffPeaksHomer (rdName1, rdName2, rdPeak, rdPvalue, rdFoldChange, rdSize, rdOther, outputdir, libraryType, cGVersion, sFasta, threads):
#(rdTag1, rdTag2, rdPeak, rdPvalue, rdFoldChange, rdSize, rdOther, outputdir, libraryType):
	totalCmd = []
	totalCmdFile = []

	#-----
	diffR=psource.resource_filename(__name__, "rscripts/Real_differentialPeaks.R")
	cmd_rs=[diffR]
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
	#--- checking if all neccessary tools exists or not??
	#___________________________________________________________________________

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

	#--- checking if all files and directory exists or not
	#___________________________________________________________________________
	rdTag1  = outputdir + '/Tagdirectories/' + rdName1 + '_expr'
	rdTag2  = outputdir + '/Tagdirectories/' + rdName2 + '_expr'

	if rdPeak != "None":
		rdFiles = [rdTag1, rdTag2] + rdPeak.split(":")
		rdPeak1 = rdPeak.split(":")[0]
		rdPeak2 = rdPeak.split(":")[1]

	else:
		if os.path.exists(outputdir + \
								'/Peaks/' + \
								rdName1 + \
								'_all-homer.Clean.bed'):

			rdPeak1 = outputdir + \
								'/Peaks/' + \
								 rdName1 + \
								 '_all-homer.Clean.bed'

		elif os.path.exists(outputdir +\
		 						'/Peaks/' +\
								 rdName1 +\
								  '_narrow-homer.Clean.bed'):

			rdPeak1 = outputdir +\
			 					'/Peaks/' +\
								 rdName1 +\
								  '_narrow-homer.Clean.bed'

		elif os.path.exists(outputdir +\
		 						'/Peaks/' +\
								rdName1 +\
								'_broad-homer.Clean.bed'):

			rdPeak1 = outputdir +\
			 					'/Peaks/' +\
								rdName1 +\
								'_broad-homer.Clean.bed'

		if os.path.exists(outputdir +\
		 						'/Peaks/' +\
								rdName2 +\
								'_all-homer.Clean.bed'):

			rdPeak2 = outputdir +\
			 					'/Peaks/' +\
								rdName2 +\
								'_all-homer.Clean.bed'

		elif os.path.exists(outputdir +\
		 						'/Peaks/' +\
								rdName2 +\
								'_narrow-homer.Clean.bed'):

			rdPeak2 = outputdir +\
								'/Peaks/' +\
								rdName2 +\
								'_narrow-homer.Clean.bed'

		elif os.path.exists(outputdir +\
								'/Peaks/' +\
								rdName2 +\
								'_broad-homer.Clean.bed'):

			rdPeak2 = outputdir +\
								'/Peaks/' +\
								rdName2 +\
								'_broad-homer.Clean.bed'

		rdFiles = [rdTag1, rdTag2, rdPeak1, rdPeak2]

	rdFile_present=0

	for rdFile in rdFiles:
		if not os.path.exists(rdFile):
			print(colored(rdFile + " does not exist. Did forget to call the peaks?"+
				" or forget to use initPeakCalling mode before running this. If you want to find "+
				"differential peaks of IDRs or your own peak file then assign it by using"+
				" --rdPeak dir/peakfile.bed or --rdPeak dir/peak1.bed,dir/peak2.",
				"green",
				attrs=['bold']
				)
			)
		rdFile_present += 1

	if rdFile_present > 0:
		exit()

	dirs=[outputdir+'/'+'comparePeaks',
			outputdir+'/'+'comparePeaks/Peaks',
			outputdir+'/'+'comparePeaks/Motifs',
			outputdir+'/'+'comparePeaks/Annotation',
			outputdir+'/'+'comparePeaks/GeneOntology',
			outputdir+'/'+'comparePeaks/Overlap',
			outputdir+'/'+'comparePeaks/ImageVplot',
			outputdir+'/'+'comparePeaks/PeakShift',
			outputdir+'/'+'comparePeaks/BlackAndWhite',
			outputdir+'/'+'comparePeaks/BulkyvsNonBulky',
			]

	for d in dirs:
		if not os.path.exists(d):
			os.makedirs(d)


	#-- normalizing the tagDirectories for differential peak calling
	#___________________________________________________________________________
	sExpr1 =outputdir + '/SpikeIn/' + rdName1 + '_expr.Flagstats.txt'
	sExpr2 =outputdir + '/SpikeIn/' + rdName2 + '_expr.Flagstats.txt'
	Expr1  =outputdir + '/Bamfiles/' + rdName1 + '_expr.Flagstats.txt'
	Expr2  =outputdir + '/Bamfiles/' + rdName2 + '_expr.Flagstats.txt'

	File_present = 0
	for File in [sExpr1, sExpr2, Expr1, Expr2]:
		if not os.path.exists(File):
			print(colored(rdFile + " does not exist. Copy or generate these Flagstats files using samtools.",
				"green",
				attrs=['bold']
				)
			)
		File_present += 1

	if File_present > 0:
		exit()

	oVal  =initPeakCalling.spike_normalization2([sExpr1, sExpr2, Expr1, Expr2],libraryType)
	print(colored("sExpr1 sExpr2 Expr1 Expr2","green",attrs=["bold"]))
	print(oVal)

	Expr1Spike_perReads = sExpr1/(sExpr1 + Expr1)
	Expr2Spike_perReads = sExpr2/(sExpr2 + Expr2)

	#-- normalization method copied from initPeakCalling.py: SpikeIn_normalized_controlReadsNew = outvalues[2] * ((outvalues[0]/(outvalues[2]+outvalues[0])) / (outvalues[1]/(readN_expr+outvalues[1])))
	SpikeIn_normalized_ReadsNew = Expr1 * ( Expr1Spike_perReads / Expr2Spike_perReads)

	cmd=['makeTagDirectory',
		'/tmp/xYz786',
		'-d', rdTag1,
		'-totalReads', SpikeIn_normalized_ReadsNew]

	universal.run_cmd(c, outputdir)

	#-- to find significant peaks

	cmd=['getDifferentialPeaks',
		rdPeak1,
		'/tmp/xYz786', rdTag2,
		' -size', str(rdSize),
		'-P', str(rdPvalue),
		'-F', str(rdFoldChange)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	cmd=['getDifferentialPeaks',
		rdPeak2,
		rdTag2, '/tmp/xYz786'
		' -size', str(rdSize),
		'-P', str(rdPvalue),
		'-F', str(rdFoldChange)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	#--- To find all genes

	cmd=['getDifferentialPeaks',
		rdPeak1,
		'/tmp/xYz786', rdTag2,
		' -size', str(rdSize),
		'-P', str(1),
		'-F', str(0)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + 'allPeaks.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	cmd=['getDifferentialPeaks',
		rdPeak2,
		rdTag2, '/tmp/xYz786'
		' -size', str(rdSize),
		'-P', str(1),
		'-F', str(0)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + 'allPeaks.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	#--- v plot real differential Peaks
	cmd = ['Rscript', diffR,
			'-f', outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + 'allPeaks.txt',
			'-p', rdPvalue,
			'-c', rdFoldChange,
			'-d', outputdir+'/'+'comparePeaks/ImageVplot',
			'-x', rdName1 + '-vs-' + rdName2
			]
	universal.run_cmd(cmd,outputdir)

	cmd = ['Rscript', diffR,
			'-f', outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + 'allPeaks.txt',
			'-p', rdPvalue,
			'-c', rdFoldChange,
			'-d', outputdir+'/'+'comparePeaks/ImageVplot',
			'-x', rdName2 + '-vs-' + rdName1
			]
	universal.run_cmd(cmd,outputdir)

	#--- finding black and white regions: for this I am using read less that or equal to 5 to match with background and greater than or equal to 15 to consider it as real Peak
	#--- better whould have been to use backround read distribution but not doing it.

	bw1 = pd.read_csv(outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
			comment='#',
			sep='\t',
			header = None)

	bw1 = bw1[(bw1.iloc[:,8] <= 5) & (bw1.iloc[:,7] >= 15)]

	bw1.to_csv(
	outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
	sep='\t',
	header= None,
	index = None
	)

	bw2 = pd.read_csv(outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1  + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
			comment='#',
			sep='\t',
			header = None)

	bw2 = bw2[(bw2.iloc[:,8] <= 5) & (bw2.iloc[:,7] >= 15)]

	bw2.to_csv(
	outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName2 + '-vs-' + rdName1  + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
	sep='\t',
	header= None,
	index = None
	)

	#--- Motifs, Annotations, gene Ontology

	annPrefixes=[rdName1 + '-vs-' + rdName2,
				rdName2 + '-vs-' + rdName1,
				rdName1 + '-vs-' + rdName2 + '-BlackAndWhite',
				rdName2 + '-vs-' + rdName1 + '-BlackAndWhite'
				]

	annpeakFiles=[outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt'
				]

	ann.Ann (annpeakFiles,
				annPrefix,
				cGVersion,
				sFasta,
				outputdir+ '/' + 'comparePeaks/',
				threads)

	#--- bulky vs non bulky, peak shift and Overlap














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
