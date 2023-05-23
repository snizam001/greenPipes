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
import pkg_resources as psource
import random
import statistics

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


def realDiffPeaksHomer (rdName1, rdName2, rdPeak, rdPvalue, rdFoldChange, rdSize, rdOther, outputdir, libraryType, cGVersion, sFasta, genomeFile, threads):
#(rdTag1, rdTag2, rdPeak, rdPvalue, rdFoldChange, rdSize, rdOther, outputdir, libraryType):
	totalCmd = []
	totalCmdFile = []

	if genomeFile == "NA":
		genomeFile = psource.resource_filename(__name__, "data/hg38.genome")
		print(colored('Using human genome version hg38 specific genome file. If your organism is different,'+
		' see option --genomeFile.',
		'red',
		attrs=['bold']
		)
		)

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

	cmd_rs=['getDifferentialPeaks','shuffleBed','getPeakTags']
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

	if sFasta == "NA":
		print(colored("--sFasta is missing",
					"green",
					attrs = ["bold"]
					)
		)

	if rdOther == "None":
		rdOther = []
	else:
		rdOther = rdOther.replace("[","").replace("]","").split(',')

	#--- checking if all files and directory exists or not
	#___________________________________________________________________________
	rdTag1  = outputdir + '/Tagdirectories/' + rdName1 + '_expr'
	rdTag2  = outputdir + '/Tagdirectories/' + rdName2 + '_expr'

	# if rdPeak != "None":
	rdFiles = [rdTag1, rdTag2] + rdPeak.split(":")
	rdPeak1 = rdPeak.split(":")[0]
	rdPeak2 = rdPeak.split(":")[1]

	# else:
	# 	if os.path.exists(outputdir + \
	# 							'/Peaks/' + \
	# 							rdName1 + \
	# 							'_all-homer.Clean.bed'):
	#
	# 		rdPeak1 = outputdir + \
	# 							'/Peaks/' + \
	# 							 rdName1 + \
	# 							 '_all-homer.Clean.bed'
	#
	# 	elif os.path.exists(outputdir +\
	# 	 						'/Peaks/' +\
	# 							 rdName1 +\
	# 							  '_narrow-homer.Clean.bed'):
	#
	# 		rdPeak1 = outputdir +\
	# 		 					'/Peaks/' +\
	# 							 rdName1 +\
	# 							  '_narrow-homer.Clean.bed'
	#
	# 	elif os.path.exists(outputdir +\
	# 	 						'/Peaks/' +\
	# 							rdName1 +\
	# 							'_broad-homer.Clean.bed'):
	#
	# 		rdPeak1 = outputdir +\
	# 		 					'/Peaks/' +\
	# 							rdName1 +\
	# 							'_broad-homer.Clean.bed'
	#
	# 	if os.path.exists(outputdir +\
	# 	 						'/Peaks/' +\
	# 							rdName2 +\
	# 							'_all-homer.Clean.bed'):
	#
	# 		rdPeak2 = outputdir +\
	# 		 					'/Peaks/' +\
	# 							rdName2 +\
	# 							'_all-homer.Clean.bed'
	#
	# 	elif os.path.exists(outputdir +\
	# 	 						'/Peaks/' +\
	# 							rdName2 +\
	# 							'_narrow-homer.Clean.bed'):
	#
	# 		rdPeak2 = outputdir +\
	# 							'/Peaks/' +\
	# 							rdName2 +\
	# 							'_narrow-homer.Clean.bed'
	#
	# 	elif os.path.exists(outputdir +\
	# 							'/Peaks/' +\
	# 							rdName2 +\
	# 							'_broad-homer.Clean.bed'):
	#
	# 		rdPeak2 = outputdir +\
	# 							'/Peaks/' +\
	# 							rdName2 +\
	# 							'_broad-homer.Clean.bed'

		# rdFiles = [rdTag1, rdTag2, rdPeak1, rdPeak2]

	rdFile_present=0
	# print(rdFiles)
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
			rdFile_present = rdFile_present + 1

	if rdFile_present > 0:
		exit()


	dirs=[outputdir+'/'+'comparePeaks',
			outputdir+'/'+'comparePeaks/Peaks',
			outputdir+'/'+'comparePeaks/Motifs',
			outputdir+'/'+'comparePeaks/Annotation',
			outputdir+'/'+'comparePeaks/GeneOntology',
			# outputdir+'/'+'comparePeaks/Overlap',
			outputdir+'/'+'comparePeaks/ImageVplot',
			outputdir+'/'+'comparePeaks/PeakShift',
			outputdir+'/'+'comparePeaks/BlackAndWhite',
			outputdir+'/'+'comparePeaks/BulkyvsNonBulky',
			]

	for d in dirs:
		if not os.path.exists(d):
			os.makedirs(d)

	print("[#-------------- Normalizing with spikeIn]")
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
	Expr1Spike_perReads =  oVal[0]/ (oVal [0] + oVal[2]) # sExpr1/(sExpr1 + Expr1)
	Expr2Spike_perReads =  oVal[1]/ (oVal [1] + oVal[3]) # sExpr2/(sExpr2 + Expr2)

	#-- normalization method copied from initPeakCalling.py: SpikeIn_normalized_controlReadsNew = outvalues[2] * ((outvalues[0]/(outvalues[2]+outvalues[0])) / (outvalues[1]/(readN_expr+outvalues[1])))
	SpikeIn_normalized_ReadsNew =  oVal[2] * ( Expr1Spike_perReads / Expr2Spike_perReads) # Expr1 * ( Expr1Spike_perReads / Expr2Spike_perReads)
#######################################33
	cmd=['makeTagDirectory',
		'./xYz786',
		'-d', rdTag1,
		'-totalReads', str(SpikeIn_normalized_ReadsNew)]

	universal.run_cmd(cmd, outputdir)
#######################################33
	#-- to find significant peaks
	print("[#-------------- Finding differential peaks]")
	cmd=['getDifferentialPeaks',
		rdPeak1,
		'./xYz786', rdTag2,
		' -size', str(rdSize),
		'-P', str(rdPvalue),
		'-F', str(rdFoldChange)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	cmd=['getDifferentialPeaks',
		rdPeak2,
		rdTag2, './xYz786'
		' -size', str(rdSize),
		'-P', str(rdPvalue),
		'-F', str(rdFoldChange)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	#--- To find all genes

	cmd=['getDifferentialPeaks',
		rdPeak1,
		'./xYz786', rdTag2,
		' -size', str(rdSize),
		'-P', str(1),
		'-F', str(0)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + 'allPeaks.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	cmd=['getDifferentialPeaks',
		rdPeak2,
		rdTag2, './xYz786'
		' -size', str(rdSize),
		'-P', str(1),
		'-F', str(0)] + rdOther

	file = outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + 'allPeaks.txt'
	with open(file,'w') as f:
		universal.run_cmd_file(cmd,f,outputdir)

	print("[#-------------- Generating V plot]")
	#--- v plot real differential Peaks
	try:
		cmd = ['Rscript', diffR,
				'-f', outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + 'allPeaks.txt',
				'-p', rdPvalue,
				'-c', rdFoldChange,
				'-d', outputdir+'/'+'comparePeaks/ImageVplot',
				'-x', rdName1 + '-vs-' + rdName2
				]
		universal.run_cmd(cmd,outputdir)
	except:
		print(colored('There is error in generated differential peaks V plot. '+
				'It is possible that none of the differential peaks were found'+
				'. Check file: ' + outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + 'allPeaks.txt',
				'red', attrs=['bold']))
	try:
		cmd = ['Rscript', diffR,
				'-f', outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + 'allPeaks.txt',
				'-p', rdPvalue,
				'-c', rdFoldChange,
				'-d', outputdir+'/'+'comparePeaks/ImageVplot',
				'-x', rdName2 + '-vs-' + rdName1
				]
		universal.run_cmd(cmd,outputdir)
	except:
		print(colored('There is error in generated differential peaks V plot. '+
				'It is possible that none of the differential peaks were found'+
				'. Check file: ' + outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName2 + '-' + 'allPeaks.txt',
				'red', attrs=['bold']))

	print("[#-------------- Finding black and white peaks: using simulation]")
	#--- finding black and white regions: for this I am using read less that or equal to 5 to match with background and greater than or equal to 15 to consider it as real Peak
	#--- better would have been to use backround read distribution and now adding it.

	if os.path.getsize(rdPeak1) != 0 and os.path.getsize(outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt') != 0:
		pk1 = pd.read_csv(rdPeak1,
				comment='#',
				sep='\t',
				header = None)
		pk1_length = pk1.shape[0]
		pk1_SimulationTimes = int(round(40000/pk1.shape[0],0))
		print(pk1_SimulationTimes)
		print(colored('Performing simulation for '+str(pk1_SimulationTimes)+ ' times to find on an average background reads in finding black and white peaks.',
				'green', attrs=['bold']))

		bw1 = pd.read_csv(outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				comment='#',
				sep='\t',
				header = None)

		# ./xYz786
		# rdPeak1
		open(outputdir+ '/' + 'comparePeaks/BlackAndWhite/random-'+rdName1+'.bed', "a").close()
		for pk1_SimulationTime in range(0,pk1_SimulationTimes):
			mySeed = random.randint(100,10000)
			cmd = ['shuffleBed',
					'-i', rdPeak1,
					'-g', genomeFile,
					'-seed', str(mySeed)]

			with open(outputdir+ '/' + 'comparePeaks/BlackAndWhite/random-'+rdName1+'.bed', "a") as f:
					universal.run_cmd_file(cmd, f, outputdir)

		cmd = ['getPeakTags',
				outputdir+ '/' + 'comparePeaks/BlackAndWhite/random-'+rdName1+'.bed',
				'./xYz786']

		with open(outputdir+ '/' + 'comparePeaks/BlackAndWhite/randomReadCount-'+rdName1+'.txt', "a") as f:
			universal.run_cmd_file(cmd, f, outputdir)

		randomReadCount = pd.read_csv(outputdir+ '/' + 'comparePeaks/BlackAndWhite/randomReadCount-'+rdName1+'.txt', sep ='\t', header = None)

		xx = randomReadCount.iloc[:,1].to_numpy()
		md1 = statistics.median(xx)
		#--- The default constant = 1.4826 (approximately = 1/qnorm(3/4)) ensures consistency and use of 3 to cover 99% area of distrbution
		mad1 = 3*(1.4826*statistics.median([abs(number-md1) for number in xx]))

		mCutoff = md1 + mad1

		print(colored('To find black and white for :' + rdName1 + ', cutoff value is: backgroup <= '+str(mCutoff)+' and peak is: 4*' + str(mCutoff),
				'green', attrs=['bold']))

		bw1 = bw1[(bw1.iloc[:,8] <= mCutoff) & (bw1.iloc[:,7] >= (4*mCutoff))]

		bw1.to_csv(
			outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
			sep='\t',
			header= None,
			index = None
			)
	else:
		print(colored('There is error in Finding black and white peaks for '+ rdName1 +
				'It is possible that none of the differential peaks were found. Check file: '+
				outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				'red', attrs=['bold']))

	if os.path.getsize(rdPeak2) != 0 and os.path.getsize(outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt') != 0:
		pk2 = pd.read_csv(rdPeak2,
				comment='#',
				sep='\t',
				header = None)
		pk2_length = pk2.shape[0]
		pk2_SimulationTimes = int(round(40000/pk2.shape[0],0))
		print(pk2_SimulationTimes)
		print(colored('Performing simulation for '+str(pk2_SimulationTimes)+ ' times to find on an average background reads in finding black and white peaks.',
				'green', attrs=['bold']))

		bw2 = pd.read_csv(outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1  + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				comment='#',
				sep='\t',
				header = None)

		open(outputdir+ '/' + 'comparePeaks/BlackAndWhite/random-'+rdName2+'.bed', "a").close()
		for pk2_SimulationTime in range(0,pk2_SimulationTimes):
			mySeed = random.randint(100,10000)
			cmd = ['shuffleBed',
					'-i', rdPeak2,
					'-g', genomeFile,
					'-seed', str(mySeed)]

			with open(outputdir+ '/' + 'comparePeaks/BlackAndWhite/random-'+rdName2+'.bed', "a") as f:
					universal.run_cmd_file(cmd, f, outputdir)

		cmd = ['getPeakTags',
				outputdir+ '/' + 'comparePeaks/BlackAndWhite/random-'+rdName2+'.bed',
				rdTag2]

		with open(outputdir+ '/' + 'comparePeaks/BlackAndWhite/randomReadCount-'+rdName2+'.txt', "a") as f:
			universal.run_cmd_file(cmd, f, outputdir)

		randomReadCount = pd.read_csv(outputdir+ '/' + 'comparePeaks/BlackAndWhite/randomReadCount-'+rdName2+'.txt', sep ='\t', header = None)

		xx = randomReadCount.iloc[:,1].to_numpy()
		md1 = statistics.median(xx)
		#--- The default constant = 1.4826 (approximately = 1/qnorm(3/4)) ensures consistency and use of 3 to cover 99% area of distrbution
		mad1 = 3*(1.4826*statistics.median([abs(number-md1) for number in xx]))

		mCutoff = md1 + mad1

		print(colored('To find black and white for :' + rdName2 + ', cutoff value is: backgroup <= '+str(mCutoff)+' and peak is: 4*' + str(mCutoff),
				'green', attrs=['bold']))

		bw2 = bw2[(bw2.iloc[:,8] <= mCutoff) & (bw2.iloc[:,7] >= (4*mCutoff))]

		bw2.to_csv(
		outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName2 + '-vs-' + rdName1  + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
		sep='\t',
		header= None,
		index = None
		)
	else:
		print(colored('There is error in Finding black and white peaks for '+ rdName1 +
				'It is possible that none of the differential peaks were found. Check file: '+
				outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				'red', attrs=['bold']))

	print("[#-------------- Motifs, Annotations, gene Ontology]")
	#--- Motifs, Annotations, gene Ontology

	annPrefixes=[rdName1 + '-vs-' + rdName2,
				rdName2 + '-vs-' + rdName1,
				rdName1 + '-vs-' + rdName2 + '-BlackAndWhite',
				rdName2 + '-vs-' + rdName1 + '-BlackAndWhite'
				]
	annPrefixes=",".join(annPrefixes)
	annpeakFiles=[outputdir+ '/' + 'comparePeaks/Peaks/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				outputdir+ '/' + 'comparePeaks/Peaks/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName1 + '-vs-' + rdName2 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt',
				outputdir+ '/' + 'comparePeaks/BlackAndWhite/' + rdName2 + '-vs-' + rdName1 + '-' + str(rdPvalue) + '-'+ str(rdFoldChange) + '.txt'
				]
	annpeakFiles=",".join(annpeakFiles)
	ann.Ann (annpeakFiles,
				annPrefixes,
				cGVersion,
				sFasta,
				outputdir+ '/' + 'comparePeaks/',
				threads)

	#--- bulky vs non bulky, peak shift and Overlap
