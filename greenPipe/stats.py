# Users can calculate some other statistics. This script will grow with time.

import os
import subprocess
from datetime import datetime
from termcolor import colored
from multiprocessing import Pool
from itertools import repeat
from greenPipe import universal
import pandas as pd
from gtfparse import read_gtf

def pauseIndex (gtfFile, exprBam, exprSpike, ctrlBam, ctrlSpike, outputdir, tssDist, statsSpike, threads):
	cmd_rs=['featureCounts']
	for cmd_r in cmd_rs:
		try:
			subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		except FileNotFoundError:
			print(colored(cmd_r+
				': It is part of subread package. It is not installed in your computer or not in the PATH.'+
				' Install or copy the executables to the default PATH',
				'green',
				attrs=['bold']
				)
			)
			exit()

	dirs=[outputdir+'/'+'Statistics', outputdir + '/Statistics/' + 'PauseIndex']
	for d in dirs:
		if not os.path.exists(d):
			os.makedirs(d)

	exprBams = exprBam.split(',')
	ctrlBams = ctrlBam.split(',')

	if len(exprBams) != len(ctrlBams):
		print(colored("The number of the bamfiles for control and experiment"+
			" is not equal (stats mode: PolII pausing index).",
			"green",
			attrs = ['bold']
			)
		)

	totalbam_count = 0
	totalbam = exprBams + ctrlBams + gtfFile
	for file in totalbam:
		if not os.path.exists(file):
			print(colored(file + ": does not exist",
				"green",
				attrs = ['bold']
				)
			)
			totalbam_count += 1
	if totalbam_count>0:
		exit()

	totalbamfiles = []
	for i in range(0,len(exprBams)):
		totalbamfiles.append(exprBams[i])
		totalbamfiles.append(ctrlBams[i])

	#-- Finding regions for the PolII pause index calculation

	g = read_gtf(gtfFile)
	print (
		"Total number of the genes in the GTF file is: " +
		str(g.shape[0])
		)

	g = g[(g["gene_type"]=="protein_coding") & (g["feature"]=="gene")]
	print(
		"After filtering protein_coding genes: " +
		str(g.shape[0])
		)

	p = g [ g["strand"] == "+"]
	print (
		"Total genes on positive strand:" +
		str(p.shape[0])
		)

	p["TSS_start"], p["TSS_end"] = p["start"] - tssDist, p["start"] + tssDist

	p["genebody_start"], p["genebody_end"] = p["TSS_end"] + 1, p["end"]

	p = p[(p["genebody_end"]-p["genebody_start"]) > 1000]

	print(
		"After filtering with genebody size (>tssDist or >" +
		str(tssDist) +
		"), number of protein_coding genes on positive strand are: " +
		str(p.shape[0])
		)

	n = g [ g["strand"] == "-"]
	print (
		"Total genes on negative strand:" +
		str(n.shape[0])
		)

	n["TSS_start"], n["TSS_end"] = n["end"] - tssDist, n["end"] + tssDist

	n["genebody_end"], n["genebody_start"] = n["TSS_start"] - 1, n["start"]

	n = n[(n["genebody_end"] - n["genebody_start"]) > 1000]

	print(
		"After filtering with genebody size (>tssDist or >" +
		str(tssDist) +
		"), number of protein_coding genes on negative strand are: " +
		p.shape[0]
		)

	g = pd.concat([p,n])
	g = g.sort_values(
		["seqname",
		"start"]
		)

	tss = g[ ["gene_name",
		"seqname",
		"TSS_start",
	 	"TSS_end",
		"strand"] ]

	genebody = g[ ["gene_name",
		"seqname",
		"genebody_start",
		"genebody_end",
		"strand"] ]

	tss.to_csv(outputdir + '/Statistics/' + 'PauseIndex' + '/' + 'tss.saf',
		sep = '\t',
		index = None,
		header = None )

	genebody.to_csv(outputdir + '/Statistics/' + 'PauseIndex' + '/' + 'genebody.saf',
		sep = '\t',
		index = None,
		header = None )

	inFile =  outputdir + '/Statistics/' + 'PauseIndex' + '/' + 'tss.saf'
	outFile =  outputdir + '/Statistics/' + 'PauseIndex' + '/' + 'tss.counts'


	c=[featureCounts,
	'-a', inFile,
	'-T', str(threads),
	'-o', outFile] + totalbamfiles + [
	'-F', 'SAF']

	universal.run_cmd(c, outputdir)

	c=['sed', '1d', outFile, '-i']
	universal.run_cmd(c, outputdir)

	tss = pd.read_csv(outFile, sep ="\t", header = None, index = None )

	inFile =  outputdir + '/Statistics/' + 'PauseIndex' + '/' + 'genebody.saf'
	outFile =  outputdir + '/Statistics/' + 'PauseIndex' + '/' + 'genebody.counts'

	c=[featureCounts,
	'-a', inFile,
	'-T', str(threads),
	'-o', outFile] + totalbamfiles + [
	'-F', 'SAF']

	universal.run_cmd(c, outputdir)

	c=['sed', '1d', outFile, '-i']
	universal.run_cmd(c, outputdir)

	genebody = pd.read_csv(outFile, sep ="\t", header = None, index = None )

#--- Prediction to find if SpikeIn reads are true number for normalization ??????????????????????????/
