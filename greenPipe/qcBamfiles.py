#-- count number of reads
#____________________________________

import pandas as pd 
import more_itertools as mit
import pathlib
import gzip
import re
from greenPipe import universal
import os
import subprocess
from termcolor import colored

def countReads (Names, fastQs, FlagstatFiles, outputdir, libraryType):

    dirs=[outputdir + '/Bamfiles_QC/']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    File_present=0

    for fastQ in fastQs:
        if not os.path.exists(fastQ):
            print(colored(fastQ + " does not exist. Path of the file is correct? check --inputdir",
                "green",
                attrs=['bold']
                )
            )
            File_present += 1

    for FlagstatFile in FlagstatFiles:
        if not os.path.exists(FlagstatFile):
            print(colored(FlagstatFile + " does not exist. Path of the file is correct? " +
                "Manually prepare flagstat files using command: samtools flagstat input.bam",
                "green",
                attrs=['bold']
                )
            )
            File_present += 1

    if File_present > 0:
        exit(File_present)

    #---- counts in the fastq files
    #__________________________________________

    fastQs_counts = []
    for i in range(0,len(fastQs)):
        print(fastQs[i])
        if '.gz' in pathlib.Path(fastQs[i]).suffixes:
            with gzip.open(fastQs[i],'r') as fqF:
                for count, line in enumerate(fqF):
                    pass
            fastQs_counts.append(count)
        else:
            with open(fastQs[i],'r') as fqF:
                for count, line in enumerate(fqF):
                    pass
            fastQs_counts.append(count)
    fastQs_counts=list(mit.chunked(fastQs_counts, 4))


    #---- counts in the bamfiles
    #__________________________________________
    align_counts = []
    for i in range(0,len(FlagstatFiles)):
        with open(FlagstatFiles[i]) as f:
            for lines in f:
                if libraryType=="pair":
                    if re.search('read1',lines):
                        outvalues=int(lines.split(' ')[0])
                else:
                    if re.search('mapped',lines) and not re.search('mate',lines):
                         outvalues=int(lines.split(' ')[0])
            align_counts.append(outvalues)
    align_counts=list(mit.chunked(align_counts, 4))


    #---- Merging all together
    #__________________________________________

    total = pd.concat([
              pd.DataFrame(Names, columns = ['Name']),
              pd.DataFrame(fastQs_counts, columns = [
                                        'OriginalExpr(FQ)', 
                                        'OriginalCtrl(FQ)',
                                        'FilterExpr(FQ)', 
                                        'FilterCtrl(FQ)',
                                        ]), 
              pd.DataFrame(align_counts, columns = [
                                        'OrgExpr(Align)', 
                                        'SpikeExpr(Align)',
                                        'OrgCtrl(Align)', 
                                        'SpikeCtrl(Align)',
                                        ])
              ], axis = 1) 

    totalDiffExpr = (total.iloc[:,3] - (total.iloc[:,5] + total.iloc[:,6])) / total.iloc[:,3]
    totalDiffCtrl = (total.iloc[:,4] - (total.iloc[:,7] + total.iloc[:,8])) / total.iloc[:,4]

    y = -1
    notice = []
    for x in totalDiffExpr.tolist():
        y += 1
        if x >= 0.1:
            notice.append('Check contamination')
            print(colored( Names[y] + ': Experiment -> Is it contaminated? Frequency of reads not aligning: ' + str(x),
                'green',
                attrs = ['bold']
                )
            )
        else:
            notice.append('-')

    noticeExpr = notice

    y = -1
    notice = []
    for x in totalDiffCtrl.tolist():
        y += 1
        if x >= 0.1:
            notice.append('Check contamination')
            print(colored( Names[y] + ': Control -> Is it contaminated? Frequency of reads not aligning: ' + str(x),
                'green',
                attrs = ['bold']
                )
            )
        else:
            notice.append('-')

    noticeCtrl = notice

    total = pd.concat([total,
                    pd.DataFrame(noticeExpr, columns=['ForExpr']),
                    pd.DataFrame(noticeCtrl,columns=['ForCtrl'])
                    ],
        axis=1)

    total.to_csv(outputdir + '/Bamfiles_QC/ReadCounts.txt',sep="\t", index = None)

def qc_bam (qcBfiles,outputdir,threads,blackListedRegions,Names):

    cmd_rs=['multiBamSummary']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                ': It is part of deeptools. It is not installed in your computer or not in the PATH.'+
                ' Install or copy the executables to the default PATH',
                'green',
                attrs=['bold']
                )
            )
            exit()

    dirs=[outputdir + '/Bamfiles_QC/']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    bamFile_present=0
    for qcBfile in qcBfiles.split(','):
        if not os.path.exists(qcBfile):
            print(colored(qcBfile + " does not exist. Did you forget to generate bamfiles?",
                "green",
                attrs=['bold']
                )
            )
            bamFile_present = 1

    if bamFile_present > 0:
        exit()

    output = outputdir + '/Bamfiles_QC/Correlation' 
    outRawCount = outputdir + '/Bamfiles_QC/CorrelationRawCounts.txt' 
    c=[
    'multiBamSummary',
    'bins',
    '-b'] + qcBfiles.split(',') + [
    '-o', output + '.gz',
    '-l'] + Names + [
    '-bl', blackListedRegions,
    '-p', threads,
    '--outRawCounts', outRawCount
    ]

    universal.run_cmd(c, outputdir)

    c=[
    'plotCorrelation',
    '-in', output + '.gz',
    '-c', 'pearson',
    '-p', 'heatmap',
    '-o', output + '.jpeg',
    '--colorMap', 'RdYlBu',
    '--plotNumbers',
    '--removeOutliers',
    '--outFileCorMatrix', output + '-mat.tab'
    ]

    universal.run_cmd(c, outputdir)
