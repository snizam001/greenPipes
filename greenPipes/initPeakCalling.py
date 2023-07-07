#-- TagDirectories or initiate_peakcalling
#____________________________________
import os
import subprocess
from datetime import datetime
from termcolor import colored
from multiprocessing import Pool
from itertools import repeat
import re

def spike_normalization (FlagstatFiles, libraryType):
    values = []
    outvalues=[]
    for i in range(0,len(FlagstatFiles)):
        with open(FlagstatFiles[i]) as f:
            for lines in f:
                if libraryType=="pair":
                    if re.search('read1',lines):
                        outvalues.append(int(lines.split(' ')[0]))
                else:
                    if re.search('mapped',lines) and not re.search('mate',lines):
                        outvalues.append(int(lines.split(' ')[0]))
    values.append(outvalues[2]/(outvalues[1]/outvalues[0]))
    values.append((outvalues[1]/outvalues[0]))
    return (values)

def spike_normalization2 (FlagstatFiles, libraryType):
    outvalues=[]
    for i in range(0,len(FlagstatFiles)):
        with open(FlagstatFiles[i]) as f:
            for lines in f:
                if libraryType=="pair":
                    if re.search('read1',lines):
                        outvalues.append(int(lines.split(' ')[0]))
                else:
                    if re.search('mapped',lines) and not re.search('mate',lines):
                         outvalues.append(int(lines.split(' ')[0]))
    return (outvalues)

def initPeakCalling (libraryType,outputdir,Names,spikeNormPeak,threads):
    totalCmd=[]
    cmd_rs=['makeTagDirectory']
    #---
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': it is part of Homer. It is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH',
                          'red', attrs=['bold']))
            exit()

    dirs=[outputdir+'/'+'Tagdirectories', outputdir+'/'+'Peaks']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    #---

    homerNormFile=outputdir + '/Peaks/' + 'Homer_normalization.txt'
    prnt=('Experiment' + '\t' +
          'Read_number_expr' + '\t' +
          'Read_number_control' + '\t' +
          'Read_number_SpikeIn_expr' + '\t' +
          'Read_number_SpikeIn_control'  + '\t' +
          'SpikeIn-ratio' + '\t' +
          'SpikeIn_normalized_controlReadsOld' + '\t' +
          'perReadsSpike_expr' + '\t' +
          'perReadsSpike_ctrl' + '\t' +
          'SpikeIn_normalized_controlReadsNew'
         )

    with open(homerNormFile,'w') as NormFile:
        NormFile.write( prnt + '\n' )

    print(colored(prnt, 'green', attrs=['bold']))

    #---
    for Name in Names:

        spikein_ctrl=outputdir + '/SpikeIn/' + Name + '_control.Flagstats.txt'
        spikein_expr=outputdir + '/SpikeIn/' + Name + '_expr.Flagstats.txt'
        genome_expr=outputdir + '/Bamfiles/' + Name + '_control.Flagstats.txt'

        for txtFile in [spikein_ctrl,spikein_expr,genome_expr]:
            if not os.path.exists(txtFile):
                print(colored(txtFile+' :does not exit'+
                          'red', attrs=['bold']))
                exit()

        values=spike_normalization ([spikein_ctrl,
                                     spikein_expr,
                                     genome_expr], libraryType)

        outvalues=spike_normalization2 ([spikein_ctrl,
                                         spikein_expr,
                                         genome_expr], libraryType)

        with open(outputdir + '/Bamfiles/' + Name + '_expr.Flagstats.txt','r') as f:
            for lines in f:
                if libraryType=="pair":
                    if re.search('read1',lines):
                        readN_expr=int(lines.split(' ')[0])
                else:
                    if re.search('mapped',lines) and not re.search('mate',lines):
                        readN_exp=int(lines.split(' ')[0])

        SpikeIn_normalized_controlReadsNew = outvalues[2] * ((outvalues[0]/(outvalues[2]+outvalues[0])) / (outvalues[1]/(readN_expr+outvalues[1])))

        prnt=(Name + '\t' +
              str(readN_expr) + '\t' +
              str(outvalues[2]) + '\t' +
              str(outvalues[1]) + '\t' +
              str(outvalues[0]) + '\t' +
              str(values[1]) + '\t' +
              str(values[0]) + '\t' +
              str(outvalues[1]/(readN_expr+outvalues[1])) + '\t' +
              str(outvalues[0]/(outvalues[2]+outvalues[0])) + '\t' +
              str(SpikeIn_normalized_controlReadsNew)
             )

        with open (homerNormFile, 'a') as NormFile:
            NormFile.write(prnt+ '\n' )


        if outvalues[1]/(readN_expr+outvalues[1]) > outvalues[0]/(outvalues[2]+outvalues[0]):
            print(colored('NOTE: SpikeIn reads per true-reads are less in control compared to your experiment.'+
                          ' Consider peak calling without control (--Control_homer False) or inititate peakcalling'+
                          ' (initPeakCalling) with --spikeNormPeak False. '+
                          'Expr = ' + str(outvalues[1]/(readN_expr+outvalues[1])) + ' vs Ctrl=' + str(outvalues[0]/(outvalues[2]+outvalues[0])),
                          'red',
                          attrs=['bold','underline']))
#---
        outExpr=outputdir + '/Tagdirectories/' + Name + '_expr'
        outCtrl=outputdir + '/Tagdirectories/' + Name + '_control'
        inExpr=outputdir + '/Bamfiles/' + Name + '_expr.bam'
        inCtrl=outputdir + '/Bamfiles/' + Name + '_control.bam'

        totalCmd.append(['makeTagDirectory', outExpr, inExpr])

        if spikeNormPeak=='True':
            #--- using old normalization method
            # totalCmd.append(['makeTagDirectory', outCtrl, inCtrl,'-totalReads',str(values[0]) ])
            #--- using new normalization method
            totalCmd.append(['makeTagDirectory', outCtrl, inCtrl,'-totalReads',str(SpikeIn_normalized_controlReadsNew) ])
        else:
            totalCmd.append(['makeTagDirectory', outCtrl, inCtrl])

    print(colored(prnt, 'green', attrs=['bold']))

    return(totalCmd)
