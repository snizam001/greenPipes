#----------- bamcoverage with without spike in
#}_______________________________________________________
import os
import subprocess
from datetime import datetime
from termcolor import colored
from greenPipe import initPeakCalling
from greenPipe import universal
import pandas as pd

def covTracks (Names,threads,outputdir,covSpike,blackListedRegions,effectiveGenomeSize,libraryType, covOtherOptions, covSpike_NormalizationFormula):

    cmd_rs=['bamCoverage']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.'+
                          " This tools is the part of deepTools. Please install it.",
                          'green',
                          attrs=['bold']))
            exit()

    dirs=[outputdir+'/'+'bamcoverage']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    experiments = ['_expr','_control']
    e_file=0
    for Name in Names:
        for experiment in experiments:
            b_file=outputdir+'/Bamfiles/'+Name+experiment+'.bam'
            if not os.path.exists(b_file):
                e_file=e_file+1
                print(colored(b_file+": does not exist",
                              'green',
                              attrs=['bold']
                             )
                     )
    if e_file > 0:
        exit()


    if covSpike == "False":
        for Name in Names:
            iCov=outputdir + '/Bamfiles/' + Name
            oCov=outputdir + '/bamcoverage/' + Name
            for experiment in experiments:
                if covOtherOptions == 'None':
                    c = ['bamCoverage',
                         '-b', iCov + experiment + '.bam',
                         '-o', oCov + experiment + '.bw',
                         '-bl', blackListedRegions,
                         '-p', str(threads),
                         '--effectiveGenomeSize', str(effectiveGenomeSize) ]
                    universal.run_cmd (c, outputdir)
                else:
                    covOtherOption = covOtherOptions.split(',')
                    c = ['bamCoverage',
                         '-b', iCov + experiment + '.bam',
                         '-o', oCov + experiment + '.bw',
                         '-bl', blackListedRegions,
                         '-p', str(threads),
                         '--effectiveGenomeSize', str(effectiveGenomeSize) ] + covOtherOption
                    universal.run_cmd (c, outputdir)
        print(colored("The coverage files are generated now. Use it to explore coverage in IGV"+
                      " or generate tracks using package https://github.com/deeptools/pyGenomeTracks",
                      "green",
                      attrs = ["bold"]
                     )
             )

    elif covSpike == "True":

        l_samples=[]
        l_samplesExpr=[]
        for Name in Names:
            l_samples.append(outputdir + '/SpikeIn/' + Name + '_expr.Flagstats.txt')
            l_samplesExpr.append(outputdir + '/Bamfiles/' + Name + '_expr.Flagstats.txt')
        vals = initPeakCalling.spike_normalization2(l_samples, libraryType)
        valsExpr = initPeakCalling.spike_normalization2(l_samplesExpr, libraryType)
        SpikeRatios = [x / y for x, y in zip(vals, valsExpr)]

        #--- change in version 3: 20 Jan, 2023
        #--- this SpikeVal was as such that spikeIn normalzation can not be compared among batches  # version 1
        if covSpike_NormalizationFormula == 1:
            SpikeVal = [ min(SpikeRatios) / SpikeRatio for SpikeRatio in SpikeRatios]  # version 1-2
        elif covSpike_NormalizationFormula == 2:
            SpikeVal = [ 0.05 / SpikeRatio for SpikeRatio in SpikeRatios]  # version 3

        print(pd.DataFrame({'Name':Names,'Spike':vals, 'Experiment':valsExpr, 'Ratio of Spike':SpikeRatios, "Normalization":  SpikeVal}))
        j=-1
        for Name in Names:
            iCov=outputdir + '/Bamfiles/' + Name
            oCov=outputdir + '/bamcoverage/' + Name
            j=j+1
            if covOtherOptions == 'None':
                c=['bamCoverage',
                   '-b', iCov+'_expr.bam',
                   '-o', oCov+'_expr.ScaledSpikeIn.bw',
                   '-bl', blackListedRegions,
                   '-p', str(threads),
                   '--effectiveGenomeSize', str(effectiveGenomeSize),
                   '--scaleFactor', str(SpikeVal[j]) ]
                universal.run_cmd (c, outputdir)
            else:
                covOtherOption=covOtherOptions.split(",")
                c=['bamCoverage',
                   '-b', iCov+'_expr.bam',
                   '-o', oCov+'_expr.ScaledSpikeIn.bw',
                   '-bl', blackListedRegions,
                   '-p', str(threads),
                   '--effectiveGenomeSize', str(effectiveGenomeSize),
                   '--scaleFactor', str(SpikeVal[j]) ] + covOtherOption
                universal.run_cmd (c, outputdir)

        print(colored("The coverage files are generated now. Use it to explore coverage in IGV"+
                      " or generate tracks using package https://github.com/deeptools/pyGenomeTracks",
                      "green",
                      attrs = ["bold"]
                     )
             )
