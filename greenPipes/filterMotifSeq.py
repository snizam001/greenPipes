import os
import pandas as pd
import subprocess
import pkg_resources as psource
import itertools
from greenPipes import universal
from greenPipes import greenCutFrq
from termcolor import colored

def motifs_peaks (genomeversion,mypeak,out_dir,size,threads,outputdir):

    #---- find motif in the given input bed or homer file
    cmd_rs=['findMotifsGenome.pl']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.'+
                          " This tools is the part of HOMER. Please install it.",
                          'red',
                          attrs=['bold']))
            exit()

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

    universal.run_cmd(c,outputdir)

def ambDNA (seq):
    #--- ambigous DNA sequence: motif sequence e.g. TGANTCA
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

def filterMotifSequence (infile,outfile1,outfile2,outfile3,pwm,sequence,tempDir,genomeversion,size,threads,seqDist,MassSpectro_Fasta,outputdir):

    print (colored("Caution: remove chrM before running this.",
                   "red",
                   attrs = ["bold"]
                  )
          )

    if not os.path.exists(infile):
        print(colored("filterMotifSequence: " + infile+' does not exist',
                      'red',
                      attrs=['bold']
                     )
             )

    peakFile=pd.read_csv(infile,
                         sep='\t',
                         header=None)

    if peakFile.shape[1] == 3:
        peakFile["Name"]   = peakFile.iloc[:,0] + ":" + peakFile.iloc[:,1].astype(str) + "-" + peakFile.iloc[:,2].astype(str)
        peakFile["Value"]  = [1]*peakFile.shape[0]
        peakFile["Strand"] = ['+']*peakFile.shape[0]

    peakFile.to_csv(infile,sep="\t",header=None,index=None)



    print ('Total number of the peaks in the ' +
           infile + ' is :' +
           str(peakFile.shape[0]))
    #-----
    pwms=pwm.split(',')

    for pw in pwms:
        print (pw)
        greenCutFrq.extrMotif([pw])
        mycmd=['findMotifsGenome.pl',
               infile,
               genomeversion,
               tempDir,
               '-size', str(size),
               '-nomotif',
               '-p', str(threads),
               '-find', 'temp.motif']

        with open(outfile1, "w+") as f:
            universal.run_cmd_file(mycmd,f,outputdir)

        pwmPeaks=pd.read_csv(outfile1,sep='\t')

        peakFile=peakFile[~peakFile.iloc[:,3].isin(pwmPeaks['PositionID'])]
        print ('Total number of the peaks after filtering peaks'+
               ' containing pwmProteinOfInterest motifs : ' + pw + " " +
               str(peakFile.shape[0]))

    #-----
    for s in sequence.split(','):
        seqs=ambDNA(s)

        print("Possible sequences for given motif is:")
        print(seqs)
        peakWithSeq=[]
        for i in range(0,peakFile.shape[0]):
            ID=peakFile.iloc[i,0]
            if 'chrM' in ID:
                print (colored("Remove chrM from peak file before running this.",
                               "red",
                               attrs = ["bold"]
                              )
                      )
                exit()
        for i in range(0,peakFile.shape[0]):
            ID=peakFile.iloc[i,0]

#            Start=int((peakFile.iloc[i,1]+peakFile.iloc[i,2])/2)-seqDist
#            End=int((peakFile.iloc[i,1]+peakFile.iloc[i,2])/2)+seqDist

            Start=int(peakFile.iloc[i,1] - seqDist)
            End=int(peakFile.iloc[i,2] + seqDist)

            Name=peakFile.iloc[i,3]
            SeqPos=MassSpectro_Fasta[ID][Start:End].seq
            SeqNeg=SeqPos[::-1].translate(SeqPos[::-1].maketrans("ATGCatgc", "TACGtacg"))
            for seq in seqs:
                if seq in SeqPos or seq in SeqNeg:
                    peakWithSeq.append(Name)

        peakFile=peakFile[~peakFile.iloc[:,3].isin(peakWithSeq)]
        print ('Total number of the peaks after filtering peaks containing'+
               ' pwmProteinOfInterest motifs + sequenceProteinOfInterest '+ s + ' :' +
               str(peakFile.shape[0]))

    peakFile.to_csv(outfile2, header=None, index=None, sep='\t')

    motifs_peaks(genomeversion, outfile2, outfile3, size, threads, outputdir)
