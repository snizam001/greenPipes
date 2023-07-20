#-- SelectingEqualReads mode: this will not run if you use 'all' mode
#_________________________________________
#if mode=="SelectingEqualReads":
import os
import subprocess
from datetime import datetime
from termcolor import colored
import re
from multiprocessing import Pool
from itertools import repeat
from greenPipes import universal


def equalRead (libraryType,outputdir,Names,RandomReadNumbers,threads):
    cmd_rs=['samtools','sambamba']
    #---
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored('samtools or sambamba is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH',
                          'red', attrs=['bold']))
            exit()
    #---
    Expr=['_control','_expr']
    bamDir=[outputdir + '/Bamfiles/',outputdir + '/SpikeIn/']
    randomNumber=786
    readN_exp=[]
    print(colored('Equal number of the reads will choosen: equalRead mode is active now',
                  'green', attrs=['bold']))

    for Name in Names:
        for j in range(0,len(bamDir)):
            infile= bamDir[j] + Name
            for k in range(0,len(Expr)):
                if not os.path.exists(infile+Expr[k]+'.bam'):
                    print(colored(infile+Expr[k]+'.bam'+
                                  ': does not exist. Did you align the reads before selecting random reads?',
                                  'red', attrs=['bold']))

    #--- identifying total number of reads
    for Name in Names:
        with open(outputdir + '/Bamfiles/' + Name + '_expr.Flagstats.txt','r') as f:
            xyz = 0
            for lines in f:
                if libraryType=="pair":
                    if re.search('read1',lines):
                        print(lines.split(' ')[0])
                        readN_exp.append(int(lines.split(' ')[0]))
                else:
                    if re.search('mapped',lines) and not re.search('mate',lines):
                        print(lines.split(' ')[0])
                        readN_exp.append(int(lines.split(' ')[0]))
    print(readN_exp)

    #--- Checking if given random reads are perfect or not?
    if RandomReadNumbers == 0:
        readN_exp = [ min(readN_exp)/readN for readN in readN_exp]
    else:
        readN_Print=0
        for i in range(0,len(readN_exp)):
            if readN_exp[i]<RandomReadNumbers:
                readN_Print=readN_Print+1
                print(colored('Sample: '+Name+' do not have sufficient reads. '+
                              'Choose number of reads wisely: random reads you suggested is = '+
                              str(RandomReadNumbers) +
                              ", but sample have only reads = " + str(readN_exp[i]),
                              'red', attrs=['bold']))
                exit()
        if readN_Print==0:
            readN_exp = [RandomReadNumbers/readN_exp for readN_exp in readN_exp]

    #--- Now selecting random number of the reads
    i=-1
    for Name in Names:
        i+=1
        for j in range(0,len(bamDir)):
            infile= bamDir[j] + Name
            for k in range(0,len(Expr)):
                if os.path.exists(infile+Expr[k]+'original.bam'):
                    print(colored('Tool finds '+infile+Expr[k]+'original.bam.'+
                                  ' Did you already selected random reads in this folder? '+
                                  'Be careful! It is possible that your are selecting reads '+
                                  'from already selected bamfiles in which reads are already less.'+
                                  ' If this is the case, convert *original.bam* files to *.bam* '+
                                  'files manually or by using --reverseName_equalRead True',
                                  'red', attrs=['bold']))
                    exit()
                cmd=['mv',
                     infile+Expr[k]+'.bam',
                     infile+Expr[k]+'original.bam']
                universal.run_cmd(cmd,outputdir)
                cmd=['mv',
                     infile+Expr[k]+'.bam.bai',
                     infile+Expr[k]+'original.bam.bai']
                universal.run_cmd(cmd,outputdir)
                cmd=['mv',
                     infile+Expr[k]+'.Flagstats.txt',
                     infile+Expr[k]+'original.Flagstats.txt']
                universal.run_cmd(cmd,outputdir)
                cmd=['sambamba',
                     'view',
                     '-h',
                     '-t', str(threads),
                     '-s', str(readN_exp[i]),
                     '-f', 'bam',
                     '--subsampling-seed' + '=' + str(randomNumber),
                       infile+Expr[k]+'original.bam',
                       '-o', infile+Expr[k]+'.bam']
                universal.run_cmd(cmd,outputdir)
                cmd=["samtools",
                     "index",
                     infile+Expr[k]+'.bam'
                    ]
                universal.run_cmd(cmd,outputdir)
                cmd=['samtools',
                     'flagstat',
                     infile+Expr[k]+'.bam']
                with open(infile+Expr[k]+'.Flagstats.txt', "w+") as f:
                     universal.run_cmd_file(cmd,f,outputdir)


def reverseName (Names,outputdir,threads):
    Expr=['_control','_expr']
    bamDir=[outputdir + '/Bamfiles/',outputdir + '/SpikeIn/']
    for Name in Names:
        for j in range(0,len(bamDir)):
            infile= bamDir[j] + Name
            for k in range(0,len(Expr)):
                for l in ['.bam','.bam.bai','.Flagstats.txt']:
                    c=['mv',
                       infile+Expr[k]+'original'+l,
                       infile+Expr[k]+l]
                    try:
                        universal.run_cmd(c,outputdir)
                    except:
                        pass
