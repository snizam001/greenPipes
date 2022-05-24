import os
import subprocess
from termcolor import colored

def contamination (libraryType,Name,inputdir,outputdir,thread):
    folders=['Contamination']
    for folder in folders:
        if not os.path.exists(outputdir+'/'+folder):
            os.makedirs(outputdir+'/'+folder)

    total_cmd=[]
    if libraryType=="pair":
        myfile=outputdir + '/Fastq/' + Name

        if not (os.path.exists(myfile + '_control_R1.fastq.gz') or os.path.exists(myfile + '_control_R2.fastq.gz') or os.path.exists(myfile + '_expr_R1.fastq.gz') or os.path.exists(myfile + '_expr_R2.fastq.gz')):
            print(colored('Input Fastq files does not exist. Check these paths:',
                             'green', attrs=['bold']))
            print(myfile + '_control_R1.fastq.gz')
            print(myfile + '_control_R2.fastq.gz')
            print(myfile + '_expr_R1.fastq.gz')
            print(myfile + '_expr_R2.fastq.gz')
            exit()
        else:
            files = [
                    myfile + '_control_R1.fastq.gz',
                    myfile + '_control_R2.fastq.gz',
                    myfile + '_expr_R1.fastq.gz',
                    myfile + '_expr_R2.fastq.gz'
                    ]
            for file in files:
                c=['fastq_screen',
                '--force',
                '--threads', str(thread),
                '--outdir', outputdir + '/Contamination/' + Name,
                file
                ]

                total_cmd.append(c)

    elif libraryType=="single":
        myfile=outputdir + '/Fastq/' + Name

        if not (os.path.exists(myfile + '_control.fastq.gz') or 
                os.path.exists(myfile + '_expr.fastq.gz')):
            print(colored('Input Fastq files does not exist. Check these paths:', 
                          'green', attrs=['bold']))
            print(myfile + '_control.fastq.gz')
            print(myfile + '_expr.fastq.gz')
            exit()
        else:
            files = [
                    myfile + '_control.fastq.gz',
                    myfile + '_expr.fastq.gz'
                    ]
            for file in files:
                c=['fastq_screen',
                '--force',
                '--threads', str(thread),
                '--outdir', outputdir + '/Contamination/' + Name,
                file
                ]

                total_cmd.append(c)
    return(total_cmd)
