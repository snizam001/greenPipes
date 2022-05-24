#if mode=='qc' or mode=='all':
#-- quality control filtering
#____________________________________
import os
import subprocess
from termcolor import colored

def qc (libraryType,Name,inputdir,outputdir,thread):
    othercommands=["trim_galore","fastqc"]
    for othercommand in othercommands:
        try:
            subprocess.call([othercommand,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(othercommand+' is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH', 
                          'green', attrs=['bold']))
            exit()
    
    folders=['FastqQC','Trim_galore']
    for folder in folders:
        if not os.path.exists(outputdir+'/'+folder):
            os.makedirs(outputdir+'/'+folder)
    #--- Analysis
    total_cmd=[]
    if libraryType=="pair":
        myfile=outputdir + '/Fastq/' + Name

        if not (os.path.exists(myfile + '_control_R1.fastq.gz') or os.path.exists(myfile + '_control_R2.fastq.gz') or 
                os.path.exists(myfile + '_expr_R1.fastq.gz') or os.path.exists(myfile + '_expr_R2.fastq.gz')):
            print(colored('Input Fastq files does not exist. Check these paths:', 
                          'green', attrs=['bold']))
            print(myfile + '_control_R1.fastq.gz')
            print(myfile + '_control_R2.fastq.gz')
            print(myfile + '_expr_R1.fastq.gz')
            print(myfile + '_expr_R2.fastq.gz')
            exit()

        total_cmd.append(['trim_galore', '--paired',
                          myfile + '_control_R1.fastq.gz', myfile + '_control_R2.fastq.gz',
                          '--output_dir', outputdir+'/'+'Trim_galore',
                          '--j', '1'])
        total_cmd.append(['trim_galore', '--paired',
                          myfile + '_expr_R1.fastq.gz', myfile + '_expr_R2.fastq.gz',
                          '--output_dir', outputdir+'/'+'Trim_galore',
                          '--j', '1'])
    elif libraryType=="single":
        myfile=outputdir + '/Fastq/' + Name

        if not (os.path.exists(myfile + '_control.fastq.gz') or 
                os.path.exists(myfile + '_expr.fastq.gz')):
            print(colored('Input Fastq files does not exist. Check these paths:', 
                          'green', attrs=['bold']))
            print(myfile + '_control.fastq.gz')
            print(myfile + '_expr.fastq.gz')
            exit()

        total_cmd.append(['trim_galore',
                          myfile + '_control.fastq.gz',
                          '--output_dir', outputdir+'/'+'Trim_galore',
                          '--j', '1'])
        total_cmd.append(['trim_galore',
                          myfile + '_expr.fastq.gz',
                          '--output_dir', outputdir+'/'+'Trim_galore',
                          '--j', '1'])

    #-- checking quality of the fastq files before and after trimming
    #____________________________________
    total_fastqc=[]
    
    myfile=outputdir + '/Fastq/' + Name
    if libraryType=="pair":
        xx=['_control_R1.fastq.gz','_control_R2.fastq.gz',
            '_expr_R1.fastq.gz','_expr_R2.fastq.gz']
    else:
        xx=['_control.fastq.gz','_expr.fastq.gz']
    for x in xx:
        total_fastqc.append(['fastqc',
                             '-o',outputdir+'/'+'FastqQC',
                             myfile + x])

    myfile=outputdir + '/Trim_galore/' + Name
    if libraryType=="pair":
        xx=['_control_R1_val_1.fq.gz','_control_R2_val_2.fq.gz',
            '_expr_R1_val_1.fq.gz','_expr_R2_val_2.fq.gz']
    else:
        xx=['_control_trimmed.fq.gz','_expr_trimmed.fq.gz']
    for x in xx:
        total_fastqc.append(['fastqc',
                             '-o',outputdir+'/'+'FastqQC',
                             myfile + x])
    
    return(total_cmd,total_fastqc)


