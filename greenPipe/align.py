#if mode=='alignment' or mode=='all':
#-- alignment
#____________________________________
import os
import subprocess
from datetime import datetime
from termcolor import colored
from greenPipe import universal


def alignment (libraryType,outputdir,Name,refgenome,spikein,threads,alignParam,gpu):
    #refgenome=
    #spikein=
    #libraryType=pair/single
    #---
    if libraryType=="pair":
      if gpu == 'False':
        myothercommand="bowtie2"
        if not os.path.exists(refgenome+".1.bt2"):
            print(colored('Bowtie2 reference genome does not exist.'+
                          'Create index of the genome file using bowtie2-build',
                          'red', attrs=['bold']))
        if not os.path.exists(spikein+".1.bt2"):
            print(colored('Bowtie2 reference genome does not exist for SpikeIn.'+
                          ' Create index of the SpikeIn fasta file using bowtie2-build',
                          'red', attrs=['bold']))

        if alignParam == "None":
            bowtie2_parameters=['--dovetail',
                                '--local',
                                '--very-sensitive-local',
                                '--no-unal',
                                '--no-mixed',
                                '--no-discordant',
                                '-I', '10',
                                '-X', '700',
                                '-p', str(threads)]
        else:
            bowtie2_parameters=alignParam.replace("[","").replace("]","").split(',') + ['-p', str(threads)]

      else:
        myothercommand="nvBowtie"
        for refFile in [refgenome + xw for xw in ['.ann','.amb','.pac','.bwt','.sa','.rpac','.rbwt','.rsa'] ]:
          r_counts = 0
          if not os.path.exists(refFile):
              print(colored(refFile+' :nvBowtie reference genome does not exist.'+
                            'Create index of the genome file using following command: '+
                            'nvBWT <input fasta file> <output index files prefix>',
                            'red', attrs=['bold']))
              r_counts += 1

          if r_counts > 0:
              exit()

        if alignParam == "None":
            nvBowtie_parameters=['--dovetail',
                                '--local',
                                '--very-sensitive-local',
                                '--no-unal',
                                '--no-mixed',
                                '--no-discordant',
                                '-I', '10',
                                '-X', '700',
                                '-p', str(threads)]
        else:
            nvBowtie_parameters=alignParam.replace("[","").replace("]","").split(',') + ['-p', str(threads)]

    elif libraryType=="single":
      if gpu == "False":
        myothercommand="bwa"
        if not os.path.exists(refgenome+".bwt"):
            print(colored('BWA reference genome does not exist. '+
                          'Create index of the genome file using bowtie2-build',
                          'red', attrs=['bold']))
        if not os.path.exists(spikein+".bwt"):
            print(colored('BWA reference genome does not exist for SpikeIn.'+
                          ' Create index of the SpikeIn fasta file using bowtie2-build',
                          'red', attrs=['bold']))

        if alignParam == "None":
            bwa_parameters=['-t', str(threads)]
        else:
            bwa_parameters=['-t', str(threads)] + alignParam.replace("[","").replace("]","").split(',')

      else:
        myothercommand="nvBowtie"
        for refFile in [refgenome + xw for xw in ['.ann','.amb','.pac','.bwt','.sa','.rpac','.rbwt','.rsa'] ]:
          r_counts = 0
          if not os.path.exists(refFile):
              print(colored(refFile+' :nvBowtie reference genome does not exist.'+
                            'Create index of the genome file using following command: '+
                            'nvBWT <input fasta file> <output index files prefix>',
                            'red', attrs=['bold']))
              r_counts += 1

          if r_counts > 0:
              exit()

        if alignParam == "None":
            nvBowtie_parameters=['--dovetail',
                                '--local',
                                '--very-sensitive-local',
                                '--no-unal',
                                '--no-mixed',
                                '--no-discordant',
                                '-I', '10',
                                '-X', '700',
                                '-p', str(threads)]
        else:
            nvBowtie_parameters=alignParam.replace("[","").replace("]","").split(',') + ['-p', str(threads)]
    #---
    try:
        subprocess.call([myothercommand,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(colored('bowtie2/bwa (if --gpu True then nvBowtie) is not installed in your computer '+
                      'or not in the PATH.'+
                      ' Install or copy the executables to the default PATH',
                      'red', attrs=['bold']))
        exit()

    try:
        subprocess.call(['samtools','--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(colored('samtools is not installed in your computer or not in the PATH.'+
                      ' Install or copy the executables to the default PATH',
                      'red', attrs=['bold']))
        exit()

    #--- creating folders
    if not os.path.exists(outputdir+'/Trim_galore/'):
        print(colored('Trim_galore folder does not exist in your output folder.'+
                      'Before running alignment, pass fastq files through quality control using qc mode of this pipeline',
                      'red', attrs=['bold']))
        exit()

    myfolders=['Bamfiles', 'SpikeIn']
    for myfolder in myfolders:
        if not os.path.exists(outputdir+'/'+myfolder):
            os.makedirs(outputdir+'/'+myfolder)

    #--- Analysis
    refGnome=[refgenome,spikein]
    bamDir=[outputdir + '/Bamfiles/',outputdir + '/SpikeIn/']
    Expr=['_control','_expr']

    #-- checking if trim_galore files are present or not
    l_exist = 0
    infile=outputdir + '/Trim_galore/' + Name
    for k in range(0,len(Expr)):
        if libraryType=="pair":
            l = [infile+ Expr[k]+'_R1_val_1.fq.gz', infile+ Expr[k]+'_R2_val_2.fq.gz']
        else:
            l = [infile+ Expr[k]+'_trimmed.fq.gz']
    for x in l:
        if not os.path.exists(x):
            l_exist = l_exist + 1
            print(colored('{} input file does not exist. Did you run qc?'.format(x),
                          'red',
                          attrs = ['bold']))
    if l_exist > 0:
        exit()

    #-- Alignment
    logfile=open(outputdir+'/'+'log.txt','a')
    logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))


    infile=outputdir + '/Trim_galore/' + Name
    for j in range(0,len(bamDir)):
        outfile= bamDir[j] + Name
        for k in range(0,len(Expr)):
          #---
          if gpu == "False":
            if libraryType=="pair":
                p1=subprocess.Popen(['bowtie2',
                                     '-x',
                                     refGnome[j] ] +
                                    bowtie2_parameters +
                                    ['-1', infile+ Expr[k]+'_R1_val_1.fq.gz',
                                     '-2', infile+ Expr[k]+'_R2_val_2.fq.gz'],
                                    stdout=subprocess.PIPE,
                                    stderr=logfile)
                p2=subprocess.Popen(['samtools',
                                     'sort',
                                     '-@', str(threads),
                                     '-O', 'BAM',
                                     '-o', outfile+Expr[k]+'.bam'],
                                    stdin=p1.stdout,
                                    stdout=subprocess.PIPE,
                                    stderr = logfile)
                stdout, stderr = p2.communicate()
                stdout, stderr
            elif libraryType=="single":
                p1=subprocess.Popen(['bwa',
                                     'mem',
                                     refGnome[j],
                                     infile+ Expr[k]+'_trimmed.fq.gz',
                                     '-o', outfile+Expr[k]+'.sam'] +
                                    bwa_parameters,
                                    stdout=subprocess.PIPE,
                                    stderr=logfile)
                p2=subprocess.Popen(['samtools',
                                     'sort',
                                     '-@', str(threads),
                                     '-O', 'BAM',
                                     '-o', outfile+Expr[k]+'.bam'],
                                    stdin=p1.stdout,
                                    stderr = logfile)
                stdout, stderr = p2.communicate()
                stdout, stderr
            c=['samtools', 'index', outfile+Expr[k]+'.bam']
            universal.run_cmd(c,outputdir)
            #----
            c=['samtools', 'flagstat', outfile+Expr[k]+'.bam']
            with open(outfile+Expr[k]+'.Flagstats.txt', "w+") as f:
                universal.run_cmd_file(c,f,outputdir)
          #-----
          elif gpu == "True":
            if libraryType=="pair":
                p1=subprocess.Popen(['nvBowtie',
                                     '-x',
                                     refGnome[j] ] +
                                    nvBowtie_parameters +
                                    ['-1', infile+ Expr[k]+'_R1_val_1.fq.gz',
                                     '-2', infile+ Expr[k]+'_R2_val_2.fq.gz'],
                                    stdout=subprocess.PIPE,
                                    stderr=logfile)
                p2=subprocess.Popen(['samtools',
                                     'sort',
                                     '-@', str(threads),
                                     '-O', 'BAM',
                                     '-o', outfile+Expr[k]+'.bam'],
                                    stdin=p1.stdout,
                                    stdout=subprocess.PIPE,
                                    stderr = logfile)
                stdout, stderr = p2.communicate()
                stdout, stderr
            elif libraryType=="single":
                p1=subprocess.Popen(['nvBowtie',
                                     '-x',
                                     refGnome[j] ] +
                                    nvBowtie_parameters +
                                    ['-U', infile+ Expr[k]+'_R1_val_1.fq.gz'],
                                    stdout=subprocess.PIPE,
                                    stderr=logfile)
                p2=subprocess.Popen(['samtools',
                                     'sort',
                                     '-@', str(threads),
                                     '-O', 'BAM',
                                     '-o', outfile+Expr[k]+'.bam'],
                                    stdin=p1.stdout,
                                    stderr = logfile)
                stdout, stderr = p2.communicate()
                stdout, stderr
            #----
            c=['samtools', 'index', outfile+Expr[k]+'.bam']
            universal.run_cmd(c,outputdir)
            #----
            c=['samtools', 'flagstat', outfile+Expr[k]+'.bam']
            with open(outfile+Expr[k]+'.Flagstats.txt', "w+") as f:
                universal.run_cmd_file(c,f,outputdir)
