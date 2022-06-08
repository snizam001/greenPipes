#if mode=='alignment' or mode=='all':
#-- alignment
#____________________________________
import os
import subprocess
from datetime import datetime
from termcolor import colored
from greenPipe import universal
      
        
def alignment (libraryType,outputdir,Name,refgenome,spikein,threads,alignParam):
    #refgenome=
    #spikein=
    #libraryType=pair/single
    #---
    if libraryType=="pair":
        myothercommand="bowtie2"
        if not os.path.exists(refgenome+".1.bt2"):
            print(colored('Bowtie2 reference genome does not exist.'+
                          'Create index of the genome file using bowtie2-build', 
                          'green', attrs=['bold']))
        if not os.path.exists(spikein+".1.bt2"):
            print(colored('Bowtie2 reference genome does not exist for SpikeIn.'+
                          ' Create index of the SpikeIn fasta file using bowtie2-build', 
                          'green', attrs=['bold']))
            
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
            bowtie2_parameters=alignParam.split(',') + ['-p', str(threads)]

    elif libraryType=="single":
        myothercommand="bwa"
        if not os.path.exists(refgenome+".bwt"):
            print(colored('BWA reference genome does not exist. '+
                          'Create index of the genome file using bowtie2-build', 
                          'green', attrs=['bold']))
        if not os.path.exists(spikein+".bwt"):
            print(colored('BWA reference genome does not exist for SpikeIn.'+
                          ' Create index of the SpikeIn fasta file using bowtie2-build', 
                          'green', attrs=['bold'])) 

        if alignParam == "None":
            bwa_parameters=['-t', str(threads)]
        else:
            bwa_parameters=['-t', str(threads)] + alignParam.split(',')

    #---
    try:
        subprocess.call([myothercommand,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(colored('bowtie2 or bwa is not installed in your computer or not in the PATH.'+
                      ' Install or copy the executables to the default PATH', 
                      'green', attrs=['bold']))
        exit()

    try:
        subprocess.call(['samtools','--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(colored('samtools is not installed in your computer or not in the PATH.'+
                      ' Install or copy the executables to the default PATH', 
                      'green', attrs=['bold']))
        exit()
    #--- creating folders
    if not os.path.exists(outputdir+'/Trim_galore/'):
        print(colored('Trim_galore folder does not exist in your output folder.'+
                      'Before running alignment, pass fastq files through quality control using qc mode of this pipeline', 
                      'green', attrs=['bold']))
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
                          'green',
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
            #----
            c=['samtools', 'index', '-@', str(threads),outfile+Expr[k]+'.bam']
            universal.run_cmd(c,outputdir)
            #----
            c=['samtools', 'flagstat', '-@', str(threads),outfile+Expr[k]+'.bam']
            with open(outfile+Expr[k]+'.Flagstats.txt', "w+") as f:
                universal.run_cmd_file(c,f,outputdir)

def rnaAlign (libraryType,outputdir,Name,refgenome,spikein,threads,alignParam):
  
  cmd_rs=['STAR']
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
  
  dirs=[outputdir+'/'+'Bamfiles']
  for d in dirs:
      if not os.path.exists(d):
          os.makedirs(d)

  #-- check if input files exists:
  if libraryType == 'pair':
    c_files = [outputdir + '/inputdir/' + Name + '_expr_R1.fastq.gz', 
              outputdir + '/inputdir/' + Name + '_expr_R2.fastq.gz'
              ]
  elif libraryType == 'single':
    c_files = [outputdir + '/inputdir/' + Name + '_expr.fastq.gz']

  e_file=0
  for c_file in c_files:
      if not os.path.exists(c_file) or os.path.getsize(c_file)==0:
          e_file=e_file+1
          print(colored(c_file+": does not exist or empty",
                        'green',
                        attrs=['bold']
                       )
               )
  if e_file > 0:
      exit()

  #-- check if references exists:
  refgenomeFiles = ["chrStart.txt",
                    "chrName.txt",
                    "chrNameLength.txt",
                    "chrLength.txt",
                    "geneInfo.tab",
                    "exonGeTrInfo.tab",
                    "transcriptInfo.tab",
                    "exonInfo.tab",
                    "sjdbList.fromGTF.out.tab",
                    "sjdbList.out.tab",
                    "sjdbInfo.txt",
                    "genomeParameters.txt",
                    "Genome",
                    "SA",
                    "SAindex"]

  for refgenomeFile in refgenomeFiles:
    if not os.path.exists(refgenome + '/' + refgenomeFile):
        print(colored(refgenomeFile + ' :does not exist or not in the path'+
                      'Create index of the genome file for STAR', 
                      'green', attrs=['bold']))
    if not os.path.exists(spikein++ '/' + refgenomeFile):
        print(colored(refgenomeFile + ' :does not exist or not in the path'+
                      ' Create index of the SpikeIn fasta file for STAR', 
                      'green', attrs=['bold']))

  if alignParam == "None":
    star_parameters=['--outSAMtype', 'BAM', 'SortedByCoordinate',
                    '--runThreadN', str(threads)]
    if 
  else:
      bowtie2_parameters=alignParam.split(',') + ['--runThreadN', str(threads)]

STAR  --runThreadN 30  --genomeDir /media/txpn/nvme/Databases/hg38/STAR 
--readFilesIn    AS-542317-LR-54102_R1.fastq.gz AS-542317-LR-54102_R2.fastq.gz
--readFilesCommand zcat      --outFileNamePrefix  ../Bamfiles/WT2_ 
--outSAMtype BAM   SortedByCoordinate 

















    elif libraryType=="single":
        myothercommand="bwa"
        if not os.path.exists(refgenome+".bwt"):
            print(colored('BWA reference genome does not exist. '+
                          'Create index of the genome file using bowtie2-build', 
                          'green', attrs=['bold']))
        if not os.path.exists(spikein+".bwt"):
            print(colored('BWA reference genome does not exist for SpikeIn.'+
                          ' Create index of the SpikeIn fasta file using bowtie2-build', 
                          'green', attrs=['bold'])) 

        if alignParam == "None":
            bwa_parameters=['-t', str(threads)]
        else:
            bwa_parameters=['-t', str(threads)] + alignParam.split(',')

    #---
    try:
        subprocess.call([myothercommand,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(colored('bowtie2 or bwa is not installed in your computer or not in the PATH.'+
                      ' Install or copy the executables to the default PATH', 
                      'green', attrs=['bold']))
        exit()

    try:
        subprocess.call(['samtools','--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(colored('samtools is not installed in your computer or not in the PATH.'+
                      ' Install or copy the executables to the default PATH', 
                      'green', attrs=['bold']))
        exit()
    #--- creating folders
    if not os.path.exists(outputdir+'/Trim_galore/'):
        print(colored('Trim_galore folder does not exist in your output folder.'+
                      'Before running alignment, pass fastq files through quality control using qc mode of this pipeline', 
                      'green', attrs=['bold']))
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
                          'green',
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
            #----
            c=['samtools', 'index', '-@', str(threads),outfile+Expr[k]+'.bam']
            universal.run_cmd(c,outputdir)
            #----
            c=['samtools', 'flagstat', '-@', str(threads),outfile+Expr[k]+'.bam']
            with open(outfile+Expr[k]+'.Flagstats.txt', "w+") as f:
                universal.run_cmd_file(c,f,outputdir)
