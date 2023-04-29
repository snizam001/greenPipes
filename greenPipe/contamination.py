import os
import subprocess
from termcolor import colored
from greenPipe import universal

def contamination (libraryType,Name,inputdir,outputdir,thread):
    folders=['Contamination']
    for folder in folders:
        if not os.path.exists(outputdir+'/'+folder):
            os.makedirs(outputdir+'/'+folder)

    cmd_rs=['fastq_screen']
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

    cmd_rs=['seqtk']
    for cmd_r in cmd_rs:
    	try:
    		subprocess.call([cmd_r,'seq','--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    	except FileNotFoundError:
    		print(colored(cmd_r+
    						': It is not installed in your computer or not in the PATH.',
    						'green',
    						attrs=['bold']
    						)
    				)
    		exit()

    total_cmd=[]
    total_cmd_seq=[]
    total_cmd_out=[]
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
                '--nohits',
                '--subset', '100000',
                '--force',
                '--threads', str(thread),
                '--outdir', outputdir + '/Contamination/' + Name,
                file
                ]

                universal.run_cmd(c, outputdir)


                #total_cmd.append(c)
                #/media/sheikh/TXPN3/Nizam/greenPipe_example/Fastq/NFYA_endo_control_R1.fastq.gz
                #NFYA_endo_control_R1.tagged_filter.fastq.gz
                c=['seqtk',
                'seq',
                '-a', outputdir + '/Contamination/'  + Name + '/' + file.split('/')[-1].replace('.fastq.gz','.tagged_filter.fastq.gz')
                ]

                universal.run_cmd_fileOpen (c,
                                            outputdir + '/Contamination/' + Name + '/' + file.split('/')[-1].replace('.fastq.gz','.tagged_filter.fa'),
                                            outputdir)

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
                '--nohits',
                '--subset', '100000',
                '--force',
                '--threads', str(thread),
                '--outdir', outputdir + '/Contamination/' + Name,
                file
                ]

                universal.run_cmd(c, outputdir)

                c=['seqtk',
                'seq',
                '-a', outputdir + '/Contamination/' + Name + '/' + file.split('/')[-1].replace('.fastq.gz','.tagged_filter.fastq.gz')
                ]

                universal.run_cmd_fileOpen (
                                            c,
                                            outputdir + '/Contamination/' + Name + '/' + file.split('/')[-1].replace('.fastq.gz','.tagged_filter.fa'),
                                            outputdir
                                            )



    return(total_cmd, total_cmd_seq, total_cmd_out)
