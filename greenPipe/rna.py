
import os
import subprocess
from datetime import datetime
from termcolor import colored
import glob
import pkg_resources as psource
from greenPipe import universal
from multiprocessing import Pool
from itertools import repeat

def rnaIntegrate (outputdir, peakFiles, rnaFiles, Names):
	# peakFiles, rnaFiles and Names should be in csv format
	#_____________________________________________________________

    cmd_rs=[] #------------------------- remaining
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
    
    dirs=[outputdir+'/'+'rnaIntegrate']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)


	rnaFiles=rnaFiles.split(',')
    Names=Names.split(',')

    if peakFiles!='NA':
    	peakFiles=peakFiles.split(',')


