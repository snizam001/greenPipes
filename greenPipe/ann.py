#-- Myannotations
#_________________

import os
import subprocess
from datetime import datetime
from termcolor import colored

def ann (outputdir,peakFile,annFiles,size,name,prefix,threads,annR):
    cmd_rs=[annR]
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
    
    dirs=[outputdir+'/'+'UserAnnotation']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
        
    c_files=[peakFile]+annFiles.split(',')
    print(c_files)
    e_file=0
    for c_file in c_files:
        if not os.path.exists(c_file):
            e_file=e_file+1
            print(colored(c_file+": does not exist",
                          'green',
                          attrs=['bold']
                         )
                 )
    if e_file > 0:
        exit()
        
    outFile=outputdir+'/UserAnnotation/'+prefix+'.txt'

    cmd=[annR,
         '-m', size,
         '-i', peakFile,
         '-d', annFiles,
         '-c', str(threads),
         '-n', name,
         '-o', outFile]
    return(cmd)

