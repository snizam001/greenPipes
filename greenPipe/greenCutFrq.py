#-- Cut-frequency around motifs
#____________________________________

import os
import pandas as pd
import pickle
import subprocess
from datetime import datetime
from termcolor import colored
from greenPipe import initPeakCalling
from greenPipe import qcTagD
import pkg_resources as psource
from io import StringIO
from greenPipe import universal

def extrMotif (maN):
    
    jFile=open(psource.resource_filename(__name__, "data/Jasper2018_homer.pickle"),'rb')
    jdata = pickle.load(jFile)
    
    if os.path.exists("temp.motif"):
        os.remove("temp.motif")
        
    for ma in maN:
        
        jMotif = pd.DataFrame(jdata[ma])[0].str.split().apply(pd.Series)
        
        jMotif.to_csv("temp.motif",
                      sep="\t",
                      header=None,
                      index=None,
                      mode="a")
        

    
def initgreenCutFrq (maN,gVersion,outputdir,threads):
    #--------
    cmd_rs=['scanMotifGenomeWide.pl','windowBed']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is not installed in your computer or not in the PATH.'+
                          " This tools is the part of HOMER/bedtools. Please install it.",
                          'green', 
                          attrs=['bold']))
            exit()
            
    #--------
    dirs=[outputdir+'/'+'motifDatabase/']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
            
    #--------
    jOut=dirs[0]+maN+".bed"
    extrMotif(maN)
    
    c = ["scanMotifGenomeWide.pl",
         "temp.motif",
         gVersion,
         "-bed",
         "-p", threads]
    
    f = dirs[0] + "/" + maN + ".bed"
    with open(f,'w') as of:
        universal.run_cmd_file (c,of,outputdir)
    #--------
    c = ['windowBed',
         '-a', f,
         '-b', f,
         '-w', '1000']
    print(colored(c, 'green', attrs=['bold']))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        process = subprocess.Popen(c,
                                   stdout=subprocess.PIPE, 
                                   stderr=logfile)
        stdout, stderr = process.communicate()
    #--
    b   = StringIO(stdout.decode('utf-8'))
    df  = pd.read_csv(b, sep="\t",header=None)
    
    df.columns=["c1",
                "s1",
                "e1",
                "seq1",
                "v1",
                "st1",
                "c2",
                "s2",
                "e2",
                "seq2",
                "v2",
                "st2"]
    
    df  = df[(df["s1"]!=df["s2"]) & (df["e1"]!=df["e2"])]
    df1 = df.iloc[:,[0,1,2,3,4,5]]
    df2 = df.iloc[:,[6,7,8,9,10,11]]
    df2.columns=df1.columns
    df  = pd.concat([df1,df2])
    df["Name"]=df.iloc[:,0].map(str)+":"+df.iloc[:,1].map(str)+"-"+df.iloc[:,2].map(str)
    bFile=pd.read_csv(f,sep="\t",header=None)
    bFile.columns=["c1",
                   "s1",
                   "e1",
                   "seq1",
                   "v1",
                   "st1"]
    
    bFile["Name"]=bFile.iloc[:,0].map(str)+":"+bFile.iloc[:,1].map(str)+"-"+bFile.iloc[:,2].map(str)
    bFile=bFile[~bFile["Name"].isin(df["Name"])]
    bFile.iloc[:,:-1].to_csv(jOut,
                             sep="\t",
                             header=None,
                             index=None)
    
    return(jOut)


def cutFreq (Names,threads,outputdir,cutMotif,cutCenter,maN):
    
    totalCmd=[]
    
    cmd_rs=['bedtools']
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
            
    dirs=[outputdir+'/'+'cut_frequency/']
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
            
    fragDir=outputdir+'/'+'Tagdirectories_qualities/'
    
    for Name in Names:

        exprFrag=fragDir + '/' + Name + '_expr.sorted.fragments.bed'
        ctrlFrag=fragDir + '/' + Name + '_control.sorted.fragments.bed'
        exprSpkFrag=fragDir + '/' + Name + '_expr.spike.sorted.fragments.bed'
        ctrlSpkFrag=fragDir + '/' + Name + '_control.spike.sorted.fragments.bed'
        
        #-------
        if not os.path.exists(exprFrag) or not os.path.exists(ctrlFrag):
            qcTagD.qcTagD_fragmentLength(outputdir,Name,threads,"True")
        
        for Frag in [exprFrag,ctrlFrag]:
            c=[ 'windowBed', 
               '-a', cutMotif, 
               '-b', Frag, 
               '-w', '1000']
            with open(Frag.replace(".bed","")+".in",'w') as of:
                universal.run_cmd_file (c,of,outputdir)
                
        #-------
        flagStatsFs=[outputdir+'/Bamfiles/'+Name+"_expr.Flagstats.txt", 
                    outputdir+'/Bamfiles/'+Name+"_control.Flagstats.txt",
                    outputdir+'/SpikeIn/'+Name+"_expr.Flagstats.txt",
                    outputdir+'/SpikeIn/'+Name+"_control.Flagstats.txt"
                   ]
        e_file=0
        for flagStatsF in flagStatsFs:
            if not os.path.exists(flagStatsF):
                e_file = e_file + 1
                print(colored(flagStatsF+": file is not present."+
                              "Did you align the files before calculating cut frequency?",
                              'green',
                              attrs = ['bold']
                             )
                     )
        if e_file > 0:
            exit()
            
        extrMotif(maN)
        
        oValues=initPeakCalling.spike_normalization2(flagStatsFs,"pair")
        
        spikeRatio=(oValues[2]/oValues[0])/(oValues[3]/oValues[1])
        
        gCutFrqR = psource.resource_filename(__name__, "rscripts/greenCutFrq.R")

        c=[gCutFrqR,
           "--bam1", exprFrag.replace(".bed",""),
           "--bam2", ctrlFrag.replace(".bed",""),
           "--motifFile", cutMotif,
           "--centreOfmotif",cutCenter,
           "--thread",str(threads),
           "--prefix",dirs[0]+Name,
           "--pwm","temp.motif",
           "--exprReads",str(oValues[0]),
           "--ctrlReads",str(oValues[1]),
           "--spikeRatio",str(spikeRatio)
          ]
        
        universal.run_cmd(c,outputdir)
        
