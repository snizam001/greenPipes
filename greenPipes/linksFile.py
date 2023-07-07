#--- Linking the files
import pandas as pd
import os

"""
greencutrun, cutrun and cuttag inputfile type:

pair-end:
control_R1 control_R2 experiment_R1 experiment_R2 givenName

single-end:
control NA experiment NA givenName

"""

def linksFile (libraryType,myin,inputdir,outputdir):
    if libraryType=="pair":
        xx=['_control_R1.fastq.gz','_control_R2.fastq.gz','_expr_R1.fastq.gz','_expr_R2.fastq.gz']
        for i in range(0,myin.shape[0]):
            for j in [0,1,2,3]:
                source = inputdir + '/' + myin.loc[i,j]
                dst = outputdir + '/Fastq/' + myin.loc[i,4] + xx[j]
                print("["+dst+"]")
                if not os.path.exists(dst):
                    os.symlink(source, dst)
                else:
                    os.unlink(dst)
                    os.symlink(source, dst)

    elif libraryType=="single":
        xx=['_control.fastq.gz','NA','_expr.fastq.gz','NA']
        for i in range(0,myin.shape[0]):
            for j in [0,2]:
                source = inputdir + '/' + myin.loc[i,j]
                dst = outputdir + '/Fastq/' + myin.loc[i,4] + xx[j]
                print("["+dst+"]")
                if not os.path.exists(dst):
                    os.symlink(source, dst)
                else:
                    os.unlink(dst)
                    os.symlink(source, dst)
