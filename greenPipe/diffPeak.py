# differential peaks, overlapping peaks
import os
import subprocess
from datetime import datetime
from termcolor import colored
from multiprocessing import Pool
from itertools import repeat
from greenPipe import universal
from greenPipe import filterMotifSeq

filterMotifSeq.motifs_peaks (genomeversion,mypeak,out_dir,size,threads,outputdir)

def annotation_peaks (mypeaks,genomeversion,outdir):
	cmd=['annotatePeaks.pl', 
	mypeaks, 
	genomeversion, 
	'-go', outdir+'_GO', 
	'-genomeOntology', 
	outdir+'_GenomeOntology']

	retun (cmd)


def differential_peaks (genomeversion,out_dir,out_dir_annotation,mypeaks,tag1,tag2,firstExpr,secondExpr,pvalueD,foldchangeD):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	if not os.path.exists(mypeak):
		print(mypeak+' does not exist')
	if not os.path.exists(tag1):
		print(tag1+' does not exist')
	if not os.path.exists(tag2):
		print(tag2+' does not exist')
	i=[firstExpr,secondExpr]
	mycmd=['getDifferentialPeaks', mypeaks, tag1, tag2, ' -size', '200', '-P', str(pvalueD), '-F', str(foldchangeD)]
	with open(out_dir+i[0]+'.vs.'+i[1]+'.differentialPeaks.txt', "w+") as f:
		run_cmd_file(mycmd,f)			
	mycmd=['getDifferentialPeaks', mypeaks, tag1, tag2, ' -size', '200', '-F', '0', '-P', '1']
	with open(out_dir+i[0]+'.vs.'+i[1]+'.TotalPeaks.txt', "w+") as f:
		run_cmd_file(mycmd,f)

def motifs_peaks (genomeversion,mypeak,out_dir,size,threads):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	mycmd=['findMotifsGenome.pl', mypeak, genomeversion, out_dir, '-size', str(size), '-nomotif', '-p', str(threads), '-mknown', '/media/txpn/nvme/Databases/Jasper2018_homer.txt'] #------------ change here the location of the file of the JASPER motifs 2018
	run_cmd(mycmd)


def overlapPeaks ():
	cmd_rs=['mergePeaks']
    for cmd_r in cmd_rs:
        try:
            subprocess.call([cmd_r,'--help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(colored(cmd_r+
                          ': It is part of Homer or bedtools. It is not installed in your computer or not in the PATH.'+
                          ' Install or copy the executables to the default PATH', 
                          'green', attrs=['bold']))
            exit()

	dirs=[]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

def diffPeaks():






#------ venn diagram of the peaks / overlapping peaks
#_______________________________________
if not os.path.exists(outputdir+'/Peaks/OverlappingPeaks/'):
	os.makedirs(outputdir+'/Peaks/OverlappingPeaks/')
myCleanbedfiles=glob.glob(outputdir + '/Peaks/'+ "*.Clean.bed")
out_dir=outputdir+'/Peaks/OverlappingPeaks/'
mycmd=['mergePeaks', '-d', str(distance)] + myCleanbedfiles + ['-venn', out_dir+'Comparison_venn', '-matrix', out_dir+'Comparison_matrix']
with open(out_dir+'Comparison.txt', "w+") as f:
	run_cmd_file(mycmd,f)
#---- cleaning of the files:1
extensions=['Comparison_venn','Comparison.txt']
for extension in extensions:	
	venn = pd.read_csv(out_dir+extension,sep='\t',header=None)
	venn = venn.astype(str).replace([outputdir+'/Peaks/','.Clean.bed','nan'],"",regex=True)
	venn.to_csv(out_dir+extension,index = None, header=None,sep='\t')
#----
venn = pd.read_csv(out_dir+'Comparison_venn',sep='\t',header=None)
if (venn.shape[1]-2) < 5:
	mycmd=['homer_venn.plot.R','Comparison_venn',out_dir+'Comparison_venn']
	run_cmd(mycmd)
else:
	print('Unable to plot Venn diagram as categories are >5')

#------ Total and differential peaks with annotations (gene ontology, genome ontology and distribution of motifs)
#_______________________________________
out_dir=outputdir+'/Peaks/DifferentialPeaks/'
out_dir_annotation=outputdir+'/Peaks/Peaks_annotations/'
out_motifs=outputdir+'/Peaks/Peaks_DNAmotifs/'
samples = list(myin.loc[:,4])
#-- motifs for total peaks
for sampl in samples:
	mypeak = outputdir + '/Peaks/' + sampl + '.Clean.bed'
	motifs_peaks (genomeversion,mypeak,out_motifs+sampl,size,threads)
#-- annotation for the total peaks
for sampl in samples:
	inpeak = outputdir + '/Peaks/' + sampl + '.Clean.bed'
	outfile = out_dir_annotation +  sampl
	annotation_peaks (inpeak,genomeversion,outfile)
#-- Differential peaks	
for i in combinations(samples,2):
	tag1 = outputdir + '/Tagdirectories/' + i[0] + '_expr'
	tag2 = outputdir + '/Tagdirectories/' + i[1] + '_expr'
	mypeak = outputdir + '/Peaks/' + i[0] + '.Clean.bed'
	differential_peaks(genomeversion,out_dir,out_dir_annotation,mypeak,tag1,tag2,i[0],i[1],pvalueD,foldchangeD)
	motifs_peaks (genomeversion,mypeak,out_motifs+i[0]+'.vs.'+i[1],size,threads)
	inpeak = out_dir+i[0]+'.vs.'+i[1]+'.differentialPeaks.txt'
	outfile = out_dir_annotation+i[0]+'.vs.'+i[1]
	annotation_peaks (inpeak,genomeversion,outfile)
	#---
	tag1 = outputdir + '/Tagdirectories/' + i[1] + '_expr'
	tag2 = outputdir + '/Tagdirectories/' + i[0] + '_expr'
	mypeak = outputdir + '/Peaks/' + i[1] + '.Clean.bed'
	differential_peaks(genomeversion,out_dir,out_dir_annotation,mypeak,tag1,tag2,i[1],i[0],pvalueD,foldchangeD)
	motifs_peaks (genomeversion,mypeak,out_motifs+i[1]+'.vs.'+i[0],size,threads)
	inpeak = out_dir+i[1]+'.vs.'+i[0]+'.differentialPeaks.txt'
	outfile = out_dir_annotation+i[1]+'.vs.'+i[0]
	annotation_peaks (inpeak,genomeversion,outfile)
#------ Really different peaks plot = data stored in XYZ.vs.YZX.differentialPeaks.txt'
#_______________________________________
xx=[]
for i in combinations(samples,2):
	xx.append(i[0]+'.vs.'+i[1])
	xx.append(i[1]+'.vs.'+i[0])
xx=",".join(xx)
mycmd=['Real_differentialPeaks.R','--files',xx,'--dir', out_dir, '-p', str(pvalueD), '-c', str(foldchangeD)]
run_cmd(mycmd)
#------ Matrix: genome ontology
#_______________________________________
matrix_genomeOntology_out=outputdir+'/Peaks/Peaks_matrix_genomeOntology.txt'
xx=[]
for i in combinations(samples,2):
	xx.append(i[0]+'.vs.'+i[1])
	xx.append(i[1]+'.vs.'+i[0])
for i in samples:
	xx.append(i)
matrix_genomeOntology=pd.read_csv(out_dir_annotation+i+ '_GenomeOntology/basic.genomeOntology.txt',sep='\t')
matrix_genomeOntology=matrix_genomeOntology.iloc[:,[0,1]]
for i in xx:
	inputfile=out_dir_annotation+i+ '_GenomeOntology/basic.genomeOntology.txt'
	if os.stat(inputfile).st_size > 0:
		mygenomeontology=pd.read_csv(inputfile,sep='\t')
		new=mygenomeontology.iloc[:,[5,9,10]]
		new=new.add_prefix(i+':')
		matrix_genomeOntology=pd.concat([matrix_genomeOntology, new], axis=1)
matrix_genomeOntology.to_csv(matrix_genomeOntology_out,index=None)
#------ Matrix: motifs
#_______________________________________
matrix_motifs_out=outputdir+'/Peaks/Peaks_matrix_motifs.txt'
xx=[]
for i in combinations(samples,2):
	xx.append(i[0]+'.vs.'+i[1])
	xx.append(i[1]+'.vs.'+i[0])	
for i in samples:
	xx.append(i)
matrix_motifs=pd.read_csv(out_motifs+i+ '/' + 'knownResults.txt',sep='\t')
matrix_motifs=matrix_motifs.iloc[:,0]
for i in xx:
	inputfile=out_motifs+i+ '/' + 'knownResults.txt'
	if os.stat(inputfile).st_size > 0:
		motifs=pd.read_csv(inputfile,sep='\t')
		new=motifs.iloc[:,[4,5,6]]
		new=new.add_prefix(i+':')
		matrix_motifs=pd.concat([matrix_motifs, new], axis=1)
matrix_motifs.to_csv(matrix_motifs_out,index=None)	
#------ Matrix: gene ontology
#_______________________________________		

