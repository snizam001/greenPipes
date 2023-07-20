# download data in current directory and then test it
#_________________________________________

# First download datasets
echo "Download Test Data from https://osf.io/ruhj9"


# Run example dataset

eval "$(conda shell.bash hook)"
conda activate greenpipes

mkdir $(pwd)/Pair-end

greenPipes --modes qc,contamination,qc,contamination,alignment,equalRead,qcExperiment,initPeakCalling,qcTagDirectories,callPeaks,PeakComparison,annotation,UserAnnotation,coverageTracks,initHeatmap,heatmap,piggyBack,cutfrequency --inputdir $(pwd)/Fastq --inputfile $(pwd)/Samplesheet.txt --compareInfile $(pwd)/compare.txt --libraryType pair --outputdir $(pwd)/Pair-end --refgenome $(pwd)/Ref/GenomeBowtei2 --spikein $(pwd)/Ref/SpikeBowtei2 --blackListedRegions $(pwd)/BlackListed.bed --sFasta $(pwd)/Ref/Genome.fa --annFiles $(pwd)/ATACseq.bed,$(pwd)/H3K4me3.bed --annName ATAC,H3K4me3 --genomeFile ./Ref/Genome.genome --covSpike True --covSpike_NormalizationFormula 2 --hCovComp coverage --covExprType gCR --cMotif True --cMaN MA0060 --cGVersion hg38 --cCenter 8 --lProt piggyBack-lProtein.txt --mPwm MA0060 --sProt CCAAT --Species human --gVer hg38  --annSize 400 --hDiffPeaks False

mkdir $(pwd)/Single-end

greenPipes --modes qc,contamination,alignment,equalRead,qcExperiment,initPeakCalling,qcTagDirectories,callPeaks,PeakComparison,annotation,UserAnnotation,coverageTracks,initHeatmap,heatmap,piggyBack,cutfrequency --inputdir $(pwd)/Fastq --inputfile $(pwd)/Samplesheet.txt --compareInfile $(pwd)/compare.txt --libraryType pair --outputdir $(pwd)/Pair-end --refgenome $(pwd)/Ref/GenomeBWA --spikein $(pwd)/Ref/SpikeBWA --blackListedRegions $(pwd)/BlackListed.bed --sFasta $(pwd)/Ref/Genome.fa --annFiles $(pwd)/ATACseq.bed,$(pwd)/H3K4me3.bed --annName ATAC,H3K4me3 --genomeFile ./Ref/Genome.genome --covSpike True --covSpike_NormalizationFormula 2 --hCovComp coverage --covExprType gCR --cMotif True --cMaN MA0060 --cGVersion hg38 --cCenter 8 --lProt piggyBack-lProtein.txt --mPwm MA0060 --sProt CCAAT --Species human --gVer hg38  --annSize 400 --hDiffPeaks False

# Check if everything is ok?
