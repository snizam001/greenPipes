#!/bin/bash

# checking if system connected with internet
wget -q --spider http://google.com

if [ $? -eq 0 ]; then
    echo ""
else
    echo "Error: system is not connected with Internet...."
		exit 0
fi

# Install mamba and create greenpipes environment
# --
currentDirectory=$(pwd)
if ! [ -x "$(command -v mamba)" ]
	then
	echo "----------------------------------------------------------------"
	echo "mamba is not present in system. Installing it before proceeding! You can install it by installing mambaforge from https://github.com/conda-forge/miniforge for your operating system."
	echo "----------------------------------------------------------------"
	exit 0
fi

# create greenpipes environment
eval "$(conda shell.bash hook)"
conda activate base

if conda info --envs | grep -q greenpipes -w
	then
		echo "---------------------------------------"
		echo "[greenpipes] environment already exists. Uninstall it before running this script by using following command:"
		echo "mamba env remove --name greenpipes"
		echo "---------------------------------------"
else
	echo "---------------------------------------"
	echo "creating python environment: greenpipes ..."
	echo "---------------------------------------"
	mamba env  create --name greenpipes -f $(pwd)/environment.yaml
fi

eval "$(conda shell.bash hook)"
conda activate greenpipes

if ! [ -x "$(command -v greenPipes)" ]
	then
	echo "---------------------------------------"
	echo "Installing greenPipes ..."
	echo "---------------------------------------"
	chmod -R +x $(pwd)/greenPipes/*
	pip install $(pwd)/
else
	echo "---------------------------------------"
	echo "greenPipes is already installed. Uninstall it before running using following command:"
	echo "pip uninstall greenPipes --yes"
	echo "---------------------------------------"
fi

# Download the homer package database
echo "---------------------------------------"
echo "Installing database for your organism of interest. "
echo "______________________________________"

mydirectory=$(conda config --show envs_dirs |sed -n 2p | sed 's/  - //g')

$mydirectory/greenpipes/bin/perl $mydirectory/greenpipes/share/homer/configureHomer.pl -list
echo "______________________________________
.
.
.
.
Choose your organism from the above Genome list. For example in case of human you can use hg38."

echo -n "genome = "
read organismName
$mydirectory/greenpipes/bin/perl $mydirectory/greenpipes/share/homer/configureHomer.pl -install $organismName

# Install fastq-screen package database
eval "$(conda shell.bash hook)"
conda activate greenpipes

echo "Installing the perl module GD and installing genomes for the fastq_screen"
$mydirectory/greenpipes/bin/cpanm install GD
mkdir $mydirectory/greenpipes/fastq_screenData
$mydirectory/greenpipes/bin/fastq_screen --get_genomes --outdir $mydirectory/greenpipes/fastq_screenData
cp  $mydirectory/greenpipes/fastq_screenData/FastQ_Screen_Genomes/fastq_screen.conf $mydirectory/greenpipes/share/fastq-screen-0.15.3-0

# Install perl library for the meme suites
eval "$(conda shell.bash hook)"
conda activate greenpipes

if ! [ -x "$(command -v make)" ]
	then
		echo "make is not available in your system. May be gcc might be also absent in this system.
		Use: sudo apt-get install build-essential in ubuntu"
		exit 1
	else
		rm XML-Parser-2.46.tar.gz
		echo "Installing different perl modules for the MEME-suite. These dependencies were not installed by the CONDA"
		$mydirectory/greenpipes/bin/cpanm install Cwd File::Which Data::Dumper Exporter Fcntl File::Basename  File::Copy  File::Path  File::Spec::Functions File::Temp  Getopt::Long  HTML::PullParser HTML::Template HTML::TreeBuilder JSON List::Util  Pod::Usage POSIX Scalar::Util  XML::Simple Sys::Info  Log::Log4perl  Math::CDF  Sys::Hostname Time::HiRes XML::Compile::SOAP11 XML::Compile::WSDL11 XML::Compile::Transport::SOAPHTTP
		mydirectory=$(conda config --show envs_dirs |sed -n 2p | sed 's/  - //g')
		wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz
		tar -zxvf XML-Parser-2.46.tar.gz
		cd XML-Parser-2.46
		$mydirectory/greenpipes/bin/perl Makefile.PL EXPATINCPATH=$mydirectory/greenpipes/include/ EXPATLIBPATH=$mydirectory/greenpipes/lib/
		make
		make test
		make install
fi

# Installing R packages
# Install perl library for the meme suites
eval "$(conda shell.bash hook)"
conda activate greenpipes

echo "Installing the different R libraries."
cd $currentDirectory
$mydirectory/greenpipes/bin/Rscript $(pwd)/Install_r.R
