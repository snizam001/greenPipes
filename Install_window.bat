@echo off
REM Install mamba and create greenpipes environment
REM --

IF NOT EXIST "%CONDA_PREFIX%\Scripts\mamba.exe" (
    echo ----------------------------------------------------------------
    echo "mamba is not present in system. Installing it before proceeding! You can install it by installing mambaforge from https://github.com/conda-forge/miniforge for your operating system."
    echo ----------------------------------------------------------------
    exit /b 0
)

REM Create greenpipes environment
call conda activate base

IF EXIST "%CONDA_PREFIX%\envs\greenpipes" (
    echo ---------------------------------------
    echo [greenpipes] environment already exists
    echo ---------------------------------------
) ELSE (
    echo ---------------------------------------
    echo creating python environment: greenpipes ...
    echo ---------------------------------------
    call mamba env create --name greenpipes -f "%CD%\environment.yml"
)

REM Install greenpipes
call conda activate greenpipes

IF NOT EXIST "%CONDA_PREFIX%\Scripts\greenPipes.exe" (
    echo ---------------------------------------
    echo Installing greenPipes ...
    echo ---------------------------------------
    chmod +x "%CD%\greenPipes\rscripts\*"
    pip install "%CD%\"
) ELSE (
    echo ---------------------------------------
    echo greenPipes is already installed ...
    echo ---------------------------------------
)

REM Download the homer package database
echo ---------------------------------------
echo Installing database for your organism of interest.
echo ---------------------------------------

FOR /F "usebackq tokens=*" %%A IN (`conda config --show envs_dirs ^| find /I "envs"`) DO SET "mydirectory=%%A"
call perl "%mydirectory%\greenpipes\share\homer\configureHomer.pl" -list
echo ---------------------------------------
...
...
...
Choose your organism from the above Genome list. For example, in the case of human, you can use hg38.

set /p "organismName=genome = "
call perl "%mydirectory%\greenpipes\share\homer\configureHomer.pl" -install %organismName%

REM Install fastq-screen package database
call conda activate greenpipes

echo Installing the perl module GD and installing genomes for the fastq_screen
call cpnam install GD
mkdir "%mydirectory%\greenpipes\fastq_screenData"
fastq_screen --get_genomes --outdir "%mydirectory%\greenpipes\fastq_screenData"
copy "%mydirectory%\greenpipes\fastq_screenData\FastQ_Screen_Genomes\fastq_screen.conf" "%mydirectory%\greenpipes\share\fastq-screen-0.15.3-0"

REM Install perl library for the meme suites
call conda activate greenpipes

IF NOT EXIST "%CONDA_PREFIX%\Library\bin\make.exe" (
    echo make is not available in your system. May be gcc might be also absent in this system.
    echo Use: sudo apt-get install build-essential in ubuntu
    exit /b 1
) ELSE (
    del XML-Parser-2.46.tar.gz
    echo Installing different perl modules for the MEME-suite. These dependencies were not installed by the CONDA
    call cpanm install Cwd File::Which Data::Dumper Exporter Fcntl File::Basename  File::Copy  File::Path  File::Spec::Functions File::Temp  Getopt::Long  HTML::PullParser HTML::Template HTML::TreeBuilder JSON List::Util  Pod::Usage POSIX Scalar::Util  XML::Simple Sys::Info  Log::Log4perl  Math::CDF  Sys::Hostname Time::HiRes XML::Compile::SOAP11 XML::Compile::WSDL11 XML::Compile::Transport::SOAPHTTP
    SET "mydirectory=%CONDA_PREFIX%"
    wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz
    tar -zxvf XML-Parser-2.46.tar.gz
    cd XML-Parser-2.46
    perl Makefile.PL EXPATINCPATH="%mydirectory%\greenpipes\include\" EXPATLIBPATH="%mydirectory%\greenpipes\lib\"
    make
    make test
    make install
)

REM Installing R packages
echo Installing the different R libraries.
Rscript "%CD%\Install_r.R"
