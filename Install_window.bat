@echo off
setlocal

REM Install mamba and create greenpipes environment
REM --

set "currentDirectory=%cd%"

where mamba > nul
if %errorlevel% neq 0 (
    echo ----------------------------------------------------------------
    echo mamba is not present in the system. Installing it before proceeding! You can install it by installing mambaforge from https://github.com/conda-forge/miniforge for your operating system.
    echo ----------------------------------------------------------------
    exit /b 0
    REM REM architecture=$(uname -m)
    REM REM wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-%architecture%.sh"
    REM REM conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'
)

REM create greenpipes environment
for /f "tokens=*" %%A in ('conda shell.bash hook') do %%A
call conda activate base

conda info --envs | find "greenpipes" > nul
if %errorlevel% equ 0 (
    echo ---------------------------------------
    echo [greenpipes] environment already exists
    echo ---------------------------------------
) else (
    echo ---------------------------------------
    echo creating python environment: greenpipes ...
    echo ---------------------------------------
    call mamba env create --name greenpipes -f %currentDirectory%\environment.yaml
)

REM Install greenpipes

REM  - bash
REM  - fish
REM  - tcsh
REM  - xonsh
REM  - zsh
REM  - powershell

for /f "tokens=*" %%A in ('conda shell.bash hook') do %%A
call conda activate greenpipes

where greenPipes > nul
if %errorlevel% neq 0 (
    echo ---------------------------------------
    echo Installing greenPipes ...
    echo ---------------------------------------
    attrib +x %currentDirectory%\greenPipes\rscripts\*
    pip install %currentDirectory%\
) else (
    echo ---------------------------------------
    echo greenPipes is already installed ...
    echo ---------------------------------------
)

REM Download the homer package database
echo ---------------------------------------
echo Installing database for your organism of interest.
echo ______________________________________

for /f "tokens=*" %%A in ('conda config --show envs_dirs ^| findstr /c:" - " ^| %SystemRoot%\System32\findstr.exe /r /n ".*"') do (
    set "mydirectory=%%A"
)
set "mydirectory=%mydirectory:~3%"

perl %mydirectory%\greenpipes\share\homer\configureHomer.pl -list
echo ______________________________________
echo.
echo Choose your organism from the above Genome list. For example, in case of human, you can use hg38.

set /p "organismName=genome = "
perl %mydirectory%\greenpipes\share\homer\configureHomer.pl -install %organismName%

REM Install fastq-screen package database
for /f "tokens=*" %%A in ('conda shell.bash hook') do %%A
call conda activate greenpipes

echo Installing the perl module GD and installing genomes for the fastq_screen
call conda install -y GD
mkdir %mydirectory%\greenpipes\fastq_screenData
fastq_screen --get_genomes --outdir %mydirectory%\greenpipes\fastq_screenData
copy %mydirectory%\greenpipes\fastq_screenData\FastQ_Screen_Genomes\fastq_screen.conf %mydirectory%\greenpipes\share\fastq-screen-0.15.3-0

REM Install perl library for the meme suites
for /f "tokens=*" %%A in ('conda shell.bash hook') do %%A
call conda activate greenpipes

where make > nul
if %errorlevel% neq 0 (
    echo make is not available in your system. Maybe gcc might be also absent in this system.
    echo Use: sudo apt-get install build-essential in Ubuntu
    exit /b 1
) else (
    del %currentDirectory%\XML-Parser-2.46.tar.gz /f
    echo Installing different perl modules for the MEME-suite. These dependencies were not installed by CONDA
    cpanm install Cwd File::Which Data::Dumper Exporter Fcntl File::Basename File::Copy File::Path File::Spec::Functions File::Temp Getopt::Long HTML::PullParser HTML::Template HTML::TreeBuilder JSON List::Util Pod::Usage POSIX Scalar::Util XML::Simple Sys::Info Log::Log4perl Math::CDF Sys::Hostname Time::HiRes XML::Compile::SOAP11 XML::Compile::WSDL11 XML::Compile::Transport::SOAPHTTP
    set "mydirectory=%mydirectory:~3%"
    powershell -Command "(New-Object Net.WebClient).DownloadFile('http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz', '%currentDirectory%\XML-Parser-2.46.tar.gz')"
    tar -zxvf %currentDirectory%\XML-Parser-2.46.tar.gz -C %currentDirectory%
    cd %currentDirectory%\XML-Parser-2.46
    perl Makefile.PL EXPATINCPATH=%mydirectory%\greenpipes\include\ EXPATLIBPATH=%mydirectory%\greenpipes\lib\
    make
    make test
    make install
)

REM Installing R packages
echo Installing the different R libraries.
cd %currentDirectory%
Rscript %currentDirectory%\Install_r.R
