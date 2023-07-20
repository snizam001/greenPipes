# Install mamba and create greenpipes environment
# --

if (-Not (Get-Command mamba -ErrorAction SilentlyContinue))
{
    Write-Host "----------------------------------------------------------------"
    Write-Host "mamba is not present in system. Installing it before proceeding! You can install it by installing mambaforge from https://github.com/conda-forge/miniforge for your operating system."
    Write-Host "----------------------------------------------------------------"
    Exit 0
    # You can download and install mamba here for Windows, similar to the commented out section in the shell script
}

# create greenpipes environment
conda activate base

if (conda info --envs | Select-String -Pattern 'greenpipes' -Quiet)
{
    Write-Host "---------------------------------------"
    Write-Host "[greenpipes] environment already exists"
    Write-Host "---------------------------------------"
}
else
{
    Write-Host "---------------------------------------"
    Write-Host "creating python environment: greenpipes ..."
    Write-Host "---------------------------------------"
    mamba env create --name greenpipes -f "$pwd/environment.yml"
}

conda activate greenpipes

if (-Not (Get-Command greenPipes -ErrorAction SilentlyContinue))
{
    Write-Host "---------------------------------------"
    Write-Host "Installing greenPipes ..."
    Write-Host "---------------------------------------"
    Get-ChildItem "$pwd/greenPipes/rscripts/*" | ForEach-Object { $_.IsReadOnly = $false }
    pip install $pwd\
}
else
{
    Write-Host "---------------------------------------"
    Write-Host "greenPipes is already installed ..."
    Write-Host "---------------------------------------"
}

# Download the homer package database
Write-Host "---------------------------------------"
Write-Host "Installing database for your organism of interest. "
Write-Host "______________________________________"

$mydirectory = (conda config --show envs_dirs).Split(" ")[1]

perl "$mydirectory/greenpipes/share/homer/configureHomer.pl" -list
Write-Host "______________________________________
.
.
.
.
Choose your organism from the above Genome list. For example in case of human you can use hg38."

Write-Host -NoNewline "genome = "
$organismName = Read-Host
perl "$mydirectory/greenpipes/share/homer/configureHomer.pl" -install $organismName

# Install fastq-screen package database
conda activate greenpipes

Write-Host "Installing the perl module GD and installing genomes for the fastq_screen"
cpnam install GD
mkdir "$mydirectory/greenpipes/fastq_screenData"
fastq_screen --get_genomes --outdir "$mydirectory/greenpipes/fastq_screenData"
Copy-Item "$mydirectory/greenpipes/fastq_screenData/FastQ_Screen_Genomes/fastq_screen.conf" "$mydirectory/greenpipes/share/fastq-screen-0.15.3-0"

# Install perl library for the meme suites
conda activate greenpipes

if (-Not (Get-Command make -ErrorAction SilentlyContinue))
{
    Write-Host "make is not available in your system. May be gcc might be also absent in this system."
    Write-Host "Use: sudo apt-get install build-essential in ubuntu"
    Exit 1
}
else
{
    Remove-Item XML-Parser-2.46.tar.gz
    Write-Host "Installing different perl modules for the MEME-suite. These dependencies were not installed by the CONDA"
    cpanm install Cwd File::Which Data::Dumper Exporter Fcntl File::Basename  File::Copy  File::Path  File::Spec::Functions File::Temp  Getopt::Long  HTML::PullParser HTML::Template HTML::TreeBuilder JSON List::Util  Pod::Usage POSIX Scalar::Util  XML::Simple Sys::Info  Log::Log4perl  Math::CDF  Sys::Hostname Time::HiRes XML::Compile::SOAP11 XML::Compile::WSDL11 XML::Compile::Transport::SOAPHTTP
    $mydirectory = (conda config --show envs_dirs).Split(" ")[1]
    wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz
    tar -zxvf XML-Parser-2.46.tar.gz
    Set-Location XML-Parser-2.46
    perl Makefile.PL EXPATINCPATH="$mydirectory/greenpipes/include/" EXPATLIBPATH="$mydirectory/greenpipes/lib/"
    make
    make test
    make install
}

# Installing R packages
Write-Host "Installing the different R libraries."
Rscript "$pwd/Install_r.R"
