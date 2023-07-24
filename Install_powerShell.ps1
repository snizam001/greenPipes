# Install mamba and create greenpipes environment
# --
$currentDirectory = Get-Location

if (-not (Get-Command mamba -ErrorAction SilentlyContinue)) {
    Write-Host "----------------------------------------------------------------"
    Write-Host "mamba is not present in the system. Installing it before proceeding! You can install it by installing mambaforge from https://github.com/conda-forge/miniforge for your operating system."
    Write-Host "----------------------------------------------------------------"
    exit 0
    # # architecture=$(uname -m)
    # # wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$architecture.sh"
    # # conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'
}

# create greenpipes environment
Invoke-Expression (& conda shell.bash hook)
conda activate base

if ((conda info --envs | Select-String -Pattern "greenpipes" -Quiet)) {
    Write-Host "---------------------------------------"
    Write-Host "[greenpipes] environment already exists"
    Write-Host "---------------------------------------"
} else {
    Write-Host "---------------------------------------"
    Write-Host "creating python environment: greenpipes ..."
    Write-Host "---------------------------------------"
    mamba env create --name greenpipes -f $($currentDirectory.Path + "\environment.yaml")
}

# Install greenpipes

#  - bash
#  - fish
#  - tcsh
#  - xonsh
#  - zsh
#  - powershell

Invoke-Expression (& conda shell.bash hook)
conda activate greenpipes

if (-not (Get-Command greenPipes -ErrorAction SilentlyContinue)) {
    Write-Host "---------------------------------------"
    Write-Host "Installing greenPipes ..."
    Write-Host "---------------------------------------"
    chmod +x $($currentDirectory.Path + "\greenPipes\rscripts\*")
    pip install $($currentDirectory.Path + "\")
} else {
    Write-Host "---------------------------------------"
    Write-Host "greenPipes is already installed ..."
    Write-Host "---------------------------------------"
}

# Download the homer package database
Write-Host "---------------------------------------"
Write-Host "Installing database for your organism of interest. "
Write-Host "______________________________________"

$mydirectory = (conda config --show envs_dirs | Select-String -Pattern "- " -Quiet).ToString().TrimStart("- ")

perl "$mydirectory\greenpipes\share\homer\configureHomer.pl" -list
Write-Host "______________________________________"
Write-Host ""
Write-Host "Choose your organism from the above Genome list. For example, in case of human, you can use hg38."

$organismName = Read-Host -Prompt "genome = "
perl "$mydirectory\greenpipes\share\homer\configureHomer.pl" -install $organismName

# Install fastq-screen package database
Invoke-Expression (& conda shell.bash hook)
conda activate greenpipes

Write-Host "Installing the perl module GD and installing genomes for the fastq_screen"
cpnam install GD
mkdir "$mydirectory\greenpipes\fastq_screenData"
fastq_screen --get_genomes --outdir "$mydirectory\greenpipes\fastq_screenData"
cp  "$mydirectory\greenpipes\fastq_screenData\FastQ_Screen_Genomes\fastq_screen.conf" "$mydirectory\greenpipes\share\fastq-screen-0.15.3-0"

# Install perl library for the meme suites
Invoke-Expression (& conda shell.bash hook)
conda activate greenpipes

if (-not (Get-Command make -ErrorAction SilentlyContinue)) {
    Write-Host "make is not available in your system. Maybe gcc might be also absent in this system."
    Write-Host "Use: sudo apt-get install build-essential in Ubuntu"
    exit 1
} else {
    Remove-Item "$($currentDirectory.Path)\XML-Parser-2.46.tar.gz" -Force
    Write-Host "Installing different perl modules for the MEME-suite. These dependencies were not installed by CONDA"
    cpanm install Cwd File::Which Data::Dumper Exporter Fcntl File::Basename File::Copy File::Path File::Spec::Functions File::Temp Getopt::Long HTML::PullParser HTML::Template HTML::TreeBuilder JSON List::Util Pod::Usage POSIX Scalar::Util XML::Simple Sys::Info Log::Log4perl Math::CDF Sys::Hostname Time::HiRes XML::Compile::SOAP11 XML::Compile::WSDL11 XML::Compile::Transport::SOAPHTTP
    $mydirectory = (conda config --show envs_dirs | Select-String -Pattern "- " -Quiet).ToString().TrimStart("- ")
    (New-Object System.Net.WebClient).DownloadFile("http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz", "$($currentDirectory.Path)\XML-Parser-2.46.tar.gz")
    Expand-Archive "$($currentDirectory.Path)\XML-Parser-2.46.tar.gz" "$($currentDirectory.Path)\XML-Parser-2.46"
    Set-Location "$($currentDirectory.Path)\XML-Parser-2.46"
    perl Makefile.PL EXPATINCPATH="$($mydirectory)\greenpipes\include\" EXPATLIBPATH="$($mydirectory)\greenpipes\lib\"
    make
    make test
    make install
}

# Installing R packages

Write-Host "Installing the different R libraries."
Set-Location $currentDirectory
Rscript "$($currentDirectory.Path)\Install_r.R"
