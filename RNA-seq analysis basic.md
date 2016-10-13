In this tutorial, we will learn how to do a basic RNA-seq analysis. We will first setup a linux server using cloud service provided by USTC. Then will follow a protocol provided by [Pertea et al (2016)](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html), which appeared at **_Nature Protocol_**:

`RNA-seq analysis using HISAT, StringTie and Ballgown, using a protocol published in the following paper at Nature Protocol.`


##1. Setup and virtual machine

Open an account at USTC cloud service, and set up a virtual machine as a server:

http://cloud.ustc.edu.cn/

Connect to the virtual machine using Xshell

##2. Software installation.

```bash
mkdir $HOME/bin
export PATH=$HOME/bin:$PATH
sudo yum groupinstall 'Development Tools'
cd ~
mkdir tools
cd tools
```

```linux
wget http://staff.ustc.edu.cn/~sma/bioinfo_tools/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
tar –xjf samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
cp samtools-bcftools-htslib-1.0_x64-linux/bin/samtools ~/bin

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
cp hisat2-2.0.4/hisat2* ~/bin
cp hisat2-2.0.4/*py ~/bin

wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
tar xvfz stringtie-1.3.0.Linux_x86_64.tar.gz
cp stringtie-1.3.0.Linux_x86_64/stringtie ~/bin
```
```linux
sudo yum install openssl-devel
sudo yum install libcurl-devel
sudo yum install httr
sudo yum install libxml2-devel
sudo yum install screen
```

Intall R packages

```linux
sudo R
```

```R
install.packages("devtools",repos="http://cran.us.r-project.org")
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown", "genefilter","dplyr","devtools"))
q()
```

##3. Download and unzip the data

```linux
screen
cd ~
mkdir my_rnaseq_exp
cd my_rnaseq_exp
wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
tar xvzf chrX_data.tar.gz
```

##4. Analysis

