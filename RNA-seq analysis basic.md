In this tutorial, we will learn how to conduct a basic RNA-seq data analysis. We will first setup a virtual machine with Centos 6.5 using the cloud service provided by USTC. Then we will follow a protocol provided by [Pertea et al (2016)](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html), which appeared at **_Nature Protocol_**:

##1. Setup a virtual machine

Open an account at USTC cloud service, and set up a virtual machine as a server:

http://cloud.ustc.edu.cn/

Connect to the virtual machine using [Xshell](http://scc.ustc.edu.cn/yhsq/dlrjxz/201303/t20130314_148192.html).

##2. Software installation.


Onced connected, type the following commands in the Xshell console.


```shell
sudo yum groupinstall 'Development Tools'
sudo yum install openssl-devel
sudo yum install libcurl-devel
sudo yum install httr
sudo yum install libxml2-devel
sudo yum install R
```
**sudo** will grant the root previledge to install softwares in the root directory.
This will install the system tools we need for the softwares we'll use. 


```shell
sudo yum install screen
```
**screen** is a tool to offers the ability to detach a long running process from a session and then attach it back at a later time. You should definitedly check it out. Here is a [tutorial](http://www.cnblogs.com/mchina/archive/2013/01/30/2880680.html) for it.


Next, we will setup two directory for bin (excutable programs) and tools (various bioinformatic tools)
```bash
mkdir $HOME/bin
export PATH=$HOME/bin:$PATH # This will let you run the programs in the bin from any other directories.

```

We will download the pre-compiled Samtools, Hisat2, and Stringtie, unzip, and copy them into the bin directory.

```shell

cd ~
mkdir tools
cd ~/tools

wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar –xjf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --without-curses
make
make install
cp samtools ~/bin


cd ~/tools
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
cp hisat2-2.0.4/hisat2* ~/bin
cp hisat2-2.0.4/*py ~/bin

wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
tar xvfz stringtie-1.3.0.Linux_x86_64.tar.gz
cp stringtie-1.3.0.Linux_x86_64/stringtie ~/bin
```
Note the we use the "~" sign as a representation of the home directory. 


Start R from the linux console and install the necessory R packages within the R console.

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
We can the detach from the screen by typing "Ctrl+A" and then "D". And we can attach to the screen again by typing 

```bash
screen -r
```

##4. Analysis

**Follow the [procedure section](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html#procedure) in the Nature Protocol paper.**

__Map the reads to the genome.__
```bash
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz -2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188104_chrX_1.fastq.gz -2 chrX_data/samples/ERR188104_chrX_2.fastq.gz -S ERR188104_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188234_chrX_1.fastq.gz -2 chrX_data/samples/ERR188234_chrX_2.fastq.gz -S ERR188234_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188245_chrX_1.fastq.gz -2 chrX_data/samples/ERR188245_chrX_2.fastq.gz -S ERR188245_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188257_chrX_1.fastq.gz -2 chrX_data/samples/ERR188257_chrX_2.fastq.gz -S ERR188257_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188273_chrX_1.fastq.gz -2 chrX_data/samples/ERR188273_chrX_2.fastq.gz -S ERR188273_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188337_chrX_1.fastq.gz -2 chrX_data/samples/ERR188337_chrX_2.fastq.gz -S ERR188337_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188383_chrX_1.fastq.gz -2 chrX_data/samples/ERR188383_chrX_2.fastq.gz -S ERR188383_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188401_chrX_1.fastq.gz -2 chrX_data/samples/ERR188401_chrX_2.fastq.gz -S ERR188401_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz -2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188454_chrX_1.fastq.gz -2 chrX_data/samples/ERR188454_chrX_2.fastq.gz -S ERR188454_chrX.sam
hisat2 -p 2 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR204916_chrX_1.fastq.gz -2 chrX_data/samples/ERR204916_chrX_2.fastq.gz -S ERR204916_chrX.sam
```
__Sort and convert the SAM files to BAM__

__To be continued.__


###Reference:

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL. 2016. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. _Nat Protoc_ **11**: 1650-1667.


