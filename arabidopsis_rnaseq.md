
# **RNA-Seq analysis for *Arabidopsis thaliana***

*by Ma Lab*


In this tutorial, we will learn how to conduct an RNA-seq data analysis for the model plant species *Arabidopsis thaliana*, using our lab server **Mercury**. We will use a protocol modified from [Pertea *et al.* (2016)](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html).

## **1. Required software**
>- [HiSat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
>- [StringTie](https://ccb.jhu.edu/software/stringtie/)
>- [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html)
>- [NCBI SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)

>**Note:** These softwares have been installed on our sever.


## **2. Overview of the analysis**
>- Build an index for the Arabidopsis genome
>- Download and map the Arabidopsis RNA-Seq reads onto Arabidopsis genome
>- Use StringTie and HTSeq to calculate the expression levels for each gene 

>**Homework for the lab members:** Use [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) or [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to calculate deferentially expressed genes

## **3. Build an index for the *Arabidopsis thaliana* genome**
The HiSat2 program requires an index file for every species in order to do the mapping. Some pre-computed indexes can be downloaded from the program's [website](https://ccb.jhu.edu/software/hisat2/index.shtml). But for plants and many other species, one has to prepare the index files using HiSat2 program. We will build the index files for the Arabidopsis genome. We need these files to build the Arabidopsis genome index.

>- The Arabidopsis genome sequences in FASTA format
>- A [GTF](http://asia.ensembl.org/info/website/upload/gff.html) file specifying the gene models
>
>**Note:** These files will be downloaded from the [Ensembl Plants](http://plants.ensembl.org/index.html) database.

Log into the **Mercury** server, download the the files from Ensembl Plants, and build the genome index. 
```bash
cd ~
# Make the arabidopsis/idx directory where the index files should be saved
mkdir arabidopsis
mkdir arabidopsis/idx
cd arabidopsis/idx

# Download the genome sequences and GTF files from Ensembl
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.34.gtf.gz

# Also copy the E coli sequences. There might be some contamination of e coli sequences during RNA-Seq library preparation. Incorporate these sequences into the index will make sure they will be removed during data analysis.
cp /home/shared/arabidopsis/ERCC* .
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.34.gtf.gz

# Append the e coli sequences and gtf to the Arabidopsis' ones.
cat Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ERCC92.fa >> genome.fa
cat Arabidopsis_thaliana.TAIR10.34.gtf ERCC92.gtf >> genome.gtf

# Extract the exons and splice sites from the gtf file and use them to build the index. Check HiSat2 manual for more details.
python3.4 /home/shared/tools/hisat2-2.0.4/hisat2_extract_exons.py genome.gtf > genome.exon
python3.4 /home/shared/tools/hisat2-2.0.4/hisat2_extract_splice_sites.py genome.gtf > genome.ss

# Build the genome index using hisat2-build
hisat2-build --exon genome.exon --ss genome.ss genome.fa Arabidopsis_thaliana.TAIR10.34_ercc_tran

# Read the HiSat2 manual to find out why should we include the exons and splice sites information in the index.

# Check the results
ls -lh
```
Your should see the following files within the folder:
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.1.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.2.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.3.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.4.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.5.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.6.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.7.ht2
>Arabidopsis_thaliana.TAIR10.34_ercc_tran.8.ht2
>
## **4. Download the RNA-Seq data files from NCBI SRA database**

We will use [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) to download the RNA-Seq data files. First add the directory  /home/shared/tools/sratoolkit.2.8.2-1-centos_linux64/bin to your PATH.
```bash
export PATH=$PATH:/home/shared/tools/sratoolkit.2.8.2-1-centos_linux64/bin
```
We will use a dataset published by [Lin *et al.* (2013)](http://www.plantphysiol.org/content/162/3/1750.full) for our analysis. The dataset, measuring Arabidopsis root's transcriptome upon iron deficiency, were deposited at NCBI SRA database (under Bioproject [PRJNA256121](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA256121)), which include the following samples:
>SRA Run id # | Sample 
>---------- | -----
> SRR1524935 | control root sample #1
> SRR1524938 | control root sample #2
> SRR1524940 | control root sample #3
> SRR1524941 |  iron-deficiency root sample #1
> SRR1524945 | iron-deficiency root sample #2
> SRR1524946 | iron-deficiency root sample #3

Download the run files for all these SRAs *via* the **prefetch** command within the SRA toolkit. **Note:** The files just need to be downloaded once. One of the lab members can do this downloading step. Avoid several people download these files in the same time. 
```bash
# Install the aspera software for fast downloading
/home/shared/tools/aspera-connect-3.6.1.110647-linux-64.sh
# Download the data files for all the runs. 
# Files will be saved at /home/shared/ncbi/public/sra/
prefetch SRR1524935
prefetch SRR1524938
prefetch SRR1524940
prefetch SRR1524941
prefetch SRR1524945
prefetch SRR1524946
# Check the downloaded files
ls /home/shared/ncbi/public
```

## **5. Mapping the files to the genome**

**Fastq-dump** will extract the FASTQ file from the downloaded SRA file. For single-end sequencing, it will generate one FASTQ file, and for paired-end sequencing it will generate two FASTQ files. The files of [PRJNA256121](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA256121) were all produced by single-end sequencing. **hisat2** will map the RNA-Seq reads onto the Arabidopsis genome using the index files we just generated.
```bash
# We will make a directory prjna256121 for the analysis.
cd ~
mkdir ~/arabidopsis/prjna256121
cd ~/arabidopsis/prjna256121
fastq-dump --split-3 SRR1524935
hisat2 -x ~/arabidopsis/idx/Arabidopsis_thaliana.TAIR10.34_ercc_tran -p 10 -U SRR1524935.fastq -S SRR1524935.sam
ls -lh
```

Mapping statistics are shown below. The keys are **reads** (should be >= 1000000) number and **overall alignment rate** (should be >= 70%). 
> 23067621 reads; of these:
>      23067621 (100.00%) were unpaired; of these:
>         658421 (2.85%) aligned 0 times
>         18506106 (80.23%) aligned exactly 1 time
>         3903094 (16.92%) aligned >1 times
>97.15% overall alignment rate

Sort and convert the mapping result sam file into bam file, using **samtools**.
```bash
samtools sort -@ 8 -o SRR1524935.bam SRR1524935.sam
# Keep the bam file, and remove the fastq and sam files
rm SRR1524935.fastq
rm SRR1524935.sam
```

Do the same operations for SRR1524938, SRR1524940, SRR1524941, SRR1524945, and SRR1524946.

## **6. Quantify expressed genes and transcripts' abundance using StringTie**
We are only interested in the known transcripts specified in the Arabidopsis GTF file.
```bash
stringtie SRR1524935.bam -p 8 -G ~/arabidopsis/idx/Arabidopsis_thaliana.TAIR10.34.gtf -e -B -A SRR1524935.gene.txt -o ballgown/SRR1524935/SRR1524935.gtf
# Check StringTie manual for the explanation of each parameter
cat SRR1524935.gene.txt | more
```
The resulted file SRR1524935.gene.txt contains the expression value for each gene, in the format of FPKM or TPM.
>**FPKM:** Fragments Per Kilobase Million
>**TPM:** Transcripts Per Kilobase Million
>
> **Note**: see [a comparison between FPKM and TPM](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/).

Do the same operations for SRR1524938, SRR1524940, SRR1524941, SRR1524945, and SRR1524946.

## **7. Quantify gene expression values using HTSeq**

```shell
samtools view SRR1524935.bam | htseq-count -m union -s no - ~/arabidopsis/idx/Arabidopsis_thaliana.TAIR10.34.gtf > SRR1524935.gene.count.txt
```
Check the result file SRR1524935.gene.count.txt for the read counts for each gene. Do the same operation for SRR1524938, SRR1524940, SRR1524941, SRR1524945, and SRR1524946.


### **Homework to do:** 
> Calculate the deferentially expressed genes between the iron-deficiency samples and the control samples using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) or [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).  These two programs use the gene count numbers as input, which were generated *via* the HTSeq program. You can also use the [ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html) program to find out induced and repressed genes.


### **References:**

Li W, Lin WD, Ray P, Lan P, Schmidt W: Genome-wide detection of condition-sensitive alternative splicing in Arabidopsis roots. *Plant Physiol* 2013, **162**:1750-1763.

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL:  Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. _Nat Protoc_ 2016, **11**: 1650-1667.


