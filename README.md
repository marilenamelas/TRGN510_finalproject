# TRGN510_finalproject
 The purpose of my final project was to perform RNAseq analysis for an individual with clear cell renal cell carcinoma, using three FFPE derived tumor samples and one adjacent normal. My samples are labelled as : Normal, Tumor1, Tumor2 and Tumor3 and all come from the same individual.  
 
PROJECT GOAL/DESCRIPTION:
RNAseq analysis (paired-end sequencing in two lanes) was performed using three FFPE derived clear cell renal cell carcinoma tumor samples VS one FFPE derived adjacent normal. Ultimate scope: QC of the samples and differential expression using the adjacent normal sample as my control. I worked both using Bash/command line and R script. The steps for my analysis are described below in detail: 
 
Work performed on the ITG server:
1) Alignment of FASTQ files: 
a) STAR aligner (reference genome hg19 was downloaded from UCSC, and the appropriate Index was built)
b) Convertion of SAM to BAM files. 
c) Convertion of BAM files to sorted_Bam.
d) Creation of Index (BAI) using the sorted BAM files.
e) Loading of BAM files on IGV to check the expression of some well known genes in the tissue.

2) Transcripts quantification:
I created the transcriptome Index using Salmon (downloaded Homo_sapiens.GRCh37.67.cdna from Ensembl). 
When the quantification was complete, I got the quant.sf files which contained the quantification results for each of my samples. 

Work performed on R studio using Bioconductor packages: 

I first used the Bioconductor package 'tximport' which imported and summarized trascript-level estimates for transcript- and gene-level analysis. I was able to create the tx2gene (tx_id, gene_id) using annotation from Ensembl (packages used: EnsDb.Hsapiens.v75, ensembldb). A DGEList object which I named 'y', was created using the 'EdgeR' package. 
Then, I did the filtering retaining only those genes which are represented at least 1cpm (counts per million) least in at least two samples. Then, I performed normalization using the trimmed mean of M-values (TM) method. Then, I set up the model (designed the matrix) and estimated the Dispersions (common dispersion, trended dispersion and tagwise dispersion). I created two groups: one with the normal sample and a second with the three tumors. From an initial QC i did (creating a plot based on multidimensional scaling), i noticed that tumor1 and tumor2 are grouped/clustered together very closely, but tumor3 is far away (which I was expecting based on the bad sequencing quality of tumor3).  As an initial attempt, I decided to perform a 'quick and dirty' analysis with all my tumors together VS the normal control. For the differential expression, I used the TopTags function to explore the results and set thresholds to identify subsets of differentially expressed genes. Finally, 
I ploted the log-fold changes of all the genes (MA plot) and highlighted those that are differentially expressed in red.
 

FILES ORGANIZATION/LOCATION:

A) All my files are inside the directory called 'finalproject' on melas@itg.usc.edu server. I have organized my fastq raw data inside different subdirectories called, a) Fastq_normal, b) Fastq_tumor_1. c) Fastq_tumor_2, and Fastq_tumor_3. I first unzipped them using 'gunzip' command.The 'quant.sf' files with the transcript quantification metrics are inside the directory 'quants' within each of the subdirectories afore mentioned for each sample individually. The detailed code/scripts which I used for the STAR alignment and transcript quantificaion using Salmon is described in the R Markdown document uploaded on github

B) Where to find my RDATA:
Inside the directory 'finalproject' there is another directory called 'R_data'. This includes a subdirectory named    'My_work_in_R' which contains two files: a) Mywork_Rdata.rmd and Tables.RData. My detailed R script is described in the R script document uploaded on github.
 
 
 BUILD WITH:
 1. Bash/command line on itg server (melas@itg.usc.edu)
 2. R script (Bioconductor packages)
 3. R markdown 
 
LICENSE:
This project is licensed under the MIT License - see the LICENSE file for details
 
ACKNOWLEDGEMENTS:
 
 Marilena Melas 
