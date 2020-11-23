# Practicals NGS Bats 
## Objectives - Characterization of the IFN response in bats
The objective is to investigate the specificities of the bats immune system, in order to understand how bat's immunity has adapted to long-term association with viruses, and why bats harbor a high number of viruses but appear mostly asymptomatic. 
To do so, we want to characterize the interferon (IFN) response in bats. We thus aim at studying the transcriptomic changes in response to IFN stimulation in bat cells, and in particular, we want to identify which interferon-stimulated genes (ISGs) are differentially expressed after IFN stimulation. 

We focus on *Myotis velifer*, a bat species for which no transcriptomic data was available yet.  

### Step 0: Data obtention
RNA-seq was extracted from bat cells exposed or not to IFN. 
Libraries were prepared and paired-end RNA-seq was performed using the Illumina NextSeq500 platform. 3 samples of each condition were sequenced.

Download the raw RNA-seq data from the IGFL sequencing platform using the command line in the script **importdataset.sh**. The sequences are under fastq format. 

### Step 1: De novo assembly 
Th first step is to assemble the transcriptome of *Myotis velifer*. Because no reference transcriptome is available for this species, we perform a *de novo* transcriptome assembly. We use all 6 sequenced samples together to assemble the transcriptome. 

- **Quality control:** We first need control the quality of the reads. We run fastqc on each of the obtained sequences, using the script **qualitycontrol_fastqc.sh**. Fastqc analyzes the quality of sequences using several statistics. Here is an overview of the fastqc report for representative sequence. 
![fastqc report](/home/rstudio/data/mydatalocal/imagesreadme/report_fastqc.png)




-  **Read trimming:** Trimmomatic **trimmomatic.sh**. Then verification with fastqc again **control_fastqc_after_trimm.sh**


- **_De novo_ assembly:** Trinity **assembly_trinity.sh**

### Step 2: Data annotation
???

### Step 3: Transcript expression quantification

- **Creating an index with salmon:** salmon index (first part of the **salmon_quantif.sh** script)

- **Quantification:** salmon quant (second part of the **salmon_quantif.sh** script)

### Step 4: Differential expression analyses
