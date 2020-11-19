# Practicals NGS Bats 
## Objectives - Characterization of the IFN response in bats
The objective is to investigate the specificities of the bats immune system, in order to understand how bat's immuntity has adapted to long-term association with viruses, and why bats harbor a high number of viruses but appear mostly asymptomatic. 
To do so, we want to characterize the interferon (IFN) response in bats. We thus aim at studying the transcriptomic changes in response to IFN stimulation in bat cells, and in particular, we want to identify which interferon-stimulated genes (ISGs) are differentially expressed after IFN stimulation. 

We focus on *Myotis velifer*, a bat species for which transcriptomic data was not available yet.  

### Step 0: Data obtention
RNA-seq was extracted from bat cells exposed or not to IFN. 
Libraries were prepared and paired-end RNA-seq was performed using the Illumina NextSeq500 platform. 3 samples of each condition were sequenced.

Download the raw RNA-seq data from the IGFL sequencing platform using the command line in the script **importdataset.sh**.

### Step 1: De novo assembly 
Th first step is to assemble the transcriptome of *Myotis velifer*. Because no reference transcriptome is available for this species, we do a *de novo* transcriptome assembly. We use all 6 sequenced samples together to assemble a transcriptome. 

- **Quality control:** fastqc **qualitycontrol_fastqc.sh**

-  **Read trimming:** Trimmomatic **trimmomatic.sh**. Then verification with fastqc again **control_fastqc_after_trimm.sh**


- **_De novo_ assembly:** Trinity **assembly_trinity.sh**

### Step 2: Data annotation
???

### Step 3: Transcript expression quantification

- **Creating an index with salmon:** salmon index (first part of the **salmon_quantif.sh** script)

- **Quantification:** salmon quant (second part of the **salmon_quantif.sh** script)

### Step 4: Differential expression analyses
