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
The first step is to assemble the transcriptome of *Myotis velifer*. Because no reference transcriptome is available for this species, we perform a *de novo* transcriptome assembly. 

**1) Quality control:** 

We first need control the quality of the reads. We run fastqc on each of the obtained sequences, using the script **qualitycontrol_fastqc.sh**. Fastqc uses several statistics to analyze the quality of sequences. Here is an overview of the fastqc report obtained for a representative sequence. 
![fastqc report](/imagesreadme/report_fastqc.png)

The quality of the first and last bases of the sequences is not satisfying. This is due to the artefacts of sequencing, that also cause the presence of adapters at the end of sequences. We need to trim the sequences before starting the assembly. 

**2) Sequences trimming:** 

Trimmomatic is a tool to perform quality trimming, adapter clipping, and eliminate reads that are too short. We use the script **trimmomatic.sh** to run Trimmomatic for each pair of sequences, with the following parameters: ILUMINACLIP:2:30:10, HEADCROP: 9, MINLEN: 100.  
Then we verify the quality of the trimmed sequenced using fastqc on the outputs from Trimmomatic:  **control_fastqc_after_trimm.sh**.
![fastqc after trimming report](/imagesreadme/report_fastqc_aftertrimm.png)

Overall, the quality is now satisfying, and we can proceed to the assembly, using the fastq sequences stored in the outputs_trimmomatic directory. 

**3) _De novo_ assembly:** 
Once the sequences have been trimmed and are of a good quality, we use Trinity to assemble the transcriptome. We do not separate the samples for the assembly, and use the paired-end sequences of the 6 samples. We only separate the left and right reads of the samples, because Trinity processes all the right reads together (FASTQ1), and all the left reads together (FASTQ2). We specify the strand-specific RNA-seq read orientation (--SS_lib_type): in our case, the type of sequencing was 'fr-firststrand', which for Trinity corresponds to "RF".
For the script, see **assembly_trinity.sh**. 

We now have an assembled transcriptome (Trinity_RF.fasta). The next step is to annotate it. 

### Step 2: Data annotation
???

### Step 3: Transcript expression quantification

We now wish to quantify the expression for each sample. To do so, we use salmon. Salmon comprises tools to perform transcript-level quantification from RNA-seq reads using selective alignment. We use salmon index and salmon quant. 

**1) Creating an index with salmon:** 

To perform quantification, salmon first requires the creation of a specific index for the assembled transcriptome (Trinity_RF.fasta). This is the first part of the **salmon_quantif.sh**.

**2) Quantification:** 

Then, we proceed to actual quantification with salmon quant. We run salmon quant on each pair of sequences (second part of the **salmon_quantif.sh** script). 

### Step 4: Differential expression analyses
???
