# Practicals NGS Bats 
## Objectives - Characterization of the IFN response in bats
The objective is to investigate the specificities of the bats immune system, in order to understand how bat's immunity has adapted to long-term association with viruses, and why bats harbor a high number of viruses but appear mostly asymptomatic. 
To do so, we want to characterize the interferon (IFN) response in bats. We thus aim at studying the transcriptomic changes in response to IFN stimulation in bat cells, and in particular, we want to identify which interferon-stimulated genes (ISGs) are differentially expressed after IFN stimulation.

We focus on *Myotis velifer*, a bat species for which no transcriptomic data was available yet. We will thus need to assemble the transcriptome of *Myotis velifer*, and because no reference transcriptome is available for this species, we perform a *de novo* transcriptome assembly. Then we will have to annotate this transcriptome. Finally, we will quantify the expression of the transcripts, and perfomr a differential expression analysis.

### Step 0: Data obtention
RNA-seq was extracted from bat cells exposed or not to IFN. 
Libraries were prepared and paired-end RNA-seq was performed using the Illumina NextSeq500 platform. For each condition, 3 samples were sequenced, so we have a total of 6 files containing reads. The sequences are under fastq format. 

To download the raw RNA-seq data from the IGFL sequencing platform, use the command line in the script **importdataset.sh**.

### Step 1: Quality control
The first step is to verify the quality of the reads, and eventually trim them before proceeding to the assembly.

**1) Quality control with fastqc:** 

We first control the quality of the reads. We run fastqc on each of the obtained sequences, using the script **qualitycontrol_fastqc.sh**. Fastqc uses several statistics to analyze the quality of sequences, and provides html files with a summary for each file. Here is an overview of the fastqc report obtained for a representative sample.
![fastqc report](/imagesreadme/report_fastqc.png)

The quality of the first and last bases of the sequences is not satisfying. This is due to the artefacts of sequencing, for which the first bases are not accurate, and to the presence of adapters at the end of sequences. We need to trim the sequences before starting the assembly.

**2) Sequences trimming with trimmomatic:** 

Trimmomatic is a tool to perform quality trimming, adapter clipping, and eliminate reads that are too short. We use the script **trimmomatic.sh** to run Trimmomatic for each pair of sequences, with the following parameters: ILUMINACLIP:2:30:10 (to use in combination with the file adapt.fasta which contains the sequences of the Illumina adapters), HEADCROP: 9 (to crop the first 9 bases of each read), MINLEN: 100 (to only keep the reads with a minimum length of 100 nucleotides).  
Then we verify the quality of the trimmed sequenced using fastqc on the outputs from Trimmomatic:  **control_fastqc_after_trimm.sh**.
![fastqc after trimming report](/imagesreadme/report_fastqc_aftertrimm.png)

Overall, the quality is now satisfying, and we can proceed to the assembly, using the fastq sequences stored in the outputs_trimmomatic directory.

### Step 2: De novo assembly 
Once the sequences have been trimmed and are of a good quality, we use Trinity to assemble the transcriptome. We do not separate the samples by condition for the assembly, and use the paired-end sequences of the 6 samples altogether. We only separate the left and right reads of the samples, because Trinity processes all the right reads together (stored in the variable FASTQ1), and all the left reads together (stored in the variable FASTQ2). We specify the strand-specific RNA-seq read orientation (--SS_lib_type): in our case, the type of sequencing was 'fr-firststrand', which for Trinity corresponds to "RF".
For the script, see **assembly_trinity.sh**. 

We now have an assembled transcriptome (Trinity_RF.fasta), and a table making the correspondence between the gene identifier and the transcript isoforms (Trinity_RF.fasta.gene_trans_map) . The next step is to annotate it. 

### Step 3: Data annotation
In order to annotate the transcriptome, we need: 1) to identify the ORFs from the assembled transcriptome 2) to predict the likely coding regions 3)  to select and create the database of known CDS against which we will use BLAST to find homologies 4) to perform a local alignment with BLAST to identify genes 

**1) Identifying the ORFs with TransDecoder.Longorfs:** 
We run TransDecoder in two steps, first TransDecoder.Longorfs identifies ORFs within transcript sequences, using the command line in the script **transdec.sh**. The input file is the assembled transcriptome obtained with Trinity, and we use the fasta.gene_trans_map for gene transmap. We also specify -s to indicate the strand specificity.

**2) Identifying the CDS with TransDecoder.Predict:** 
We then use TransDecoder.Predict to identify candidate coding regions within transcript sequences, with the corresponding command line in the script **transdec.sh**. We specify –single best only to only keep the best candidate for each transcript. 

We obtain trinity.fasta.transdecoder.cds, a file with all the CDS predicted from the assembled transcriptome. In order to know to what known coding sequences these correspond, we look for homologies with Basic Local Alignment Search Tool. 

**3) Make a blast dadabase with makeblastdb:**
To annotate the CDS from *Myotis velifer* transcriptome, we compare them to human CDS. We thus need to build a blast database containing human CDS.
We download the known human cds using the wget command from the script **get_human_cds.sh**.

Then we use makeblastdb, using the file with human CDS as input file, and specifying -dbtype 'nucl' -parse_seqids, because we deal with nucleotide sequences, under a fasta format (for command line, see the first part of the **blast.sh** script). 

**4) Perform alignment with blastn:**
We then perform the alignment using blastn (because both our query and our database are nucleotides), see the scrip **blast.sh**. 
We use the output from TransDecoder as the query, and the output from makeblastdb as the database for reference. We define the manimal e-value as 1e-4, and max_target_seqs as 1 to only keep one hit by contig. We choose -outfmt 6 to visualize the output as a table summarizing the main features of each hit. 

### Step 4: Transcript expression quantification
We now wish to quantify the expression for each sample. We use salmon to align the reads on the transcript and quantify how many reads align on the transcript. Salmon comprises several tools to perform transcript-level quantification from RNA-seq reads using selective alignment, and we use salmon index and salmon quant. 

**1) Creating an index with salmon index:**

To perform quantification, salmon first requires the creation of a specific index for the assembled transcriptome (Trinity_RF.fasta). This is the first part of the **salmon_quantif.sh**.

**2) Quantification with salmon quant:**

Then, we proceed to actual quantification with salmon quant. We run salmon quant on each pair of sequences (second part of the **salmon_quantif.sh** script). We use --validateMappings and --gcBias options. 

Of note, only 40% of the reads aligned to the transcriptome, when we expected more around 90%. This is probably because the left and right reads sometimes overlapped, so salmon running in "paired-end" mode does not handle it well. We could try running it as "single-end".

We obtained the count tables quant.sf for each sample. We have counts for each isoform, but now want the counts at the level of the gene, to be able to do differential expression analysis. 

### Step 5: Differential expression analyses

**1) Building a count table:**
