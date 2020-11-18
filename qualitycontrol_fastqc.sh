#!/bin/bash
#Script to do quality control of sequences before assembly
 
#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (quality control with fastqc)
data_fastqc="outputs_fastqc"
mkdir -p $data_fastqc
cd $data_fastqc
 
#Get the data (the sequences to control)
home_fastq="/home/rstudio/data/mydatalocal/data_from_ftp/igfl_data_seq/Projet_31_20_UE_NGS_2020/FASTQ"
fastq=$home_fastq/"*.gz"

#loop for to apply fastqc to each sequence
for sample in $fastq
do
  echo $sample
  fastqc $sample --outdir $data/$data_fastqc
done 