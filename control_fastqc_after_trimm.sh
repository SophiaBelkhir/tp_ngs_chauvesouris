#!/bin/bash
#Script to do quality control of the outputs of Trimmomatic (after trimming)
 
#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (quality control with fastqc after trimmommatic)
trimm_fastqc="fastqc_after_trimmomatic"
mkdir -p $trimm_fastqc
cd $trimm_fastqc
 
#Get the data (the sequences to control)
home_data="/home/rstudio/data/mydatalocal/data_from_ftp/outputs_trimmomatic"
fastq=$home_data/"*_paired.fq.gz"

#loop for to apply fastqc to each sequence
for sample in fastq
do
  echo $sample
  fastqc $sample --outdir $data/$trimm_fastqc
done 