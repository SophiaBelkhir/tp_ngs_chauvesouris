#!/bin/bash
#Script to trim the sequences with Trimmomatic

#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (trimming with Trimmomatic)
data_trimmomatic="trimmomatic_data"
mkdir -p $data_trimmomatic
cd $data_trimmomatic

#Get the data (the sequences to trim)
home_fastq="/home/rstudio/data/mydatalocal/data_from_ftp/igfl_data_seq/Projet_31_20_UE_NGS_2020/FASTQ"

#Create the outputs files for trimmomatic
#loop for to create the paired and unpaired output files for each pair of sequence
library_seq="Lib1_31_20_S1
Lib2_31_20_S2
Lib3_31_20_S3
Lib4_31_20_S4
Lib5_31_20_S5
Lib6_31_20_S6"

#now run trimmomatic on each sequence using a loop
for seq in $library_seq
do
  java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 -phred33 $home_fastq/${seq}_R1_001.fastq.gz $home_fastq/${seq}_R2_001.fastq.gz $data/$data_trimmomatic/${seq}_R_unpaired.fq.gz $data/$data_trimmomatic/${seq}_F_paired.fq.gz $data/$data_trimmomatic/${seq}_R_unpaired.fq.gz $data/$data_trimmomatic/${seq}_F_unpaired.fq.gz ILLUMINACLIP:$data/adapt.fasta:2:30:10 HEADCROP:9 MINLEN:100
done
