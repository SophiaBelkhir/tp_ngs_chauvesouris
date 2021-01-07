#!/bin/bash
#Script to do de novo assembly with Trinity
#to run on sequences which quality has been controlled with fastqc and cleaned with trimmomatic

#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (de novo assembly with Trinity)
output_trinity="trinity_data"
mkdir -p $output_trinity
cd $output_trinity

#Get the data (the sequences to assemble)
home_fastq=$data/"outputs_trimmomatic"
FASTQ1=$(ls $home_fastq/*1_paired.fq.gz |paste -d "," -s) #all the left reads
FASTQ2=$(ls $home_fastq/*2_paired.fq.gz |paste -d "," -s) #all the right reads

#Run Trinity 
Trinity --seqType fq --max_memory 14G --left $FASTQ1 --right $FASTQ2 --CPU 4 --SS_lib_type RF --output $data/$output_trinity
