#! /bin/bash
#Script to run Transdecoder to detect the ORFs 

#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (ORF detection with TransDecoder)
output_transdec="transdec_data"
mkdir -p $output_transdec

#launch Transdecoder

#Step 1: extract the long open reading frames
TransDecoder.LongOrfs -t $data/trinity_data/Trinity_RF.fasta --gene_trans_map $data/trinity_data/Trinity_RF.fasta.gene_trans_map -S -O $output_transdec

#Step 2: predict the coding regions
TransDecoder.Predict -t $data/trinity_data/Trinity_RF.fasta --single_best_only -O $output_transdec