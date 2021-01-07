#!/bin/bash
#Script to dowload human CDS from ftp ensembl

#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory to store the file 
output_blast="blast_data"
mkdir -p $output_blast
touch $output_blast/Homo_sapiens.GRCh38.cds.fa.gz

wget -O $output_blast/Homo_sapiens.GRCh38.cds.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz

gunzip $output_blast/Homo_sapiens.GRCh38.cds.fa.gz