#!/bin/bash
#Script to build database and run BLAST on it
#to run on the Transdecoder output (that detected the ORFs)

#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (blast)
output_blast="blast_data"
mkdir -p $output_blast
cd $output_blast

#get necessary files
HumanCDS=Homo_sapiens.GRCh38.cds.fa  #for the database, human cds
query=$data/transdec_data/Trinity_RF.fasta.transdecoder.cds  #the query, our output from transdecoder

#build reference database with makeblastdb
/softwares/ncbi-blast-2.10.1+/bin/makeblastdb -in $HumanCDS -dbtype 'nucl' -out 'DatabaseHumanCDS' -parse_seqids 

#blast our query (Trinity_RF.fasta.transdecoder.cds) against the reference database we juist built
/softwares/ncbi-blast-2.10.1+/bin/blastn -db "DatabaseHumanCDS" -query $query -evalue 1e-4 -outfmt 6 -out "DataBats.blast" -max_target_seqs 1
