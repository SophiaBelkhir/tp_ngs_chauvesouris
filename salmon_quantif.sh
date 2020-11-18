#!/bin/bash
#Script to create index with salmon for quantification
#Create working directory
data="/home/rstudio/data/mydatalocal/data_from_ftp"
mkdir -p $data
cd $data

#Create directory for the current task (salmon index and quantif)
output_index="outputs_salmon"
mkdir -p $output_index
cd $output_index

#Run salmon index
salmon index -t $data/trinity_data/Trinity_RF.fasta -i $data/$output_index -p 4

#For the quantification 

#Library of each pair of reads, to use in the loop for salmon quant
library_reads="Lib1_31_20_S1
Lib2_31_20_S2
Lib3_31_20_S3
Lib4_31_20_S4
Lib5_31_20_S5
Lib6_31_20_S6"

#run salmon quant for each sample
for lib in $library_reads
do
  salmon quant -i $data/$output_index -l A -1 $data/outputs_trimmomatic/${lib}_1_paired.fq.gz -2 $data/outputs_trimmomatic/${lib}_2_paired.fq.gz --validateMappings -o $data/$output_index
done
