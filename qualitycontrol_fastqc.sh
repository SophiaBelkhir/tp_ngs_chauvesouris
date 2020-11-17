 #!/bin/bash
 #Quality control of sequences before assembly
 #Go to the right working directory
 cd /home/rstudio/data/mydatalocal/data_from_ftp/
 #create a new folder to put the output of fastq
 mkdir outputs_fastqc

#loop for to apply fastqc to each sequence
for sample in /home/rstudio/data/mydatalocal/data_from_ftp/igfl_data_seq/Projet_31_20_UE_NGS_2020/FASTQ/*.gz
do
  echo $sample
  fastqc $sample --outdir /home/rstudio/data/mydatalocal/data_from_ftp/outputs_fastqc/
done 