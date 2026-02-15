#!/bin/sh -x

bam_folder="/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/bam_files_SARS-CoV-2"

NAME1=$1  	#SRR
NAME2=$2	#sample name

echo $NAME2	#name

cd $bam_folder
#align with hisat2 and create sorted indexed bam files
#SARS-CoV-2 index!!!!!!!!!!!!!!!
#--mm - Use memory-mapped I/O to load the index, rather than typical file I/O. 
#Memory-mapping allows many concurrent bowtie processes on the same computer to share the same memory image of the index (i.e. you pay the memory overhead just once).

#-F 4
#Remove reads with FLAG 4 = unmapped
#hisat2 \
#  -x /gpfs0/tals/projects/data/Genomes/SARS-CoV-2_Delta/SARS-CoV_index \
#  --mm \
#  -1 ${NAME1}_R1_001.fastq.gz \
#  -2 ${NAME1}_R2_001.fastq.gz \
#| samtools view -b -F 4 \
#| samtools sort -o $bam_folder/$NAME2.virus.sorted.bam

#samtools index "$bam_folder/$NAME2.virus.sorted.bam"

samtools view "$bam_folder/$NAME2.virus.sorted.bam" | wc -l >> "$bam_folder/viral_RNA_count.txt"
echo $NAME2".virus.sorted.bam" >> "$bam_folder/viral_RNA_count.txt"


