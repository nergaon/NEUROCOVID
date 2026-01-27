#!/bin/sh -x

bam_folder="/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/bam_files"
leafcutter_folder="/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/leafcutter_0.2.9"
leafcutter_scripts="/gpfs0/tals/projects/software/leafcutter_0.2.9"	

NAME1=$1  	#SRR
NAME2=$2	#sample name

echo $NAME2	#name

cd $bam_folder
#align with hisat2 and create sorted indexed bam files
#human index!!!!!!!!!!!!!!!
#--mm - Use memory-mapped I/O to load the index, rather than typical file I/O. 
#Memory-mapping allows many concurrent bowtie processes on the same computer to share the same memory image of the index (i.e. you pay the memory overhead just once).
hisat2 -x /gpfs0/tals/projects/data/Genomes/hg38/grch38/genome --add-chrname --mm -1 $NAME1"_R1_001.fastq.gz" -2 $NAME1"_R2_001.fastq.gz" -S "$bam_folder/$NAME2.sam"
chmod 770 "$bam_folder/$NAME2.sam"
samtools view -b "$bam_folder/$NAME2.sam" -o "$bam_folder/$NAME2.bam"
chmod 770 "$bam_folder/$NAME2.bam"
rm -f "$bam_folder/$NAME2.sam"
samtools sort "$bam_folder/$NAME2.bam" -o "$bam_folder/$NAME2.sort.bam"
chmod 770 "$bam_folder/$NAME2.sort.bam"
rm -f "$bam_folder/$NAME2.bam"
samtools index "$bam_folder/$NAME2.sort.bam"
chmod 550 "$bam_folder/$NAME2.sort.bam"
chmod 550 "$bam_folder/$NAME2.sort.bam.bai"
samtools flagstat "$bam_folder/$NAME2.sort.bam" >> "$bam_folder/flagStat.txt"
echo $NAME2".sort.bam" >> "$bam_folder/flagStat.txt"

#convert bam files to junc files
#chmod 770 *
bamfile=$NAME2".sort.bam"
juncfile=$NAME2.junc
echo Converting $bamfile to $juncfile
$leafcutter_scripts/bin/regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $bamfile -o $leafcutter_folder/$bamfile.junc
chmod 770 *

# Program:        regtools
# Version:        0.5.2
# Usage:          regtools junctions extract [options] indexed_alignments.bam
# Options:
                # -a INT  Minimum anchor length. Junctions which satisfy a minimum
                         # anchor length on both sides are reported. [8]
                # -m INT  Minimum intron length. [70]
                # -M INT  Maximum intron length. [500000]
                # -o FILE The file to write output to. [STDOUT]
                # -r STR  The region to identify junctions
                         # in "chr:start-end" format. Entire BAM by default.
                # -s INT  Strand specificity of RNA library preparation
                         # (0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). REQUIRED
                # -t STR  Tag used in bam to label strand. [XS]