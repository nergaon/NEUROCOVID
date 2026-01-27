#!/bin/sh
#COMMAND: qsub -cwd -V -q tals.q execute_all_HN2.sh /gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/bam_files/srr_names.txt
#parameters_file_name - text file with input names. Begin at first line (no headers, etc.)

bam_folder="/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/bam_files/"
#remove the list and programs from the previous run
rm -f $bam_folder/*bam*
rm -f *sh.*
rm -f $bam_folder/flagStat.txt
while IFS=$'\t' read -r line || [[ -n "$line" ]]; do
	#r=$(qstat -u nergaon|wc -l)
	#while [ $r -ge 390 ]; #if i have more than 390 jobs the program will run only 390 and wait. ge-greater or equal
	#do
	#	r     =$(qstat -u nergaon|wc -l)
	#	sleep 30
	#done
	arrIN=(${line//$'\t'/ })
	arg1=${arrIN[0]/$'\r'/}			#SRR
	arg2=${arrIN[1]/$'\r'/}			#name 
	echo $arg1
	qsub -cwd -V -q tals.q execute_single_bam_junc_HN1.sh $arg1 $arg2 #create bam and junc files
done < "$1"
