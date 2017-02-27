#!/bin/bash
#This script will automatically find single-end reads that are named according to Illumina convention and generate gene count tables for a specified reference genome and annotation file
#The script will first align the reads to the reference genome of the host cell, after which the unmapped reads will be mapped to the parasite reference genome and a count table will be generated
#pipeline tested with Java 1.8.0_45, Pysam 0.9.1.4, HTSeq/0.6.1, Python-2.7.9, samtools-1.3.1, Trimmomatic 035, bwa 0.7.15
#please make sure these tools are installed and available in $PATH
#index host cell reference genome with 'bowtie2-build' and give the index the basename 'hostcell_bowtie'
#index tryp reference genome by running 'bwa index <tryprefgenome.fasta>'
#adapt following parameters accordingly
gff_file="SL_tryp.gff3"
refgenome_human_bowtie="tryprefgenome.fasta"

arr=$(ls *R1_001.fastq)
for i in $arr
do
	fastqc $i
	java -jar trimmomatic-035.jar SE \
	-phred33 -trimlog r.log \
	$i \
	${i/'_R1'/'_trimmed_R1'}  \
	LEADING:20 TRAILING:20 ILLUMINACLIP:adapter.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36
	fastqc ${i/'_R1'/'_trimmed_R1'} 
	tophat2 -p 20 -o ${i/'_R1_001.fastq'/''} hostcell_bowtie ${i/'_R1'/'_trimmed_R1'}
done

arr=$(find . -maxdepth 4 -iregex .*unmapped.bam)
for i in $arr
do
	samtools sort -n -o ${i/'.bam'/'_L.namesorted.bam'} $i
	samtools fastq ${i/'.bam'/'_L.namesorted.bam'} > ${i/'.bam'/'_L_R1.fastq'}
	bwa mem -t 20 tryprefgenome.fasta ${i/'.bam'/'_L_R1.fastq'} | samtools view -F 2048 -bS > ${i/'.bam'/'_Trypaligned.bam'}
	samtools sort -o ${i/'.bam'/'_Trypaligned.sorted.bam'} ${i/'.bam'/'_Trypaligned.bam'}
	samtools index ${i/'.bam'/'_Trypaligned.sorted.bam'}
	samtools flagstat ${i/'.bam'/'_Trypaligned.sorted.bam'} > ${i/'.bam'/'_Trypaligned.sorted.txt'}
	htseq-count -r name -f bam -t gene -i ID -s no ${i/'.bam'/'_Trypaligned.sorted.bam'} SL_tryp.gff3 > ${i/'.bam'/'_Trypaligned.counttable'}
done
