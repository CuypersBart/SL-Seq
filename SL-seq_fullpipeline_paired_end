#This script will automatically find paired end reads that are named according to illumina conventional and generate gene count tables for a specified reference genome and annotation file
#The script will first align the reads to the reference genome of the host cell, after which the unmapped reads will be mapped to the parasite reference genome and a count table will be generated
#pipeline tested with TopHat 2.1.0, Java 1.8.0_45, Pysam 0.9.1.4, HTSeq/0.6.1, Python-2.7.9, samtools-1.3.1, Trimmomatic 035, bwa 0.7.15
#please make sure these tools are installed and available in $PATH
#index host cell reference genome with 'bowtie2-build' and give the index the basename 'hostcell_bowtie'
#index tryp reference genome by running 'bwa index <tryprefgenome.fasta>'
#adapt following parameters accordingly
gff_file="SL_tryp.gff3"
refgenome_human_bowtie="tryprefgenome.fasta"

arr=$(ls *R1_001.fastq.gz)
for i in $arr
do
	fastqc $i
	i2=${i/'R1_001.fastq.gz'/'R2_001.fastq.gz'}
	fastqc $i2
	java -jar trimmomatic-035.jar SE \
	-phred33 \
	$i2 \
	${i2/'_R2'/'_cropped_trimmed_R2'}  \
	HEADCROP:21
	fastqc ${i2/'_R2'/'_cropped_trimmed_R2'}
	java -jar trimmomatic-035.jar PE \
	-threads 5 -phred33 \
	$i \
	${i2/'_R2'/'_cropped_trimmed_R2'} \
	${i/'_R1'/'_trimmed_R1'} \
	${i/'_R1'/'_unpaired_R1'} \
	${i2/'_R2'/'_trimmed_R2'} \
	${i2/'_R2'/'_unpaired_R1'} \
	LEADING:20 TRAILING:20 ILLUMINACLIP:adapter.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36
	fastqc ${i/'_R1'/'_trimmed_R1'}
	fastqc ${i2/'_R2'/'_trimmed_R2'} 
	tophat2 -p 20 -o ${i2/'_R2'/''} hostcell_bowtie ${i/'_R1'/'_trimmed_R1'} ${i2/'_R2'/'_trimmed_R2'}	
done

arr=$(find . -maxdepth 4 -iregex .*unmapped.bam)
for i in $arr
do
	i2=${i/'bam'/'1.fq'}
	i3=${i/'bam'/'2.fq'}
	i3bis=${i/'bam'/'unpaired.fq'}
	i4=${i/'bam'/'tryp.bam'}
	i5=${i/'bam'/'tryp.sorted.bam'}
	samtools sort -n -o ${i/'.bam'/'_L.namesorted.bam'} $i
	samtools fastq -1 ${i/'.bam'/'_L_R1.fastq'} -2 ${i/'.bam'/'_L_R2.fastq'} -0 ${i/'.bam'/'_L_neither.fastq'} -s ${i/'.bam'/'_L_singleton.fastq'}  ${i/'.bam'/'_L.namesorted.bam'} 
	bwa mem -t 20 $refgenome_human_bowtie ${i/'.bam'/'_L_R1.fastq'} ${i/'.bam'/'_L_R2.fastq'} | samtools view -bS > ${i/'.bam'/'_Trypaligned.bam'}
	samtools sort -o ${i/'.bam'/'_Trypaligned.sorted.bam'} ${i/'.bam'/'_Trypaligned.bam'}
	samtools index ${i/'.bam'/'_Trypaligned.sorted.bam'}
	samtools flagstat ${i/'.bam'/'_Trypaligned.sorted.bam'} > ${i/'.bam'/'_Trypaligned.sorted.txt'}
	htseq-count -r name -f bam -t gene -i ID -s no ${i/'.bam'/'_Trypaligned.sorted.bam'} $gff_file > ${i/'.bam'/'_Trypaligned.counttable'}
done
