#! /bin/bash

for fq in reads/*.fq.gz; do
    bam=${fq/fq.gz/raw.bam}
	bam=${bam/reads/mapped}

    # mapping fastq reads to the reference genome
	bowtie2 -p 10 --local ~/genome/osa_MSU7/idx_bowtie2/MSU7 -U $fq | samtools view -bS - > $bam

    # remove PCR duplicates
    samtools rmdup $bam - | samtools sort -@ 10 -m 1G - -o ${bam/raw.bam/bam}

    # index bam file
    samtools index ${bam/raw.bam/bam}
done;

