#! /bin/bash

# mpileup to generate bcf file
samtools mpileup -uIDV -f ~/genome/osa_MSU7/MSU7_genome.fa H3K27me3_93xN.bam H3K27me3_Nx93.bam | bcftools view -bvcg - > H3K27me3.bcf 
samtools mpileup -uIDV -f ~/genome/osa_MSU7/MSU7_genome.fa H3K9me2_93xN.bam H3K9me2_Nx93.bam | bcftools view -bvcg - > H3K9me2.bcf 


# transforming and filtering bcf file
for bcf in *.bcf; do 
	snp=${bcf/bcf/snp}
    mbias=${bcf/bcf/mbias}
    mbiasbed=${bcf/bcf/mbias.bed}
    snpbiasgene=${bcf/bcf/snpbias.gene.anno}
    snpbiasexon=${bcf/bcf/snpbias.exon.anno}

    # transform bcf file to snp file
	bcftools view $bcf | awk -F'\t' '
	BEGIN{OFS="\t"} /^#CHROM/{print $1,$2,$4 "|" $5,$10,$11}
	!/^#/{
		split($10,X,":"); split($11,Y,":");
		print $1, $2, $4 "|" $5, (X[3]-X[4]) "|" X[4], (Y[3]-Y[4]) "|" Y[4]
	}
	' > $snp

    # calulation of maternal bias
	cat $snp | awk '
    BEGIN{OFS="\t"}
	/^/{print}
	!/^/{print $1, $2, $3, $4/(2-$4)*2-1, (1-$5)/(1+$5)*2-1}
	' > $mbias

    # filtering maternal bias of significant SNPs
	cat $mbias | awk '
	BEGIN{OFS="\t"}
	$4 > .5 && $5 > .5{
		print $1, $2, $2+1, $1 ":" $2 "." $3, int($4*1000)+$5/10
	}
	$4 < -.5 && $5 < -.5{
		print $1, $2, $2+1, $1 ":" $2 "." $3, int($4*1000)+$5/10
	}' > $mbiasbed

    # annotation of maternal bias in gene and exon
	bedtools closest -D b -a $mbiasbed -b $HOME/genome/osa_MSU7/MSU7_gene.gene.bed > $snpbiasgene
	bedtools closest -D b -a $mbiasbed -b $HOME/genome/osa_MSU7/MSU7_gene.exon.bed > $snpbiasexon
done