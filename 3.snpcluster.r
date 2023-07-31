library(GenomicRanges)
library(ggplot2)
library(Vennerable)


readSeqinfo = function(file, source) {
	df = read.delim(file, sep='\t', head=F, as.is=T)
	Seqinfo(seqnames=df[,1], seqlengths=df[,2], genome=source)
}

seqinfo = readSeqinfo('~/genome/osa_TIGR7/TIGR7_genome.chrsize', 'TIGR7')

txdb = read.table('~/genome/osa_TIGR7/TIGR7_genes.gff', stringsAsFactors=F)
genedb = txdb[txdb$V3 == 'gene',]
genegr = GRanges(seqnames=genedb$V1, strand=genedb$V7, seqinfo=seqinfo,
			ranges=IRanges(start=genedb$V4, end=genedb$V5))
names(genegr) = gsub('ID=(.*);Name.+', '\\1', genedb$V9)

tedb = read.table('~/genome/osa_TIGR7/TIGR7_TE.gff', stringsAsFactors=F)
tename = gsub(';.*', '', gsub('.+Name=(.*)$', '\\1', tedb$V9))
tename = gsub('%29', ')', gsub('%28', '(', tename))
#ir = tedb$V1 %in% c('ChrSy', 'ChrUn')
tegr = GRanges(seqnames=tedb$V1, strand='*', seqinfo=seqinfo,
			ranges=IRanges(start=tedb$V4, end=tedb$V5))
names(tegr) = tename
#tegr = tegr[!ir]


####################### FUNCTIONs ###########################
readSNPbias2GRange = function(snpfile, seqinfo) {
	snp = read.table(snpfile)
	colnames(snp) = c('Chrom', 'Pos', 'SNP', 'm9311xNip', 'mNipx9311', 'Bias')
	grg = GRanges()
	
	grg = GRanges(seqnames=snp$Chrom, strand='*', snp=snp$SNP, bias=snp$Bias, 
		seqinfo=seqinfo, ranges=IRanges(start=snp$Pos, width=1))
		
	grg
}


################################################################################
grsnps = sapply(c('H3K27me3', 'H3K9me2'), function(s) {
	readSNPbias2GRange(paste0(s, '.peaks.mpbias.snp'), seqinfo)
}, simplify=F)

grgsnp = sapply(grsnps, function(gr) {
	hits = findOverlaps(gr, genegr)
	grg = genegr[subjectHits(hits)]
	grg$snp = gr[queryHits(hits)]$snp
	grg$bias = gr[queryHits(hits)]$bias
	grg
})

gsnpinfo = sapply(grgsnp, function(grg) {
	tapply(grg$bias, INDEX=names(grg), FUN=function(x) {
		type = 'M&P'
		if(all(x=='M'))
			type = 'M'
		else if(all(x=='P'))
			type = 'P'
		type
	})
})

for (s in c('H3K27me3', 'H3K9me2')) {
	grg = genegr[names(gsnpinfo[[s]])]
	grg$Bias = gsnpinfo[[s]]
	write.table(as.data.frame(grg), file=paste0(s,'.mpbias.gene.txt'), sep='\t')
}