options(stringsAsFactors=F)
library(WriteXLS)
rm(list=ls())
setwd("/sc/orga/projects/STARNET/ariella/Clint/")

files = list.files(".","annotation.xls")
# files = list.files(".","TCF21.xls")
# files = list.files(".","SMAD.xls")

## merging the files together and adding some more annotations -- 
cis_snps = "/sc/orga/projects/STARNET/ariella/Clint/snpsInterest"
trans_snps = "/sc/orga/projects/STARNET/ariella/Clint/transInterestedIn"
snp_annotation = read.delim("/sc/orga/projects/STARNET/ariella/Clint/SNP_annotations")

load("/sc/orga/projects/STARNET/ariella/cis_trans_causlity/rerun/geneNames_all.Rdata")

# snpFileName = "/sc/orga/scratch/cohaia01/MGH/SnpsInterestedIn"
# geneAnnotationFN = "/sc/orga/scratch/cohaia01/MGH/Expression/human gene annotations_v2.txt"

cis = read.delim(cis_snps)
# gene_info = read.delim(geneAnnotationFN)
cis$marker =paste(cis$chr,cis$pos,sep=":")

trans = read.delim(trans_snps)
# gene_info = read.delim(geneAnnotationFN)
trans$marker =paste(trans$chr,trans$pos,sep=":")

snps = rbind(cis,trans)

output = lapply(files,function(file){ #for(file in files){
	x = read.delim(paste(strsplit(file,".xls")[[1]][1],".txt",sep=""))
	gene_anno = as.matrix(sapply(1:nrow(x),function(X){geneNames[which(geneNames$id == x$gene[X])[1],'name']}))
	snp_anno = t(sapply(1:nrow(x),function(X){snps[which(snps$marker == x$snps[X])[1],]}))
	colnames(gene_anno) = "geneName"
	colnames(snp_anno) = colnames(snps)
	annotated = cbind(x,gene_anno,snp_anno)
	return(annotated)
	
})

names(output) = files
AOR_cis = output[[1]]
AOR_trans = output[[2]]

MAM_cis = output[[3]]
MAM_trans = output[[4]]

AOR_cis = AOR_cis[,-c('Gene','X..5L..','gwas_cis_pvalue','gwas_snp')]
MAM_cis = MAM_cis[,-c('Gene','X..5L..','gwas_cis_pvalue','gwas_snp')]

# AOR_cis = output[[1]]
# AOR_cis = AOR_cis[,-c('Gene','X..5L..','gwas_cis_pvalue','gwas_snp')]
# MAM_cis = output[[2]]
# MAM_cis = MAM_cis[,-c('Gene','X..5L..','gwas_cis_pvalue','gwas_snp')]

AOR_cis$a1 =sapply(AOR_cis$snps,function(x){snp_annotation$a1[snp_annotation$pos == x]})
AOR_cis$a2 =sapply(AOR_cis$snps,function(x){snp_annotation$a2[snp_annotation$pos == x]})
AOR_cis$exp_freq_a1 = sapply(AOR_cis$snps,function(x){snp_annotation$exp_freq_a1[snp_annotation$pos == x]})
AOR_cis$MSSM_MAF =  sapply(AOR_cis$snps,function(x){snp_annotation$MAF[snp_annotation$pos == x]})

MAM_cis$a1 =sapply(MAM_cis$snps,function(x){snp_annotation$a1[snp_annotation$pos == x]})
MAM_cis$a2 =sapply(MAM_cis$snps,function(x){snp_annotation$a2[snp_annotation$pos == x]})
MAM_cis$exp_freq_a1 = sapply(MAM_cis$snps,function(x){snp_annotation$exp_freq_a1[snp_annotation$pos == x]})
MAM_cis$MSSM_MAF =  sapply(MAM_cis$snps,function(x){snp_annotation$MAF[snp_annotation$pos == x]})

AOR_trans$a1 =sapply(AOR_trans$snps,function(x){snp_annotation$a1[snp_annotation$pos == x]})
AOR_trans$a2 =sapply(AOR_trans$snps,function(x){snp_annotation$a2[snp_annotation$pos == x]})
AOR_trans$exp_freq_a1 = sapply(AOR_trans$snps,function(x){snp_annotation$exp_freq_a1[snp_annotation$pos == x]})
AOR_trans$MSSM_MAF =  sapply(AOR_trans$snps,function(x){snp_annotation$MAF[snp_annotation$pos == x]})

MAM_trans$a1 =sapply(MAM_trans$snps,function(x){snp_annotation$a1[snp_annotation$pos == x]})
MAM_trans$a2 =sapply(MAM_trans$snps,function(x){snp_annotation$a2[snp_annotation$pos == x]})
MAM_trans$exp_freq_a1 = sapply(MAM_trans$snps,function(x){snp_annotation$exp_freq_a1[snp_annotation$pos == x]})
MAM_trans$MSSM_MAF =  sapply(MAM_trans$snps,function(x){snp_annotation$MAF[snp_annotation$pos == x]})


colnames(AOR_cis)[colnames(AOR_cis)=="cond_pval_MaxGvnGwas"] = "cond_pval_MaxGvnSNP"
colnames(AOR_cis)[colnames(AOR_cis)=="cond_rsqd_MaxGvnGwas"] = "cond_rsqd_MaxGvnSNP"
colnames(AOR_cis)[colnames(AOR_cis)=="cond_pval"] = "cond_pval_SNPgvnMax"
colnames(AOR_cis)[colnames(AOR_cis)=="cond_rsqd"] = "cond_rsqd_SNPgvnMax"

colnames(MAM_cis)[colnames(MAM_cis)=="cond_pval_MaxGvnGwas"] = "cond_pval_MaxGvnSNP"
colnames(MAM_cis)[colnames(MAM_cis)=="cond_rsqd_MaxGvnGwas"] = "cond_rsqd_MaxGvnSNP"
colnames(MAM_cis)[colnames(MAM_cis)=="cond_pval"] = "cond_pval_SNPgvnMax"
colnames(MAM_cis)[colnames(MAM_cis)=="cond_rsqd"] = "cond_rsqd_SNPgvnMax"

WriteXLS(c("AOR_cis","MAM_cis","AOR_trans","MAM_trans"),ExcelFileName = "/sc/orga/projects/STARNET/ariella/Clint/annotatedSNPs_updated.xls",SheetNames = c("AOR_cis","MAM_cis","AOR_trans","MAM_trans"),row.names =F)

# WriteXLS(c("AOR_cis","MAM_cis"),ExcelFileName = "/sc/orga/projects/STARNET/ariella/Clint/TCF21_annotatedSNPs.xls",SheetNames = c("AOR_cis","MAM_cis"),row.names =F)

 # WriteXLS(c("AOR_cis","MAM_cis"),ExcelFileName = "/sc/orga/projects/STARNET/ariella/Clint/SMAD3_annotatedSNPs.xls",SheetNames = c("AOR_cis","MAM_cis"),row.names =F)