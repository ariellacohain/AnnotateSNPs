# grep -w "PDGFD" /sc/orga/projects/STARNET/ref/gencode.v19.annotation.transcript_level.gff
# head -1 /sc/orga/projects/STARNET/ariella/genotype.dose.11 > /sc/orga/projects/STARNET/ariella/Clint/11_103668962.dose
# grep -w "11:103668962" /sc/orga/projects/STARNET/ariella/genotype.dose.11  >> /sc/orga/projects/STARNET/ariella/Clint/11_103668962.dose

# grep -w "CDKN2B" /sc/orga/projects/STARNET/ref/gencode.v19.annotation.transcript_level.gff
# head -1 /sc/orga/projects/STARNET/ariella/genotype.dose.9 > /sc/orga/projects/STARNET/ariella/Clint/9_22103341.dose
# grep -w "9:22103341" /sc/orga/projects/STARNET/ariella/genotype.dose.9  >> /sc/orga/projects/STARNET/ariella/Clint/9_22103341.dose

# grep -w "LMOD1" /sc/orga/projects/STARNET/ref/gencode.v19.annotation.transcript_level.gff
# head -1 /sc/orga/projects/STARNET/ariella/genotype.dose.1 > /sc/orga/projects/STARNET/ariella/Clint/1_201886769.dose
# grep -w "1:201886769" /sc/orga/projects/STARNET/ariella/genotype.dose.1  >> /sc/orga/projects/STARNET/ariella/Clint/1_201886769.dose

# grep -w "IL6R" /sc/orga/projects/STARNET/ref/gencode.v19.annotation.transcript_level.gff
# head -1 /sc/orga/projects/STARNET/ariella/genotype.dose.1 > /sc/orga/projects/STARNET/ariella/Clint/1_154404336.dose
# grep -w "1:154404336" /sc/orga/projects/STARNET/ariella/genotype.dose.1  >> /sc/orga/projects/STARNET/ariella/Clint/1_154404336.dose

rm(list=ls())
options(stringsAsFactors=F)
library(MatrixEQTL)
library(ggplot2)
library(scales)

setwd("/sc/orga/projects/STARNET/ariella/Clint/")

# args <- commandArgs(trailingOnly=TRUE)  

# args <- c("PDGFD","11:103668962", "/sc/orga/projects/STARNET/ariella/Clint/11_103668962.dose",103777914,104035107)

# args <- c("CDKN2B","9:22103341", "/sc/orga/projects/STARNET/ariella/Clint/9_22103341.dose",22002902,22009362)
args <- c("LMOD1","1:201886769", "/sc/orga/projects/STARNET/ariella/Clint/1_201886769.dose", 201865580,201915715)
# args <- c("IL6R","1:154404336", "/sc/orga/projects/STARNET/ariella/Clint/1_154404336.dose",154377669,154441926)



gene = args[[1]]
snp = args[[2]]
chr = strsplit(snp,":",fixed=T)[[1]][1]
pos = strsplit(snp,":",fixed=T)[[1]][2]
gene_start =args[[4]]
gene_end = args[[5]]


load("/sc/orga/projects/STARNET/ariella/cis_trans_causlity/rerun/geneNames_all.Rdata")
ENS_gene = geneNames$id[geneNames$name == gene]

# expressionFN = paste("/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.",tissue,".exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN",sep="")
aor = read.delim("/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.AOR.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN",sep=" ",row.names=1)
mam = read.delim("/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.MAM.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN",sep=" ",row.names=1)

# aor_raw = 
# mam_raw = 

geno = read.delim(args[[3]])
geno_info = geno[,1:9]
geno_data = geno[,-c(1:9)]
geno_samples = cbind(1:ncol(geno_data),sapply(colnames(geno_data),function(x){strsplit(x,"X")[[1]][2]}))
colnames(geno_samples) = c("geno.idx","sample")

#matching AOR samples:
aor_samples = cbind(1:ncol(aor),sapply(colnames(aor),function(x){strsplit(x,"_",fixed=T)[[1]][2]}))
colnames(aor_samples) = c("gene.idx","sample")

samples_aor = merge(geno_samples,aor_samples,by=2)
samples_aor$gene.idx =as.numeric(as.character(samples_aor$gene.idx))
samples_aor$geno.idx =as.numeric(as.character(samples_aor$geno.idx))

aor_geno = round(unlist(geno_data[,samples_aor$geno.idx]))
aor_gene = unlist(aor[rownames(aor)==ENS_gene,samples_aor$gene.idx])

#matching MAM samples:
mam_samples = cbind(1:ncol(mam),sapply(colnames(mam),function(x){strsplit(x,"_",fixed=T)[[1]][2]}))
colnames(mam_samples) = c("gene.idx","sample")

samples_mam = merge(geno_samples,mam_samples,by=2)
samples_mam$gene.idx =as.numeric(as.character(samples_mam$gene.idx))
samples_mam$geno.idx =as.numeric(as.character(samples_mam$geno.idx))

mam_geno = round(unlist(geno_data[,samples_mam$geno.idx]))
mam_gene = unlist(mam[rownames(mam)==ENS_gene,samples_mam$gene.idx])

df = data.frame(expression = c(mam_gene,aor_gene), genotype = factor(c(mam_geno,aor_geno)),tissue = c(rep("MAM",length(mam_geno)),rep("AOR",length(aor_geno))))

pdf(paste(gene,"_boxplotsLM.pdf",sep=""))

x = ggplot(df, aes(x=genotype,y=expression,colour = tissue)) 	+ geom_point(shape=1) + theme_bw() + xlab(paste(snp,"Genotype Group")) +
	 ylab(paste(gene,"Normalized Expression")) + 
	 ggtitle(paste(gene, snp)) + 
	 geom_smooth(aes(group=1),method='lm')+ facet_wrap(~tissue) 

show(x)
dev.off()


pdf(paste(gene,"_boxplots.pdf",sep=""))

x = ggplot(df, aes(x=genotype,y=expression,colour = tissue)) + geom_boxplot() + theme_bw() + xlab(paste(snp,"Genotype Group")) + ylab(paste(gene,"Normalized Expression"))+ ggtitle(paste(gene, snp, "Boxplots"))

show(x)
dev.off()


## make Plots like SNP in region -- 
# args <- c("/sc/orga/projects/STARNET/ariella/Clint/snpsInterest","/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.MAM.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN",
# "/sc/orga/projects/STARNET/ariella/Clint/MAM_annotation_SMAD.xls")

# snpFileName = args[[1]] # "/sc/orga/projects/STARNET/ariella/Clint/snpsInterest" 
pdf(paste(gene,"_scatterplotInRegion.pdf",sep=""))
for(tissue in c("AOR","MAM")){
	expressionFileName = paste("/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.",tissue,".exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN",sep="")
	genotypeExpressionFileName = "/sc/orga/projects/STARNET/ariella/genotype.dose."

	threads = 10

	# snps = read.delim(snpFileName)
	# split_snps_by_chr = split(snps,snps$chr)

	geneAnnotationFN = "/sc/orga/projects/STARNET/ariella/MatrixEQTLgenes.txt"

	# SNP_file_name = paste(genotypeExpressionFileName[1],chr,genotypeExpressionFileName[2],sep="")
	snpSkipCol = 9
	snp.delimiter = "\t" 
	snp_omitChar = "NA"

	# snpPosFN = paste(genotypeExpressionFileName[1],chr,"_snpPos",sep="");
	# snpPosNameFN =  paste(genotypeExpressionFileName[1],chr,"_snpNames",sep="")

	expression_file_name = expressionFileName
	exprSkipCol = 1
	exprSkipRow = 1  
	expr_omitChar  = "NA"
	expr.delimiter = " "
	pvOutputThreshold_tra = 0;  ## NOT CALLING TRANS EQTLS!
	pvOutputThreshold_cis = 0.1;


	SNP_file_name = paste(genotypeExpressionFileName[1],chr,sep="")
	  
	snpPosFN = paste("/sc/orga/projects/STARNET/ariella/chr",chr,"_snps",sep="");
	snpPosNameFN =  paste("/sc/orga/projects/STARNET/ariella/chr",chr,"_snps.name",sep="");

	##Setting up Matrix eQTL for running all linear models:
	# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
	useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
	errorCovariance = numeric();
	cisDist = 1e6;

	## Load genotype data
	snpsME = SlicedData$new();
	snpsME$fileDelimiter = snp.delimiter;      # the TAB character
	snpsME$fileOmitCharacters = snp_omitChar; # denote missing values;
	snpsME$fileSkipRows = 1;          # one row of column labels
	snpsME$fileSkipColumns = snpSkipCol;       # one column of row labels
	snpsME$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	snpsME$LoadFile(SNP_file_name);

	## Load gene expression data
	geneME = SlicedData$new();
	geneME$fileDelimiter = expr.delimiter;      # the TAB character
	geneME$fileOmitCharacters = expr_omitChar; # denote missing values;
	geneME$fileSkipRows = exprSkipRow;          # one row of column labels
	geneME$fileSkipColumns = exprSkipCol;       # one column of row labels
	geneME$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	geneME$LoadFile(expression_file_name)

	##making sure samples are the same in SNP and Gene Expression:
	sample_snps= as.data.frame(cbind(colnames(snpsME))) 
	sample_snps = cbind(1:ncol(snpsME),sample_snps)
	names(sample_snps)=c("snp.idx","sample")
	sample_genes= as.data.frame(cbind(unlist(lapply(colnames(geneME),function(x){strsplit(x,"_",fixed=T)[[1]][2]}))))
	sample_genes = cbind(1:ncol(geneME),sample_genes)
	names(sample_genes)=c("gene.idx","sample")

	samples = merge(sample_snps,sample_genes,by=2)
	samples$gene.idx =as.numeric(as.character(samples$gene.idx))

	#taking only the samples that are in both data types:
	snpsME$ColumnSubsample( c(samples$snp.idx))
	geneME$ColumnSubsample( c(samples$gene.idx))

	## setting up the snp positions and gene positions
	snpspos = read.table(snpPosFN, header = T, stringsAsFactors = FALSE);
	snpspos.name = read.table(snpPosNameFN,header=TRUE,stringsAsFactors=FALSE)

	snpsposMarker = cbind(snpspos.name[,1],snpspos)
	names(snpsposMarker) = c("marker","chr","pos")
	snpsposMarker$chr= as.numeric(snpsposMarker$chr)
	snpsposMarker$pos= as.numeric(snpsposMarker$pos)

	if(ncol(snpsposMarker)==4){
	  snpsposMarker = snpsposMarker[,c('marker','chr','pos')]
	}

	#geneAnnotations Setting up for eQTL calling:
	gene_annotations1 = read.delim(geneAnnotationFN)
	genes = rownames(geneME)
	gene_annotation  = gene_annotations1[is.element(gene_annotations1$geneid,genes), ]
	colnames(gene_annotation) = c('reporter_id','chr','start_coord','end_coord')

	gene_annotation$start_coord = as.numeric(as.character(gene_annotation$start_coord))
	gene_annotation$end_coord = as.numeric(as.character(gene_annotation$end_coord))

	gene_annotation = gene_annotation[!is.na(gene_annotation$start_coord),]

	pvOutputThreshold_cis = 0.05
	pvOutputThreshold_tra = 0

	# Running Matrix EQTL: 
	me = Matrix_eQTL_main(
	  snps = snpsME, #matrix eqtl formmatted sliced data of snp of interest
	  gene = geneME, #all genes in matrix eqtl format
	  useModel = useModel, #linear model
	  errorCovariance = errorCovariance,
	  verbose = TRUE,
	  output_file_name     = NULL,#output_file_name_cis, #output_file_name_tra,
	  pvOutputThreshold     = pvOutputThreshold_tra,
	  output_file_name.cis = NULL,#output_file_name_tra,
	  pvOutputThreshold.cis = pvOutputThreshold_cis,
	  snpspos = snpsposMarker, # position of the snp of interest 
	  genepos = gene_annotation, # positions of all the genes
	  cisDist = cisDist,
	  pvalue.hist = "qqplot",
	  min.pv.by.genesnp = FALSE,
	  noFDRsaveMemory = FALSE
	  );

	eqtls = me$cis$eqtls

	data = eqtls[eqtls$gene == ENS_gene,]
	data$positions = sapply(data$snps,function(X){as.numeric(strsplit(X,":",fixed=T)[[1]][2])})
	data$color = "grey"
	data$color[data$snp == snp]="red"
	data$color[data$pvalue == min(data$pvalue) ]="blue"
	if(sum(data$snp == snp & data$pvalue == min(data$pvalue))==1){
		data$color[data$pvalue == min(data$pvalue)]="purple"
	}


	plot(data$positions,-log10(data$pvalue),pch =rep(20,nrow(data)), col = alpha(data$color,0.8),
		xlab=paste("Chromosome",chr,"(bp)"),ylab="-log10 p value",
		main = paste(gene,"eQTLs in Cis",tissue))
	# points(data$positions[data$color == "blue"],-log10(data$pvalue[data$color == "blue"]),col="blue",pch=20,cex =1)
	# points(data$positions[data$color == "red"],-log10(data$pvalue[data$color == "red"]),col="red",pch=20,cex =1)

	# abline(h=-log10(10e-3),col="black",lty = 2)
	abline(v = gene_start,col="black",lty=2)
	abline(v = gene_end,col="black",lty=2)
	legend("topright",c(paste("Snp of Interest"),"Max Snp"),col=c("red","blue"),pch =c(20,20),bty="n")

}
dev.off()
