options(stringsAsFactors=F)
library(WriteXLS)

files = list.files(".","_trans_annotation.xls")

## merging the files together and adding some more annotations -- 

snpFileName = "/sc/orga/scratch/cohaia01/MGH/SnpsInterestedIn"
geneAnnotationFN = "/sc/orga/scratch/cohaia01/MGH/Expression/human gene annotations_v2.txt"

snps =  read.delim(snpFileName)
gene_info = read.delim(geneAnnotationFN)
snps$marker =paste(snps$chr,snps$pos,sep=":")

output = lapply(files,function(file){ #for(file in files){
	x = read.delim(paste(strsplit(file,".xls")[[1]][1],".txt",sep=""))
	gene_anno = t(sapply(1:nrow(x),function(X){gene_info[which(gene_info$reporter_id == x$gene[X])[1],]}))
	snp_anno = t(sapply(1:nrow(x),function(X){snps[which(snps$marker == x$snps[X])[1],]}))
	colnames(gene_anno) = colnames(gene_info)
	colnames(snp_anno) = colnames(snps)
	annotated = cbind(x,gene_anno,snp_anno)
	return(annotated)
	
})

names(output) = files


liver = output[[1]]
omental = output[[2]]
subq = output[[3]]

l = liver[!is.na(liver$reporter_id),]
o = omental[!is.na(omental$reporter_id),]
s = subq[!is.na(subq$reporter_id),]

WriteXLS(c("liver"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/Liver_annotatedSNPs_trans.xls",SheetNames = c("Liver"),row.names =F)
WriteXLS(c("omental"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/Omental_annotatedSNPs_trans.xls",SheetNames = c("Omental"),row.names =F)
WriteXLS(c("subq"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/Subq_annotatedSNPs_trans.xls",SheetNames = c("Subq_Fat"),row.names =F)

WriteXLS(c("l"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/Liver_annotatedSNPs_trans2.xls",SheetNames = c("Liver"),row.names =F)
WriteXLS(c("o"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/Omental_annotatedSNPs_trans2.xls",SheetNames = c("Omental"),row.names =F)
WriteXLS(c("s"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/Subq_annotatedSNPs_trans2.xls",SheetNames = c("Subq_Fat"),row.names =F)

liv = as.data.frame(l)
ome = as.data.frame(o)
sub = as.data.frame(s)

write.table(liv,file = "/sc/orga/scratch/cohaia01/MGH/Liver_annotatedSNPs_trans2.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(ome,file = "/sc/orga/scratch/cohaia01/MGH/Omental_annotatedSNPs_trans2.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(sub,file = "/sc/orga/scratch/cohaia01/MGH/Subq_annotatedSNPs_trans2.txt",sep="\t",row.names=F,col.names=T,quote=F)



WriteXLS(c("l","o","s"),ExcelFileName = "/sc/orga/scratch/cohaia01/MGH/annotatedSNPs_trans.xls",SheetNames = c("Liver","Omental","Subq_Fat"),row.names =F)

