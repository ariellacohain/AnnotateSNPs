options(stringsAsFactors=F)
library(WriteXLS)

files = list.files(".","annotation.xls")

## merging the files together and adding some more annotations -- 

snpFileName = "/sc/orga/scratch/cohaia01/MGH/SnpsInterestedIn"
geneAnnotationFN = "/sc/orga/scratch/cohaia01/MGH/Expression/human gene annotations_v2.txt"

snps =  read.delim(snpFileName)
gene_info = read.delim(geneAnnotationFN)

for(file in files){
	x = read.delim(paste(strsplit(file,".xls")[[1]][1],".txt",sep=" "))
	gene_anno = lapply(1:nrow(x),function(X){gene_infowhich(gene_info$reporter_id == x$gene[X],]})
	snp_anno = sapply(1:nrow(x),function(X){snps[snps$ == x$snps[X],]})

	annotated = cbind(x,gene_anno,snp_anno)
	write.table(annotated, file = paste(strsplit(".xls")[[1]][1],"2.txt",sep=""))
}

results = ldply(output,data.frame) #do.call(rbind,output)


liver = read.delim()
omental = read.delim()
subq = read.delim()

WriteXLS(c("liver","omental","subq"),ExcelOutputFN = "/sc/orga/scratch/cohaia01/MGH/annotatedSNPs.xls",SheetNames = c("Liver","Omental","Subq_Fat"),row.names =F)

