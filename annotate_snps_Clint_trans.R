# Rscript ~/thirdparty/AnnotateSNPs/annotate_snps.R /sc/orga/scratch/cohaia01/MGH/SnpsInterestedIn /sc/orga/scratch/cohaia01/MGH/Expression/MGH/subq_PCA.txt /sc/orga/scratch/cohaia01/MGH/subq_annotation.xls

options(stringsAsFactors=F)
library(foreach)
library(multicore)
library(doMC)
options(cores = multicore:::detectCores())
library(MatrixEQTL)
library(WriteXLS)

args <- commandArgs(trailingOnly=TRUE)  
#args <- c("/sc/orga/projects/STARNET/ariella/Clint/transInterestedIn","/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.MAM.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN",
  # "/sc/orga/projects/STARNET/ariella/Clint/MAM_trans_annotation.xls")

snpFileName = args[[1]] # "/sc/orga/projects/STARNET/ariella/Clint/snpsInterest" 
expressionFileName = args[[2]] #"/sc/orga/projects/STARNET/expression/normalized.cases/intermediates/STARNET.AOR.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN"
# genotypeExpressionFileName = c("/sc/orga/scratch/cohaia01/MGH/1000GinputeMay2013/chr",".dose.txt")
genotypeExpressionFileName = "/sc/orga/projects/STARNET/ariella/genotype.dose."
ExcelOutputFN = args[[3]] #"/sc/orga/projects/STARNET/ariella/Clint/AOR_annotation.xls"
OutputTableName = paste(strsplit(ExcelOutputFN,".xls",fixed=T)[[1]][1],".txt",sep="")
threads = 10

snps = read.delim(snpFileName)
split_snps_by_chr = split(snps,snps$chr)

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
pvOutputThreshold_tra = 0.05;  ## NOT CALLING TRANS EQTLS!
pvOutputThreshold_cis = 0.05;


#%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%

getConditionalResults<-function(cis_eqtl_data,snpNameInterest,geneME = geneME,snpME = snpInterest){
    # cis_eqtl_data = me_cis[me_cis$gene==cis.gene,]
    # snpInterest = paste(snp.of.interest[2],snp.of.interest[3],sep=":")

    if(sum(cis_eqtl_data$snps == snpNameInterest)==0){
      cat("GWAS SNP NOT PRESENT\n")
      return(NA)
    }else{
        cis_gene  = unique(cis_eqtl_data$gene)
        gwas_pval = cis_eqtl_data$pvalue[cis_eqtl_data$snp==snpNameInterest]
        most_sigSnp = cis_eqtl_data$snps[cis_eqtl_data$pvalue==min(cis_eqtl_data$pvalue,na.rm=T)]
        gene_expr = unlist(geneME$FindRow(cis_gene)$row)[1,]
        gwas_snp_expr = unlist(snpME$FindRow(snpNameInterest)$row)[1,]
        most_sig_snp_expr = unlist(snpME$FindRow(most_sigSnp[1])$row)[1,]
        # if(is.element(most_sigSnp,snpNameInterest)){
        #     res = cbind(cis_gene,gwas_pval, snpNameInterest, most_sigSnp[1],min(cis_eqtl_data$pvalue,na.rm=T), 0,1)
        #     colnames(res) = c("cis_gene","gwas_cis_pvalue","gwas_snp","most_sig_snp","sig_snp_pval","cond_pval","cond_rsqd")
        #     return(res)
        #   }else{
        lm_conditional = lm(gene_expr~most_sig_snp_expr,na.action = na.exclude)
        resid.vec = naresid(lm_conditional$na.action, lm_conditional$residuals)
        a = try(lm(resid.vec[!is.na(resid.vec) & !is.na(gwas_snp_expr) & resid.vec != "X" & gwas_snp_expr != "X"]~as.numeric(gwas_snp_expr[!is.na(resid.vec) & !is.na(gwas_snp_expr) & resid.vec != "X" & gwas_snp_expr != "X"])), TRUE)  

        if (!inherits(a, "try-error"))
        {
          a = summary(a)
          # cat(counter, " ", esnps.mod[counter,"gene"], "\t", esnps.mod[counter,"dbsnp_id"], "\t", a$r.squared, "\t", 1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1], "\n", sep="")
          cond_rsqd = a$r.squared
          cond_pval = 1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1]
        }else{
          cond_pval = NA
          cond_rsqd = NA
        }

        lm_conditional_MaxGvnGwas = lm(gene_expr~gwas_snp_expr,na.action = na.exclude)
        resid.vec_MaxGvnGwas = naresid(lm_conditional_MaxGvnGwas$na.action, lm_conditional_MaxGvnGwas$residuals)
        a_MaxGvnGwas = try(lm(resid.vec_MaxGvnGwas[!is.na(resid.vec_MaxGvnGwas) & !is.na(most_sig_snp_expr) & resid.vec_MaxGvnGwas != "X" & most_sig_snp_expr != "X"]~as.numeric(most_sig_snp_expr[!is.na(resid.vec_MaxGvnGwas) & !is.na(most_sig_snp_expr) & resid.vec_MaxGvnGwas != "X" & most_sig_snp_expr != "X"])), TRUE)  

        if (!inherits(a_MaxGvnGwas, "try-error"))
        {
          a_MaxGvnGwas = summary(a_MaxGvnGwas)
          # cat(counter, " ", esnps.mod[counter,"gene"], "\t", esnps.mod[counter,"dbsnp_id"], "\t", a$r.squared, "\t", 1-pf(a$fstatistic, a$df[1]-1, a$df[2])[1], "\n", sep="")
          cond_rsqd_MaxGvnGwas = a_MaxGvnGwas$r.squared
          cond_pval_MaxGvnGwas = 1-pf(a_MaxGvnGwas$fstatistic, a_MaxGvnGwas$df[1]-1, a_MaxGvnGwas$df[2])[1]
        }else{
          cond_pval_MaxGvnGwas = NA
          cond_rsqd_MaxGvnGwas = NA
        }

      
        res = cbind(cis_gene, gwas_pval, snpNameInterest, most_sigSnp[1],min(cis_eqtl_data$pvalue,na.rm=T), cond_pval,cond_rsqd, cond_pval_MaxGvnGwas,cond_rsqd_MaxGvnGwas)
        colnames(res) = c("cis_gene","gwas_cis_pvalue","gwas_snp","most_sig_snp","sig_snp_pval","cond_pval","cond_rsqd","cond_pval_MaxGvnGwas","cond_rsqd_MaxGvnGwas")
        
        ## calculating the means and standard errors for each of the genotype group:
        means = sapply(unique(round(gwas_snp_expr)),function(x){mean(gene_expr[round(gwas_snp_expr) == x],na.rm=T)})
        names(means) = paste("mean_",unique(round(gwas_snp_expr)),sep="")
        standardErr = sapply(unique(round(gwas_snp_expr)),function(x){sd(gene_expr[round(gwas_snp_expr) == x],na.rm=T)/sqrt(sum(!is.na(gene_expr[round(gwas_snp_expr) == x])))})
        names(standardErr) = paste("stderr_",unique(round(gwas_snp_expr)),sep="")
        res = c(res,means,standardErr)
        names(res) = c("cis_gene","gwas_cis_pvalue","gwas_snp","most_sig_snp","sig_snp_pval","cond_pval","cond_rsqd","cond_pval_MaxGvnGwas","cond_rsqd_MaxGvnGwas",
          paste("mean_",unique(round(gwas_snp_expr)),sep=""),paste("stderr_",unique(round(gwas_snp_expr)),sep=""))

        return(res)
        # }
    }

}


eqtlsOfInterest <- function(eqtls, snpList,geneME,snpME){
  e = eqtls[is.element(eqtls$snps,snpList),]
  if(nrow(e)<1){
    return("NA")
  }else{
    all_cond =c()
    for(i in 1:nrow(e)){
      eqtls_cis = eqtls[eqtls$gene == e$gene[i],]
      cond  = getConditionalResults(eqtls_cis,e$snps[i],geneME,snpME= snpME)
      all_cond = rbind(all_cond,cond)
    } 
    return(cbind(e,all_cond))
  }
  ## add condtional testing for each of those --

}


registerDoMC(threads) 

getMeansStdErrs<-function(eqtls,geneME,snpsME ){

  eqtls$mean_0 = NA
  eqtls$mean_1 = NA
  eqtls$mean_2 = NA
  eqtls$sterr_0 = NA
  eqtls$sterr_1 = NA
  eqtls$sterr_2 = NA
  for(snp in unique(eqtls$snps)){
    genotype = round(unlist(snpsME$FindRow(snp)$row)[1,])
    
    for(gene in unique(eqtls$gene[eqtls$snps == snp])){
      geneExpr = unlist(geneME$FindRow(gene)$row)[1,]  
      eqtls$mean_0[eqtls$snps == snp & eqtls$gene == gene] = mean(geneExpr[genotype == 0 ],na.rm=T)
      eqtls$mean_1[eqtls$snps == snp & eqtls$gene == gene] = mean(geneExpr[genotype == 1 ],na.rm=T)
      eqtls$mean_2[eqtls$snps == snp & eqtls$gene == gene] = mean(geneExpr[genotype == 2 ],na.rm=T)

      eqtls$sterr_0[eqtls$snps == snp & eqtls$gene == gene] = sd(geneExpr[genotype == 0 ],na.rm=T)/sqrt(sum(!is.na(geneExpr[genotype == 0 ])))
      eqtls$sterr_1[eqtls$snps == snp & eqtls$gene == gene] = sd(geneExpr[genotype == 1 ],na.rm=T)/sqrt(sum(!is.na(geneExpr[genotype == 1 ])))
      eqtls$sterr_2[eqtls$snps == snp & eqtls$gene == gene] = sd(geneExpr[genotype == 2 ],na.rm=T)/sqrt(sum(!is.na(geneExpr[genotype == 2 ])))
    }
  }

  return(eqtls)
  

}


{
output = foreach(chr  = names(split_snps_by_chr),combine= rbind)%dopar%{

  # cat("---on run",chr,"---\n")
  # tmp = paste(genotypeExpressionFileName[1],chr,genotypeExpressionFileName[2],sep="")
  # tmp2 = paste(genotypeExpressionFileName[1],chr,"_snpPos",sep="")
  # tmp3 = paste(genotypeExpressionFileName[1],chr,"_snpNames",sep="")
  # run_eqtlResults(chr, geneAnnotationFN,tmp,tmp2, tmp3 )

    # chr = names(split_snps_by_chr)[1]
    snpsInterest  = split_snps_by_chr[[as.character(chr)]]
    snpsInterest$marker = paste(snpsInterest$chr,snpsInterest$pos,sep=":")
    # geneAnnotationFN = "/sc/orga/scratch/cohaia01/MGH/Expression/human gene annotations_v2.txt"

    SNP_file_name = paste(genotypeExpressionFileName[1],chr,sep="")
    # snpSkipCol = 7 
    # snp.delimiter = "\t" 
    # snp_omitChar = "NA"

    snpPosFN = paste("/sc/orga/projects/STARNET/ariella/chr",chr,"_snps",sep="");
    snpPosNameFN =  paste("/sc/orga/projects/STARNET/ariella/chr",chr,"_snps.name",sep="");

    # expression_file_name = expressionFileName
    # exprSkipCol = 1
    # exprSkipRow = 1  
    # expr_omitChar  = "NA"
    # expr.delimiter="\t"
    # pvOutputThreshold_tra = 0;  ## NOT CALLING TRANS EQTLS!
    # pvOutputThreshold_cis = 10e-4;


    
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
    geneME$LoadFile(expression_file_name);


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

    trans_snp=c()
    for(snp in snpsInterest$marker){
      trans_snp_tmp =unlist(snpsME$FindRow(snp)$row)[1,]
      trans_snp=rbind(trans_snp,trans_snp_tmp)
    }
    rownames(trans_snp) = snpsInterest$marker
    snpInterest = SlicedData$new()  
    snpInterest$CreateFromMatrix(trans_snp)


    # Running Matrix EQTL: 
    me = Matrix_eQTL_main(
      snps = snpInterest, #matrix eqtl formmatted sliced data of snp of interest
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
   


    eqtls = me$trans$eqtls

    getMeansStdErrs(eqtls,geneME,snpsME)

    # chr_snps_of_interest = data.frame(split_snps_by_chr[as.character(chr)][[1]])
    # chr_snps_of_interest$marker = paste(chr_snps_of_interest$chr,chr_snps_of_interest$pos,sep=":")



    # eqtlsOfInterest(eqtls,chr_snps_of_interest$marker, geneME, snpsME)

    # gene_anno = t(sapply(1:nrow(e),function(x){gene_annotation[gene_annotation$reporter_id == e$gene[x],]}))
    # colnames(gene_anno) = colnames(gene_annotation)

    # cbind(e,gene_anno)

}
}
library(plyr)

results = ldply(output,data.frame) #do.call(rbind,output)

WriteXLS("results",ExcelFileName = ExcelOutputFN,
  SheetNames = "annotated_SNPs",row.names =F)

write.table(results,sep="\t",file = OutputTableName,row.names=F,quote=T)
