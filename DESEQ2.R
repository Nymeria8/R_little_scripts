#Run DESeq2 for differential expression analysis

library("DESeq2")

#matrix preparation (untreated /treated)
x<-read.table("/home/joanafino/Desktop/Link to Vitis_project/final_counts_genes/tables/GacSHT7.table", header=TRUE, row.names=1) #read the column matrix with the counts - with the gene name and condition name- for 3 replicates
samples <- data.frame(row.names=c("treated1","treated2","treated3","untreated1","untreated2","untreated3"), condition=as.factor(c(rep("treat",3),rep("untreat",3)))) #replicate names
dds<-DESeqDataSetFromMatrix(x, samples, formula(~ condition))#format matrix to deseq

#calculo dos DES
#normalizationFactors(dds) <- normFactors #faz os normalization factors not beter than size factors
dds <- estimateSizeFactors(dds) # use insted of the normfacotrs
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

#tratamento dos resultados
res<-results(dds)# submissao dos results
resSig <- res[which(res$padj < 0.05),]#cutoff de FDR=0.05
resSig2<- resSig[which(resSig$log2FoldChange > 1.5),]#cutoff de fold change
resSig3<- resSig[which(resSig$log2FoldChange < -1.5),]#cutoff de fold change

#DEs Statistics
nrow(res) #total DEs without fdr
nrow(resSig)#total DEs with fdr
nrow(resSig2)+nrow(resSig3)#total DEs with fold change cutoff
nrow(resSig2)#total upregulated DEs
nrow(resSig3)#total downregulated DEs

#Write file with the total DEs without fdr
write.csv (res,file="total")
#write file with the total cuted off DEs
write.csv (rbind(resSig2,resSig3),file="DEs")