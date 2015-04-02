#Run ebseq for differential expression analysis

library(EBSeq)

#Read the matrix of columns of values of counts - with the name of the genes and conditions - for 3 replicates
x<-read.table("/home/joanafino/Desktop/Link to Vitis_project/final_counts_genes/CGacT7.table", header=TRUE, row.names=1) 
rawData<-as.matrix(x)
conditions<-as.factor(c("C1","C1","C1","C2","C2","C2"))#Conditions and replicates of the table

Sizes=MedianNorm(rawData)#size factor estimation

EBOut=EBTest(Data=rawData, Conditions=conditions,sizeFactors=Sizes, maxround=5)#DE calculation

PP=GetPPMat(EBOut)#PPDE e PPEE calculation

GeneFC=PostFC(EBOut)#Fold change calculation
output <- matrix(unlist(GeneFC[1]))
b=cbind(PP,output)
resSig2<- b[b[, "PPDE"] >= .95,]#cutoff FDR=0.05

#write csv
write.csv (resSig2,file="/home/joanafino/Desktop/teste2")#Escreve resultados para ficheiro