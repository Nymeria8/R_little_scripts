#PCA snps e cenas (a descobrir)
library(gdsfmt)
library ("SNPRelate")

#input
vcf.fn<-"~/Desktop/Qsuber_indiv90_miss80_maf06.recode.vcf"#load the data
snpgdsVCF2GDS(vcf.fn,"test.gds",method="copy.num.of.ref")#import to class. use method= biallelic.only for biallelic snps

#summary of the vcf imported
snpgdsSummary("test.gds")

#open the translated vcf imported
genofile<-snpgdsOpen("test.gds")

#import population names
pop_code <- scan("~/Desktop/pops.txt", what=character())

#PCA
pca<-snpgdsPCA(genofile)

# variance proportion (%)
pc.percent<-pca$varprop*100
head(pc.percent)


sample.id<-read.gdsn(index.gdsn(genofile,"sample.id"))


tab<-data.frame(sample.id= pca$sample.id,
                pop=factor(pop_code)[match(pca$sample.id, sample.id)],
                EV1= pca$eigenvect[,2],# first eignvector
                EV2= pca$eigenvect[,3],# the second eigenvector
                stringsAsFactors=FALSE)

plot(tab$EV2, tab$EV1,  col=as.integer(tab$pop), xlab="eigenvector 2",ylab="eigenvector 1", pch=c(0,1,2,3,4,5))
legend("topleft", legend=levels(tab$pop), pch=c(0,1,2,3,4,5), col=1:nlevels(tab$pop))

cbind(sample.id, pop_code)

