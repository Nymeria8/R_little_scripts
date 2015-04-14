#VCF to PCA

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

#get sample ids
sample.id<-read.gdsn(index.gdsn(genofile,"sample.id"))

#data frame of eigenvectors
tab<-data.frame(sample.id= pca$sample.id,
                pop=factor(pop_code)[match(pca$sample.id, sample.id)],
                EV1= pca$eigenvect[,1],# first eignvector (can be changed)
                EV2= pca$eigenvect[,2],# the second eigenvector (can be changed)
                stringsAsFactors=FALSE)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)#creates space outside plot area
plot(tab$EV2, tab$EV1,  col=as.integer(tab$pop), xlab="eigenvector 2",ylab="eigenvector 1", pch=as.integer(tab$pop))
legend("topright", legend=levels(tab$pop),inset=c(-0.15,0), pch=1:nlevels(tab$pop), col=1:nlevels(tab$pop))#inset forces legend outside area

snpgdsClose(genofile)