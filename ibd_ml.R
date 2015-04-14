#isolation by distance - maximum likelihood estimation
#ibs analysis

library ("SNPRelate")
library(lattice)

#input
vcf.fn<-"~/Desktop/PSP_c88m4p3_indvmiss090_maf005_missing075.recode.vcf"#load the data
snpgdsVCF2GDS(vcf.fn,"testt.gds",method="copy.num.of.ref")#import to class. use method= biallelic.only for biallelic snps

#summary of the vcf imported
snpgdsSummary("testt.gds")

#open the translated vcf imported
genofile<-snpgdsOpen("testt.gds")

#import population names
pop_code <- scan("~/Desktop/Pop_fileas.txt", what=character())

#get the sample of interest
YRI.id<-sample.id[pop_code=="Var"]

# Estimate IBD coefficients
set.seed(1000)
snp.id <- sample(snpset.id, 1000)# random 5000 
SNPsibd <- snpgdsIBDMLE(genofile, sample.id=YRI.id, snp.id=snp.id,maf=0.05, missing.rate=0.05)
ibd.coeff<-snpgdsIBDSelection(SNPsibd)
plot(ibd.coeff$k0, ibd.coeff$k1,xlim=c(0,1),ylim=c(0,1),xlab="k0",ylab="k1",main="YRI samples (MoM)")
lines(c(0,1),c(1,0),col="red",lty=2)

#ibs analysis
ibs<-snpgdsIBS(genofile,num.thread=2)
L<-order(pop_code)
levelplot(ibs$ibs[L, L],col.regions= terrain.colors)

snpgdsClose(genofile)
