#VCF to PCA

library ("SNPRelate")

args <- commandArgs(trailingOnly = TRUE)

#Define inputs:
vcf.fn <- args[1] # "~/Desktop/sequenom/VCF/Qsuber_Sequenom.vcf"
gds.fn <- "/tmp/file.gds"
pops.fn <- args[2] # "~/Desktop/sequenom/VCF/pop_names_sorted.txt"
output_plot <- args[3] # "~/Desktop/sequenom/PCA/Qsuber_sequenom_PCA.svg"

#input
snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref") # import to class. use method= biallelic.only for biallelic snps

#summary of the vcf imported
snpgdsSummary(gds.fn)

#open the translated vcf imported
genofile<-snpgdsOpen(gds.fn)

#import population names
pop_code <- scan(pops.fn, what=character())
pop_code = gsub("_", " ", pop_code) # Replace any "_" with " ".

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
                EV1=pca$eigenvect[,1], # first eignvector (can be changed)
                EV2=pca$eigenvect[,2], # the second eigenvector (can be changed)
                stringsAsFactors=FALSE)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Creates space outside plot area

svg(filename=output_plot)

par(xpd=NA)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 0, 0, 9))

plot(tab$EV1, tab$EV2, col=as.integer(tab$pop),
     xlab="Eigenvector 1", ylab="Eigenvector 2",
     pch=as.integer(tab$pop))


legend("topright", legend=levels(tab$pop),
       inset=c(-0.45,0),  # Inset forces legend outside area
       pch=1:nlevels(tab$pop), col=1:nlevels(tab$pop))

dev.off()
snpgdsClose(genofile)