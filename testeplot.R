#Clustering by time

library(Mfuzz)

#load table
exprsFile<-file.path('profile.txt')
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1,as.is=TRUE))

#create expression set
eSet_xxx <- new("ExpressionSet", exprs = exprs)

#filtering missing data
tmp <- filter.std(eSet_xxx,min.std=0)

#standarizing: Since the clustering is performed in Euclidian space, the expression values of genes were
#standardised to have a mean value of zero and a standard deviation of one. This step ensures
#that vectors of genes with similar changes in expression are close in Euclidean space
mydataStandard<-standardise(tmp)

#m estimation
m1<-mestimate(mydataStandard)

#number of centers (arbitrary)
c1<-8 

#make clusters
cl <- mfuzz(mydataStandard,c=c1,m=m1)

#plot clusters
#mfrow: number of rows and columns
mfuzz.plot(mydataStandard,cl=cl,mfrow=c(2,4),time.labels=c("C","Gac","SH"))
dev.copy(svg, "~/Desktop/T7_shared.svg")#nome do ficheiro de output
dev.off()

#write file
write.csv (cl[3],file="/home/joanafino/Desktop/T7clusters_shared")
