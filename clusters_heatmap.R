#Heatmap using two files, one with values and another with the row names

library("gplots")

#table with values
x<-read.table("/home/joanafino/Desktop/Link to Vitis_project/final_counts_genes/fpkms/Cvalues")
rawData<-as.matrix(x)#number of colums in the table

#label colums
col<-c("FS1-1","FS1-2", "FS1-3", "FS2-1","FS2-2", "FS2-3", "FS3-1","FS3-2", "FS3-3") #Column names
colnames(rawData)<-col

#filter by values >1
test<-do.call(rbind, apply(rawData, 1, function(r)if (! (all(r[1:9]<1)))r))

#log reduction (if needed)
aa<-apply(test, 1:2, function(x) log(x,2))

#stadarize data (mean=0 sd=1)
scaled.dat<-scale(aa, center=FALSE)

#confirm mean
apply(scaled.dat, 2, mean)

#confirm sd
apply(scaled.dat, 2, sd)

# perform hiearchical clustering and plot heatmap
heatmap.2(scaled.dat, density.info="none", trace="none", ,col=colorRampPalette(c("blue","black", "yellow")),symm=F,symkey=F, symbreaks=F, na.color="white", na.rm=TRUE)

