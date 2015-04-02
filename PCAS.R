#Make a 2D/3D PCA

x<-read.table("soex_met.txt")#file input name
col<-c("C1T3","C2T3","C3T3","C1T5", "C2T5", "C3T5","C1T7","C2T7","C3T7")#name of the colums to be ploted
colnames(x)<-col

#PCA
xx<-prcomp(t(x))
dados<-xx$rotation
percentage<-xx$sdev^2/sum(xx$sdev^2)
percentage# explained variance

#gene/metabolite contribution (para fazer a percentagem soma-se a coluna e faz-se a percentagem m relacao a esse total) - write file
dados<-xx$rotation
write(dados[,1], file = "~/Desktop/data", ncolumns=3, append = FALSE, sep = "\t")

#plot coordinates 2D PCA
cordx=c(-0.1,0.1)
cordy=c(-0.1,0.1)

#plot 2D PCA
biplot(xx,var.axes=FALSE,cex=c(0.7),ylabs=NULL)#xlim=cordx, ylim=cordy)

#PCA 3D
library(rgl)
plot3d(xx$x,xlab="PC1",ylab="PC2",zlab="PC3",type="h")
spheres3d(xx$x, radius=0.5,col=rainbow(length(xx$x[,1])))
grid3d(side="z", at=list(z=0))
text3d(xx$x, text=rownames(xx$x), adj=1.3)