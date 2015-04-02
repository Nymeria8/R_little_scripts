#Vens diagram creator with 3 datasets

library("grid")
library("VennDiagram")

area1<-57#size of dataset 1
area2<-113#size of dataset 2
area3<-914#size of dataset 3


#sizes of intersection
n12<-3
n23<-2
n13<-2
n123<-1

#names of the sets
cat=c("dataset1","datset2", "dataset3")#names

#graphic parameters
linecolor=c("red","seagreen", "blue")#color of the line
fillcolor=c("red", "lightgreen", "blue")#fill color
trans=rep(0.7, 3)#transparency

#Call function
draw.triple.venn(area1,area2,area3,n12,n23,n13,n123,category=cat, col=linecolor, fill=fillcolor, alpha =trans, cat.cex=rep(3,3))