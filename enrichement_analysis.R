#Enrichment Test by The Hypergeometric Distribution and FDR correction
#can be used with diferent categories
#read a table with the values, where:
#1 column - b - frequency of the term in the list
#2 column - d - frequency of the term in the genome
#3 column - e - frequency of the other terms in the genome
#4 column - a - number of differentially expressed genes in the list


#Read table with the values and rename the colums
x<-read.table("values")
col<-c("b","d","e","a")
colnames(x)<-col

#make hypergeometric distribution with phyper mode
results<-phyper(x$b,x$d,x$e,x$a,lower.tail = FALSE)

#correct the values with bonferroni/FDR
Qbal <- p.adjust(results, method="bonferroni")

#makes a matrix with the corrected and uncorrected values
s<-cbind(results,Qbal)

#gets the names of the rows and anex then to the matrix
r<-read.table("names")
row<-as.matrix(r)
rownames(s)<-row

#filter the values by test<0.05
filtred<- s[s[,2]< 0.05,]

#Write to file the filtred matrix
write.csv (filtred,file="output")

#Print the enriched categories
filtred
