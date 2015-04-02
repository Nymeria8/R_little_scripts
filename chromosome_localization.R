#Chromossome and Gene intron/exon organization in the vitis vinifera genome. Easly changed for other organisms from bioMart

library("GenomeGraphs")

#import genome
mart <- useMart("plants_mart_25", "vvinifera_eg_gene")#constant update

#atributes
ch=7 # chromossome
str=3460966 #start position
ed=3467045 #end position
st="+"#strand
cl="purple4"#color

#make chromossome location and intron exon organization
minStrand <- makeGeneRegion( chromosome =ch, start=str, end=ed, strand = st, biomart = mart)
ideogram <- makeIdeogram(chromosome =ch) 
setPar(minStrand, "protein_coding", cl)
setPar(ideogram, "size", "2")
gdPlot(list(ideogram, minStrand))
minStrand
