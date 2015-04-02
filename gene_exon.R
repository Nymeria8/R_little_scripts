#Gene intron/exon organization in the vitis vinifera genome. Easly changed for other organisms from bioMart
library("GenomeGraphs")

#import genome
mart <- useMart("plants_mart_25", "vvinifera_eg_gene")#constant update

#Extract gene of iterest based on the ensemble gene id
gene <- makeGene(id = "VIT_04s0008g00330", type="ensembl_gene_id", biomart = mart)

showDisplayOptions(gene)#gene parameters (show)

#set colors
setPar(gene,"color","azure3")
setPar(gene, "idColor", "black")
setPar(gene,"plotId", "TRUE")
#plot gene
gdPlot(gene)
#print arguments of gene
gene