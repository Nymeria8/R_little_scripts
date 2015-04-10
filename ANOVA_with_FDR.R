#ANOVA caluculation in a table with fdr

#tell where the data come from
datafilename="~/Desktop/metabolits"#table localization. rows are replicates and columns metabolits(the names should only contain numbers and leters)
data.ex1=read.table(datafilename, header=T)   #read the data into a table

#Input the names of the metabolits, again
d<-c("betahydroxypyruvate","betaine","glycine","serine","Smethylmethionine","Smethylcysteine","dihydoxyphenylalanine LDOPA","phenylalanine","quinate","shikimate","tryptophan","tyrosine","tyramine","phenethylamine isobar with 1phenylethanamine","2aminoadipate","2oxoadipate","alanine","asparagine","aspartate","glutarate pentanedioate","homoserine","lysine","methionine","methionine sulfoxide","pipecolate","Sadenosylhomocysteine SAH","threonine","6oxopiperidine2carboxylic acid","2hydroxyadipate","13diaminopropane","2aminobutyrate","4acetamidobutanoate","4hydroxybutyrate GHB","arginine","carboxyethylGABA","gammaaminobutyrate GABA","glutamate","glutamine","histidine","Nacetylproline","Nacetylputrescine","proline","pyroglutamine","stachydrine","trans4hydroxyproline","homostachydrine","Ndeltaacetylornithine","4guanidinobutanoate","isoleucine","23dihydroxyisovalerate","leucine","levulinate 4oxovalerate","valine","2isopropylmalate","5methylthioadenosine MTA","putrescine","spermidine","5oxoproline","glutathione oxidized GSSG","13dihydroxyacetone","fructose6phosphate","glucose","glucose6phosphate G6P","glycerate","Isobar fructose 16diphosphate glucose 16diphosphate myoinositol 14 or 13diphosphate","pyruvate","alphaketoglutarate","citrate","fumarate","isocitrate","malate","succinate","sedoheptulose7phosphate","glycolate hydroxyacetate","tartarate","2ketogulonate","arabonate","erythritol","erythronate","gluconate","ribitol","ribonate","ribose","ribulose","threitol","xylonate","xylose","xylulose","chiroinositol","inositol 1phosphate I1P","myoinositol","myoinositol hexakisphosphate","scylloinositol","myoinositol tetrakisphosphate 1346 or 3456 or 1345","3deoxyoctulosonate","erythrulose","fructose","galactinol","galactitol dulcitol","galactose","mannitol","mannose6phosphate","raffinose","sorbitol","sucrose","citramalate","2hydroxyglutarate","2hydroxypalmitate","2hydroxystearate","3hydroxypropanoate","azelate nonanedioate","cisvaccenate 181n7","linoleate 182n6","linolenate [alpha or gamma; 183n3 or 6]","oleate 181n9","undecanedioate","13HODE  9HODE","12propanediol","1palmitoylglycerol 1monopalmitin","glycerol","1linoleoylglycerophosphoinositol","1palmitoylglycerophosphate","1palmitoylglycerophosphocholine 160","1palmitoylglycerophosphoethanolamine","1palmitoylglycerophosphoinositol","glycerol 3phosphate G3P","glycerophosphorylcholine GPC","phosphoethanolamine","2linoleoylglycerophosphoinositol","1palmitoylglycerophosphoglycerol","1linolenoylglycerophosphoinositol","2linolenoylglycerophosphoinositol","choline phosphate","ethanolamine","phytosphingosine","betasitosterol","campesterol","stigmasterol","pantothenate","nicotinamide","nicotinamide adenine dinucleotide NAD","nicotinamide riboside","nicotianamine","nicotinate ribonucleoside","trigonelline N'methylnicotinate","methylphosphate","phosphate","carnitine","riboflavin Vitamin B2","flavin mononucleotide FMN","5ketogluconate","dehydroascorbate","threonate","alphatocopherol","betatocopherol","deltatocopherol","gammatocopherol","gammatocotrienol","pheophorbide A","adenine","adenosine","adenosine 5'monophosphate AMP","N6carbamoylthreonyladenosine","allantoin","guanine","guanosine","inosine","urate","xanthosine","cytidine","pseudouridine","uridine","betaalanine","gammaglutamylglutamine","gammaglutamylisoleucine","gammaglutamylleucine","gammaglutamylmethionine","gammaglutamylphenylalanine","gammaglutamyltryptophan","gammaglutamylvaline","abscisate","gibberellate","cyanoalanine","loganin","246trihydroxybenzoate","benzyl alcohol","benzoylOglucose","hydroquinone betaDglucopyranoside","galactarate mucic acid","kaempferol 7Oglucoside","apigenin7oglucoside","catechin","epicatechin","eriodictyol","dihydroquercetin","naringenin","naringenin7Oglucoside","procyanidin B1","procyanidin B2","quercetin","quercetin3galactoside","quercetin3oglucoside","rutin","kaempferol3rhamnoside","ferulate","gallate","resveratrol","salicylate","salidroside","salicin","caffeate","procyanidin trimer","catechin gallate","oleanolate","trizma acetate")


#ANOVA calculation for each metabolite
b<-c()#temporary file
#for loop to use each colum of the file individually
for(i in data.ex1[2:length(data.ex1)]){
  a<-cbind(data.ex1[1],i)#temporary matrix with FS/metabolit
  aov.ex1 = aov(i~FS,data=a)#anova calculation. summary(aov.ex1) gives degrees of freedom, sum sq, mean sq, F value and p-value.
  p_no_correction<-summary(aov.ex1)[[1]][["Pr(>F)"]][1]#select only p-value
  b<-c(b,p_no_correction) #add p-value to list
}

#fdr ajustment
Ajust <- p.adjust(b, method="fdr") #adjust the p-values with fdr. fdr can be replaced by bonferroni

#matrix formatation and write
matrix_final<-cbind(b,Ajust)#makes a matrix with the non corrected and corrected p-values
rownames(matrix_final)<-d #uses the metabolites names as row names
colnames(matrix_final)<- c("p-value", "fdr") #insert column names
write.table(matrix_final,file="anova_total_output") # write the total matrix to file as csv. can be opened in exel

#fdr cutoff
fdr_cut <- matrix_final[matrix_final[2,] < 0.05,]
write.table(fdr_cut,file="anova_fdr_output") # write the cuted of file as csv. can be opened in exel