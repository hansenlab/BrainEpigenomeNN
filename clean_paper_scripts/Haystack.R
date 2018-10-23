load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/assays-and-features.rda")
setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses")
library(GenomicRanges)
library(biomaRt)
library(XML)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

genome <- as(keepSeqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg19), paste0("chr", 1:22)), "GRanges")

promoters_use=reduce(unique(unstrand(promoters)))
export(promoters_use,con="promoters_by_gene.bed",format="bed") #background for promoter overlap
not_promoters=setdiff(genome,promoters_use)
export(not_promoters,con="not_promoters.bed",format="bed") #background for not_promoter overlap

hyper_DMRs_overlap_promoter=subsetByOverlaps(dmrs_NAvsBA9pos_hyper,promoters_use,ignore.strand=T) #2044
export(hyper_DMRs_overlap_promoter,con="hyper_DMRs_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hyper_DMRs_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename promoters_by_gene.bed
#new test
haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/hyper_DMRs_overlap_promoter_NAvBA9pos.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed


hypo_DMRs_overlap_promoter=subsetByOverlaps(dmrs_NAvsBA9pos_hypo,promoters_use,ignore.strand=T) #420
export(hypo_DMRs_overlap_promoter,con="hypo_DMRs_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hypo_DMRs_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename promoters_by_gene.bed

hyper_DAPs_overlap_promoter=subsetByOverlaps(daps_hyper,promoters_use,ignore.strand=T) #1155
export(hyper_DAPs_overlap_promoter,con="hyper_DAPs_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hyper_DAPs_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename promoters_by_gene.bed

hypo_DAPs_overlap_promoter=subsetByOverlaps(daps_hypo,promoters_use,ignore.strand=T) #826
export(hypo_DAPs_overlap_promoter,con="hypo_DAPs_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hypo_DAPs_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename promoters_by_gene.bed

non_promoter_hypo_DMRs=subsetByOverlaps(dmrs_NAvsBA9pos_hypo,promoters_use,invert=TRUE) #2891
export(non_promoter_hypo_DMRs,con="Hypo_DMRs_NAvBA9pos_not_in_promoter.bed",format="bed")
#haystack_motifs Hypo_DMRs_NAvBA9pos_not_in_promoter.bed hg19 --bed_bg_filename not_promoters.bed

non_promoter_hyper_DMRs=subsetByOverlaps(dmrs_NAvsBA9pos_hyper,promoters_use,invert=TRUE) #7540
export(non_promoter_hyper_DMRs,con="Hyper_DMRs_NAvBA9pos_not_in_promoter.bed",format="bed")
#haystack_motifs Hyper_DMRs_NAvBA9pos_not_in_promoter.bed hg19 --bed_bg_filename not_promoters.bed

non_promoter_hyper_DAPs=subsetByOverlaps(daps_hyper,promoters_use,invert=TRUE) #9529
export(non_promoter_hyper_DAPs,con="Hyper_DAPs_NAvBA9pos_not_in_promoter.bed",format="bed")
#haystack_motifs Hyper_DAPs_NAvBA9pos_not_in_promoter.bed hg19 --bed_bg_filename not_promoters.bed

non_promoter_hypo_DAPs=subsetByOverlaps(daps_hypo,promoters_use,invert=TRUE) #7816
export(non_promoter_hypo_DAPs,con="Hypo_DAPs_NAvBA9pos_not_in_promoter.bed",format="bed")
#haystack_motifs Hypo_DAPs_NAvBA9pos_not_in_promoter.bed hg19 --bed_bg_filename not_promoters.bed

#qrsh -l cegs,mem_free=100G,h_vmem=101G,h_fsize=200G
#cd /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/


#############  new Haystack groups  ################
#  take A = union(hypo dmrs and hyper daps) and B = union(hyper dmrs and hypo daps).
# Run Haystack on each of A and B and plot result

hyper_DMRs_hypo_DAPs_overlap_promoter=subsetByOverlaps(union(dmrs_NAvsBA9pos_hyper,daps_hypo),promoters_use,ignore.strand=T) #2618
export(hyper_DMRs_hypo_DAPs_overlap_promoter,con="hyper_DMRs_hypo_DAPs_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hyper_DMRs_hypo_DAPs_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename promoters_by_gene.bed

hypo_DMRs_hyper_DAPs_overlap_promoter=subsetByOverlaps(union(dmrs_NAvsBA9pos_hypo,daps_hyper),promoters_use,ignore.strand=T) #1435
export(hypo_DMRs_hyper_DAPs_overlap_promoter,con="hypo_DMRs_hyper_DAPs_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hypo_DMRs_hyper_DAPs_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename promoters_by_gene.bed

hyper_DMRs_hypo_DAPs_not_overlap_promoter=subsetByOverlaps(union(dmrs_NAvsBA9pos_hyper,daps_hypo),not_promoters,ignore.strand=T) #14463
export(hyper_DMRs_hypo_DAPs_not_overlap_promoter,con="hyper_DMRs_hypo_DAPs_not_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hyper_DMRs_hypo_DAPs_not_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename not_promoters.bed

hypo_DMRs_hyper_DAPs_not_overlap_promoter=subsetByOverlaps(union(dmrs_NAvsBA9pos_hypo,daps_hyper),not_promoters,ignore.strand=T) #11734
export(hypo_DMRs_hyper_DAPs_not_overlap_promoter,con="hypo_DMRs_hyper_DAPs_not_overlap_promoter_NAvBA9pos.bed",format="bed")
#haystack_motifs hypo_DMRs_hyper_DAPs_not_overlap_promoter_NAvBA9pos.bed hg19 --bed_bg_filename not_promoters.bed

##############  new Haystack groups  ###############

expressed=readRDS("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/rna-seq_expression_with_gene_symbols.rds")

expressed[expressed$gene_symbol == "TCF3", ]$exp

names=strsplit(degs$gene_id,"[.]")
names=sapply(names,"[[",1)
degs$symbol=names
x=getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene","description"),
  values= names,
  mart= grch37)
colnames(x)=c("symbol","external_gene_name" ,"entrezgene","description") 

y=match(degs$symbol,x$symbol)
symbol=x[y,]
degs$symbol=symbol$external_gene_name

#### 

hyper_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hyper_DMRs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_promoter=hyper_DMRs_promoter$motiftable
hyper_DMRs_promoter=hyper_DMRs_promoter[,1:7]
hyper_DMRs_promoter[,2]=toupper(hyper_DMRs_promoter[,2])
hyper_DMRs_promoter[,2]=replace(hyper_DMRs_promoter[,2],hyper_DMRs_promoter[,2]=="JUN(VAR.2)","JUN")
hyper_DMRs_promoter[,2]=replace(hyper_DMRs_promoter[,2],hyper_DMRs_promoter[,2]=="JDP2(VAR.2)","JDP2")

hypo_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hypo_DMRs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_promoter=hypo_DMRs_promoter$motiftable
hypo_DMRs_promoter=hypo_DMRs_promoter[,1:7]
hypo_DMRs_promoter[,2]=toupper(hypo_DMRs_promoter[,2])
hypo_DMRs_promoter[,2]=replace(hypo_DMRs_promoter[,2],hypo_DMRs_promoter[,2]=="JUN(VAR.2)","JUN")
hypo_DMRs_promoter[,2]=replace(hypo_DMRs_promoter[,2],hypo_DMRs_promoter[,2]=="JDP2(VAR.2)","JDP2")

hypo_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hypo_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_promoter=hypo_DAPs_promoter$motiftable
hypo_DAPs_promoter=hypo_DAPs_promoter[,1:7]
hypo_DAPs_promoter[,2]=toupper(hypo_DAPs_promoter[,2])
hypo_DAPs_promoter[,2]=replace(hypo_DAPs_promoter[,2],hypo_DAPs_promoter[,2]=="JUN(VAR.2)","JUN")
hypo_DAPs_promoter[,2]=replace(hypo_DAPs_promoter[,2],hypo_DAPs_promoter[,2]=="JDP2(VAR.2)","JDP2")

hyper_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hyper_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_promoter=hyper_DAPs_promoter$motiftable
hyper_DAPs_promoter=hyper_DAPs_promoter[,1:7]
hyper_DAPs_promoter[,2]=toupper(hyper_DAPs_promoter[,2])
hyper_DAPs_promoter[,2]=replace(hyper_DAPs_promoter[,2],hyper_DAPs_promoter[,2]=="JUN(VAR.2)","JUN")
hyper_DAPs_promoter[,2]=replace(hyper_DAPs_promoter[,2],hyper_DAPs_promoter[,2]=="JDP2(VAR.2)","JDP2")

#########
hypo_DAPs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hypo_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_not_promoter=hypo_DAPs_not_promoter$motiftable
hypo_DAPs_not_promoter=hypo_DAPs_not_promoter[,1:7]
hypo_DAPs_not_promoter[,2]=toupper(hypo_DAPs_not_promoter[,2])
hypo_DAPs_not_promoter[,2]=replace(hypo_DAPs_not_promoter[,2],hypo_DAPs_not_promoter[,2]=="JUN(VAR.2)","JUN")
hypo_DAPs_not_promoter[,2]=replace(hypo_DAPs_not_promoter[,2],hypo_DAPs_not_promoter[,2]=="JDP2(VAR.2)","JDP2")

hypo_DMRs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hypo_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_not_promoter=hypo_DMRs_not_promoter$motiftable
hypo_DMRs_not_promoter=hypo_DMRs_not_promoter[,1:7]
hypo_DMRs_not_promoter[,2]=toupper(hypo_DMRs_not_promoter[,2])
hypo_DMRs_not_promoter[,2]=replace(hypo_DMRs_not_promoter[,2],hypo_DMRs_not_promoter[,2]=="JUN(VAR.2)","JUN")
hypo_DMRs_not_promoter[,2]=replace(hypo_DMRs_not_promoter[,2],hypo_DMRs_not_promoter[,2]=="JDP2(VAR.2)","JDP2")

hyper_DAPs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hyper_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_not_promoter=hyper_DAPs_not_promoter$motiftable
hyper_DAPs_not_promoter=hyper_DAPs_not_promoter[,1:7]
hyper_DAPs_not_promoter[,2]=toupper(hyper_DAPs_not_promoter[,2])
hyper_DAPs_not_promoter[,2]=replace(hyper_DAPs_not_promoter[,2],hyper_DAPs_not_promoter[,2]=="JUN(VAR.2)","JUN")
hyper_DAPs_not_promoter[,2]=replace(hyper_DAPs_not_promoter[,2],hyper_DAPs_not_promoter[,2]=="JDP2(VAR.2)","JDP2")

hyper_DMRs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hyper_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_not_promoter=hyper_DMRs_not_promoter$motiftable
hyper_DMRs_not_promoter=hyper_DMRs_not_promoter[,1:7]
hyper_DMRs_not_promoter[,2]=toupper(hyper_DMRs_not_promoter[,2])
hyper_DMRs_not_promoter[,2]=replace(hyper_DMRs_not_promoter[,2],hyper_DMRs_not_promoter[,2]=="JUN(VAR.2)","JUN")
hyper_DMRs_not_promoter[,2]=replace(hyper_DMRs_not_promoter[,2],hyper_DMRs_not_promoter[,2]=="JDP2(VAR.2)","JDP2")

y=intersect(expressed$gene_symbol,hyper_DAPs_promoter[,2])
hyper_DAPs_promoter$Expressed="NO"
hyper_DAPs_promoter$Expressed[which(hyper_DAPs_promoter[,2] %in% y)] = "YES"
hyper_DAPs_promoter$DEG="NO"
hyper_DAPs_promoter$DEG[which(hyper_DAPs_promoter[,2] %in% degs$symbol)]="YES"

y=intersect(expressed$gene_symbol,hypo_DAPs_promoter[,2])
hypo_DAPs_promoter$Expressed="NO"
hypo_DAPs_promoter$Expressed[which(hypo_DAPs_promoter[,2] %in% y)] = "YES"
hypo_DAPs_promoter$DEG="NO"
hypo_DAPs_promoter$DEG[which(hypo_DAPs_promoter[,2] %in% degs$symbol)]="YES"

y=intersect(expressed$gene_symbol,hyper_DMRs_promoter[,2])
hyper_DMRs_promoter$Expressed="NO"
hyper_DMRs_promoter$Expressed[which(hyper_DMRs_promoter[,2] %in% y)] = "YES"
hyper_DMRs_promoter$DEG="NO"
hyper_DMRs_promoter$DEG[which(hyper_DMRs_promoter[,2] %in% degs$symbol)]="YES"

y=intersect(expressed$gene_symbol,hypo_DMRs_promoter[,2])
hypo_DMRs_promoter$Expressed="NO"
hypo_DMRs_promoter$Expressed[which(hypo_DMRs_promoter[,2] %in% y)] = "YES"
hypo_DMRs_promoter$DEG="NO"
hypo_DMRs_promoter$DEG[which(hypo_DMRs_promoter[,2] %in% degs$symbol)]="YES"

#########
y=intersect(expressed$gene_symbol,hyper_DAPs_not_promoter[,2])
hyper_DAPs_not_promoter$Expressed="NO"
hyper_DAPs_not_promoter$Expressed[which(hyper_DAPs_not_promoter[,2] %in% y)] = "YES"
hyper_DAPs_not_promoter$DEG="NO"
hyper_DAPs_not_promoter$DEG[which(hyper_DAPs_not_promoter[,2] %in% degs$symbol)]="YES"

y=intersect(expressed$gene_symbol,hypo_DAPs_not_promoter[,2])
hypo_DAPs_not_promoter$Expressed="NO"
hypo_DAPs_not_promoter$Expressed[which(hypo_DAPs_not_promoter[,2] %in% y)] = "YES"
hypo_DAPs_not_promoter$DEG="NO"
hypo_DAPs_not_promoter$DEG[which(hypo_DAPs_not_promoter[,2] %in% degs$symbol)]="YES"

y=intersect(expressed$gene_symbol,hyper_DMRs_not_promoter[,2])
hyper_DMRs_not_promoter$Expressed="NO"
hyper_DMRs_not_promoter$Expressed[which(hyper_DMRs_not_promoter[,2] %in% y)] = "YES"
hyper_DMRs_not_promoter$DEG="NO"
hyper_DMRs_not_promoter$DEG[which(hyper_DMRs_not_promoter[,2] %in% degs$symbol)]="YES"

y=intersect(expressed$gene_symbol,hypo_DMRs_not_promoter[,2])
hypo_DMRs_not_promoter$Expressed="NO"
hypo_DMRs_not_promoter$Expressed[which(hypo_DMRs_not_promoter[,2] %in% y)] = "YES"
hypo_DMRs_not_promoter$DEG="NO"
hypo_DMRs_not_promoter$DEG[which(hypo_DMRs_not_promoter[,2] %in% degs$symbol)]="YES"


#Make table of these results
install.packages("openxlsx", dependencies=TRUE)
library(openxlsx)
library(readr)

wb=createWorkbook()
addWorksheet(wb,"hyper_DMRs_overlap_promoter")
addWorksheet(wb,"hypo_DMRs_overlap_promoter")
addWorksheet(wb,"hypo_DAPs_overlap_promoter")
addWorksheet(wb,"hyper_DAPs_overlap_promoter")
writeDataTable(wb,"hyper_DMRs_overlap_promoter",x=hyper_DMRs_promoter)
writeDataTable(wb,"hypo_DMRs_overlap_promoter",x=hypo_DMRs_promoter)
writeDataTable(wb,"hyper_DAPs_overlap_promoter",x=hyper_DAPs_promoter)
writeDataTable(wb,"hypo_DAPs_overlap_promoter",x=hypo_DAPs_promoter)
saveWorkbook(wb, "Supplementary Table X – Haystack_Promoter_Overlap.xlsx", overwrite = TRUE)


wb=createWorkbook()
addWorksheet(wb,"hyper_DMRs_not_overlap_promoter")
addWorksheet(wb,"hypo_DMRs_not_overlap_promoter")
addWorksheet(wb,"hypo_DAPs_not_overlap_promoter")
addWorksheet(wb,"hyper_DAPs_not_overlap_promoter")
writeDataTable(wb,"hyper_DMRs_not_overlap_promoter",x=hyper_DMRs_not_promoter)
writeDataTable(wb,"hypo_DMRs_not_overlap_promoter",x=hypo_DMRs_not_promoter)
writeDataTable(wb,"hyper_DAPs_not_overlap_promoter",x=hyper_DAPs_not_promoter)
writeDataTable(wb,"hypo_DAPs_not_overlap_promoter",x=hypo_DAPs_not_promoter)

saveWorkbook(wb, "Supplementary Table X – Haystack_Not_Promoter_Overlap.xlsx", overwrite = TRUE)

Haystack_promoter=as.data.frame(cbind(hypo_DMRs_promoter,"HYPO_DMR"))
colnames(Haystack_promoter)[10]="Category"
Haystack_promoter2=as.data.frame(cbind(hyper_DMRs_promoter,"HYPER_DMR"))
colnames(Haystack_promoter2)[10]="Category"

Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hyper_DAPs_promoter,"HYPER_DAP"))
colnames(Haystack_promoter2)[10]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hypo_DAPs_promoter,"HYPO_DAP"))
colnames(Haystack_promoter2)[10]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)


Haystack_not_promoter=as.data.frame(cbind(hyper_DAPs_not_promoter,"HYPER_DAP"))
colnames(Haystack_not_promoter)[10]="Category"
Haystack_not_promoter2=as.data.frame(cbind(hypo_DAPs_not_promoter,"HYPO_DAP"))
colnames(Haystack_not_promoter2)[10]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)
Haystack_not_promoter2=as.data.frame(cbind(hypo_DMRs_not_promoter,"HYPO_DMR"))
colnames(Haystack_not_promoter2)[10]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)
Haystack_not_promoter2=as.data.frame(cbind(hyper_DMRs_not_promoter,"HYPER_DMR"))
colnames(Haystack_not_promoter2)[10]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)

Haystack_promoter$Type="Overlap_Promoter"

Haystack_not_promoter$Type="Not_Overlap_Promoter"

Haystack_list=rbind(Haystack_promoter,Haystack_not_promoter)
Haystack_list2=Haystack_list
for (i in 1:nrow(Haystack_list)){
if("MAX::MYC" %in% Haystack_list[i,1]){
	Haystack_list[i,5]=replace(Haystack_list[i,8],Haystack_list[i,8]=="NO","YES")
}}
Haystack_list2$log=-10*log10(as.numeric(as.character(Haystack_list[,7])))  #log of qvalue
#changed Inf to 2000 (highest number used)

rnames <- unique(expressed_Haystack[,1])                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(cbind(as.numeric(as.character(x[,2])),as.numeric(as.character(a[,2])),as.numeric(as.character(y[,2])),as.numeric(as.character(z[,2]))))  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames  

x=hyper_DMRs_not_promoter[,c(2,7)][which(hyper_DMRs_not_promoter$Expressed=="YES"),]
y=hypo_DMRs_not_promoter[,c(2,7)][which(hypo_DMRs_not_promoter$Expressed=="YES"),]
z=hyper_DAPs_not_promoter[,c(2,7)][which(hyper_DAPs_not_promoter$Expressed=="YES"),]
a=hypo_DAPs_not_promoter[,c(2,7)][which(hypo_DAPs_not_promoter$Expressed=="YES"),]

q=merge(x,a,by="Motif Name",all.x=TRUE,all.y=TRUE)
r=merge(y,z,by="Motif Name",all.x=TRUE,all.y=TRUE)
s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)
colnames(s)=c("TF","Hyper DMR","Hypo DAP","Hypo DMR","Hyper DAP")
s[,5]=-10*log10(as.numeric(as.character(s[,5])))

rnames <- s[,1]  
                          # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:ncol(s)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames  



x=hyper_DMRs_promoter[,c(2,5)][which(hyper_DMRs_promoter$Expressed=="YES"),]
y=hypo_DMRs_promoter[,c(2,5)][which(hypo_DMRs_promoter$Expressed=="YES"),]
z=hyper_DAPs_promoter[,c(2,5)][which(hyper_DAPs_promoter$Expressed=="YES"),]
a=hypo_DAPs_promoter[,c(2,5)][which(hypo_DAPs_promoter$Expressed=="YES"),]

q=merge(x,a,by="Motif Name",all.x=TRUE,all.y=TRUE)
r=merge(y,z,by="Motif Name",all.x=TRUE,all.y=TRUE)
s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)
colnames(s)=c("TF","Hyper DMR","Hypo DAP","Hypo DMR","Hyper DAP")
s[,5]=log10(as.numeric(as.character(s[,5])))
s[,4]=log10(as.numeric(as.character(s[,4])))
s[,3]=log10(as.numeric(as.character(s[,3])))
s[,2]=log10(as.numeric(as.character(s[,2])))


rnames <- s[,1]  
                          # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:ncol(s)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames  
mat_data[is.na(mat_data)] <- 0
mat_data2=log(mat_data)
mat_data2[is.infinite(mat_data2)] <- 0

expressed_Haystack=Haystack_list2[which(Haystack_list2$Expressed=="YES"),]
#expressed_Haystack2=expressed_Haystack
#expressed_Haystack=expressed_Haystack2
#expressed_Haystack$Ratio=as.numeric(as.character(expressed_Haystack$Ratio))

#expressed_Haystack3=as.matrix(expressed_Haystack)
rownames(expressed_Haystack)=NULL
Promoter=expressed_Haystack[1:101,]
Not_Promoter=expressed_Haystack[102:662,]
x=Not_Promoter[,c(2,12)][which(Not_Promoter$Category=="HYPO_DMR"),]
y=Not_Promoter[,c(2,12)][which(Not_Promoter$Category=="HYPER_DAP"),]
z=Not_Promoter[,c(2,12)][which(Not_Promoter$Category=="HYPER_DMR"),]
a=Not_Promoter[,c(2,12)][which(Not_Promoter$Category=="HYPO_DAP"),]

q=merge(x,a,by="Motif Name",all.x=TRUE,all.y=TRUE)
r=merge(y,z,by="Motif Name",all.x=TRUE,all.y=TRUE)
s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)
colnames(s)=c("TF","Hypo DMR","Hyper DAP","Hyper DMR","Hypo DAP")

rownames(s)=NULL

rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:ncol(s)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0

library(pheatmap)   
library(gplots)

mat_data2=mat_data[,c(1,3)]
mat_data2 = mat_data[ rowSums(mat_data)!=0, ] 
pdf("test_new8.pdf")
test=heatmap.2(mat_data2,
  #cellnote = mat_data, 
  main = "Promoter", 
  #notecol="none",     
  #density.info="none", 
  lhei = c(2, 8),
  trace="none", 
  na.color="white",cexRow=.1, cexCol=1,  
  margins =c(9,25),    
  col= viridis::viridis(n = 12),  
  scale="none",  
  breaks=breaks3,
  dendrogram="row")  
  #Colv="NA")  
dev.off()
Promoter
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  12.04   22.34   32.49   56.38   58.79  381.50

breaks <- seq(10, 200, length.out = 30)
breaks=c(0,breaks,seq(300,2000,length.out=10))
breaks=unique(breaks)
breaks2=seq(10,60,length.out=30)
breaks2=c(0,breaks2)
breaks2=unique(breaks2)

breaks3 <- seq(10, 500, length.out = 10)
breaks3=c(0,breaks3,seq(600,2000,length.out=2))
breaks3=unique(breaks3)


#if (nrow(table) > 100) stop("Too many rows for heatmap, who can read?!")
#fontsize_row = 10 - nrow(mat_data) / 15
#pdf("test2.pdf")
#pheatmap(mat_data,  main="NOT Promoters", cluster_cols=F, 
#         fontsize_row=fontsize_row, border_color=NA)
#dev.off()

#col=whiteblue(256),




colfunc<-colorRampPalette(c("white",my_palette,"black"))
my_palette
col_breaks = c(seq(0,1,length=1),  
  seq(1.01,2.5,length=150),           
  seq(2.51,3.5,length=100)) 
png("test.png",    # create PNG for the heat map        
  width = 5*300,        # 5 x 300 pixels
  height = 225*300,
  res=300,
  pointsize = 4)        # smaller font size
my_palette <- c("black",colorRampPalette(c("royalblue","springgreen","yellow","orange", "red")) (n=99))
breaks <- seq(1, 1.6, length.out = 90)
breaks=c(0,breaks,seq(1.61,3.46,length.out=10))
breaks=unique(breaks)
mat_data2=na.omit(mat_data)

#just DMRs
mat_data3=mat_data[,c(1,3)]
mat_data3[is.na(mat_data3)] <- 0
for (i in 1:nrow(mat_data3)){
if (sum(mat_data3[i,1],mat_data3[i,2])==0){mat_data3<-mat_data3[-i,]}}

mat_data3 = mat_data3[ rowSums(mat_data3)!=0, ] 


mat_data4=mat_data3[c(1:67,82:178),]

mat_data3[is.na(mat_data3)] <- 0
pdf("test7.pdf")
pdf("test_new.pdf")
test=heatmap.2(mat_data3,
  #cellnote = mat_data, 
  main = "DAPS", 
  #notecol="none",     
  #density.info="none", 
  lhei = c(2, 8),
  trace="none", 
  na.color="white",cexRow=0.2, cexCol=1,  
  margins =c(9,25),    
  col= viridis::viridis(n = 99),  
  scale="none",  
  breaks=breaks )
  #dendrogram="row",    
  #Colv="NA")  
dev.off()




h = heatmap(mat_data,keep.dendro=TRUE)
row.clusters = as.hclust( h$Rowv )
cutree(row.clusters,k=4)
row.clusters = hclust(dist(mat_data))

test=heatmap.2(mat_data)
y=cutree(as.hclust(test$rowDendrogram), k=4)
pdf("test.pdf")
heatmap(mat_data, Colv=F, scale='none')

rd=dist(expressed_Haystack)
h=hclust(rd)
pdf("test2.pdf")
heatmap(mat_data, Rowv=as.dendrogram(h), Colv=NA)
dev.off()
### Heatmap
# Load of libraries

library(reshape2)
library(ggplot2)

Promoter_Haystack=expressed_Haystack[which(expressed_Haystack$Type=="Overlap_Promoter"),]
Not_promoter_Haystack=expressed_Haystack[which(expressed_Haystack$Type=="Not_Overlap_Promoter"),]
# Elaboration of heatmap (white - steelblue)
test=rev(Promoter_Haystack)
pdf("test.pdf")
ggplot(Promoter_Haystack, aes(Promoter_Haystack[,1],Category)) +
  geom_tile(aes(fill = Ratio), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of TFs ") +
  xlab("Category") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Enrichment Ratio")
dev.off()

ggplot(df_heatmap, aes(patient, genes )) +
  geom_tile(aes(fill = expression_level), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of genes ") +
  xlab("List of patients") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")










write.csv(unique(expressed_Haystack[,2]),file="TF_list.csv")

AP1=c("JUND","JDP2","FOS","FOSL1","BATF3","JUNB","FOSL2","JUN")
AP1=expressed_Haystack[which(expressed_Haystack[,2] %in% AP1),]
AP1$group="AP-1"
Circ=c("BHLHE40","BHLHE41","DBP","CLOCK","NFIL3","ARNTL","ID2","NPAS2")
Circ=expressed_Haystack[which(expressed_Haystack[,2] %in% Circ),]
Circ$group="Circadian"
CRE=c("CREB3","CREB5","CREB3L1","CREBL2","CREB3L2","CREM","ATF7","ATF1","ATF3")
CRE=expressed_Haystack[which(expressed_Haystack[,2] %in% CRE),]
CRE$group="cAMP response element"
Ebox=c("MNT","MLX","HEY1","MLXIPL","MLXIP","MIXL1","NEUROD2","BHLHE22","MYC","MYCN","HES1","HES2","MAX")
Ebox=expressed_Haystack[which(expressed_Haystack[,2] %in% Ebox),]
Ebox$group="E-box"

EGR=c("EGR1","EGR2","EGR3","EGR4")
EGR=expressed_Haystack[which(expressed_Haystack[,2] %in% EGR),]
EGR$group="EGR"

E2F=c("E2F3","E2F4","E2F7","E2F8")
E2F=expressed_Haystack[which(expressed_Haystack[,2] %in% E2F),]
E2F$group="E2F"

ETS=c("ETV3","ETV5","FEV","ETV4","ETV1","FLI1","ERG","ELK3","ETS1","SPDEF","ELK4","SPIB","GABPA")
ETS=expressed_Haystack[which(expressed_Haystack[,2] %in% ETS),]
ETS$group="ETS"

homeo=c("MEIS2","MEIS3","MEIS1","PDX1","EMX2","MEOX2","CRX","EMX1","FOXG1","FOXO1","FOXO3","FOXO6","POU4F1","POU5F1B",
"VSX1","VAX1","DLX1","DLX2","DLX4","DLX6","HOXD3","HOXB3","MEOX1","POU6F1","HESX1","LBX2","NOTO","LHX6","GSX2",
"EN2","LHX2","LHX4","GBX2","POU6F2","MSX2","MSX1","POU4F2","HMX1","LHX8","RAX2","LMX1B","TGIF2","ALX4","PAX2","PAX6","POU2F3")
homeo=expressed_Haystack[which(expressed_Haystack[,2] %in% homeo),]
homeo$group="Homeobox"

MAF=c("MAFF","MAX","MAFK","NRL","MAFB")
MAF=expressed_Haystack[which(expressed_Haystack[,2] %in% MAF),]
MAF$group="MAF"

MEF2=c("MEF2A","MEF2C","MEF2D")
MEF2=expressed_Haystack[which(expressed_Haystack[,2] %in% MEF2),]
MEF2$group="MEF2"

NFAT=c("NFAT5","NFATC1","NFATC2","NFATC3")
NFAT=expressed_Haystack[which(expressed_Haystack[,2] %in% NFAT),]
NFAT$group="NFAT"

NF1=c("NFIA","NFIX","NFIL3","NFEL2","NFE2L2")
NF1=expressed_Haystack[which(expressed_Haystack[,2] %in% NF1),]
NF1$group="Nuclear Factor"

other=c("XBP1","TCF7L2","TCF21","TCFL5","MTF1","OLIG1","OLIG2","TFAP4","NHLH1","GCM2","NR3C1",
"NR3C2","HSF4","BCL6","ARID3B","BCL6B","CENPB","EBF1","HSF2","IRF8","NRF1","SOX4","GMEB1","MYB","NR2E3")
other=expressed_Haystack[which(expressed_Haystack[,2] %in% other),]
other$group="other"

PAR=c("HLF","TEF")
PAR=expressed_Haystack[which(expressed_Haystack[,2] %in% PAR),]
PAR$group="PAR bZIP"

RFX=c("RFX2","RFX3","RFX4","RFX5")
RFX=expressed_Haystack[which(expressed_Haystack[,2] %in% RFX),]
RFX$group="RFX"

STAT=c("STAT1","STAT3","STAT4","STAT6")
STAT=expressed_Haystack[which(expressed_Haystack[,2] %in% STAT),]
STAT$group="STAT"

ZN=c("ZIC4","ZIC1","ZNF740")
ZN=expressed_Haystack[which(expressed_Haystack[,2] %in% ZN),]
ZN$group="ZnFinger"

TEA=c("TEAD1","TEAD3")
TEA=expressed_Haystack[which(expressed_Haystack[,2] %in% TEA),]
TEA$group="TEA"


groups=rbind(AP1,Circ,CRE,Ebox,EGR,ETS,homeo,MAF,NF1,other,PAR,RFX,STAT,TEA,MEF2,NFAT,E2F,ZN)
rownames(groups)=NULL
groups=groups[order(groups[,2]),]

Haystack_groups=groups[order(groups$group),]

DEG_Haystack_groups=Haystack_groups[which(Haystack_groups$DEG=="YES"),]
DEG_Haystack_groups=DEG_Haystack_groups[order(DEG_Haystack_groups$group),]


y=match(DEG_Haystack_groups[,2],degs$symbol)
logFC=degs[y,]
logFC=as.data.frame(logFC)
logFC=logFC[,7]
DEG_Haystack_groups$logFC=logFC
write.csv(Haystack_groups,file="Haystack_groups.csv")
write.csv(DEG_Haystack_groups,file="DEG_Haystack_groups.csv")


#Plot expression of DEGs based on promoter overlap or not promoter overlap
DEG_Haystack_groups_promoter=DEG_Haystack_groups[which(DEG_Haystack_groups$Type=="Overlap_Promoter"),]
rownames(DEG_Haystack_groups_promoter)=NULL
DEG_Haystack_groups_promoter=DEG_Haystack_groups_promoter[rev(rownames(DEG_Haystack_groups_promoter)),]

gene=factor(DEG_Haystack_groups_promoter[,2],levels=unique(DEG_Haystack_groups_promoter[,2]))
library(ggplot2)
theme_set(theme_bw(base_size=12))
pdf("DEG_TF_plot_overlap_promoter2.pdf",width = 3, height = 3,useDingbats=FALSE)
ggplot(DEG_Haystack_groups_promoter, aes(x=gene, y=DEG_Haystack_groups_promoter[,13])) + ylim(-4,4)+
  geom_point(stat='identity') + ylab("log2(FC)") + xlab("")+geom_hline(yintercept = 0)+
  coord_flip() 
dev.off()


DEG_Haystack_groups_not_promoter=DEG_Haystack_groups[which(DEG_Haystack_groups$Type=="Not_Overlap_Promoter"),]
rownames(DEG_Haystack_groups_not_promoter)=NULL
DEG_Haystack_groups_not_promoter=DEG_Haystack_groups_not_promoter[rev(rownames(DEG_Haystack_groups_not_promoter)),]


gene2=factor(DEG_Haystack_groups_not_promoter[,2],levels=unique(DEG_Haystack_groups_not_promoter[,2]))

pdf("DEG_TF_plot_not_overlap_promoter2.pdf",width = 3, height = 9,useDingbats=FALSE)
ggplot(DEG_Haystack_groups_not_promoter, aes(x=gene2, y=DEG_Haystack_groups_not_promoter[,13])) + ylim(-6,6)+
  geom_point(stat='identity') + ylab("log2(FC)") + xlab("")+geom_hline(yintercept = 0)+
  coord_flip() 
dev.off()

DEG_Haystack_groups=DEG_Haystack_groups[rev(rownames(DEG_Haystack_groups)),]


gene=factor(DEG_Haystack_groups[,2],levels=unique(DEG_Haystack_groups[,2]))

pdf("DEG_TF_plot.pdf",width = 3, height = 9,useDingbats=FALSE)
ggplot(DEG_Haystack_groups, aes(x=gene, y=DEG_Haystack_groups[,13])) + ylim(-7,7)+
  geom_point(stat='identity') + ylab("log2(FC)") + xlab("")+geom_hline(yintercept = 0)+
  coord_flip() 
dev.off()
#########################

# checking list of expressed TFs to see how many are differentially methylated and 
# how many of THOSE are differentially expressed (and in what directions)
overlap_genes
overlap=match(overlap_genes,expressed$gene)
overlap=expressed[overlap,]
overlap=overlap$gene_symbol  #ok now overlap with TFs

check=match(expressed_Haystack[,2],overlap)
check2=expressed_Haystack[,2]
check2=na.omit(unique(overlap[check]))


enriched in both hypo-hyper DMRs:
ID2, NPAS2, ARNTL,MAX, MEIS1, RFX2, RFX3, PAX2, CREB3L1, TCFL5, ETV5,ELK3,E2F7,


expressed_Haystack[,2]


enriched in both DAPs:
JDP2, RFX2, RFX3, MEIS1, NFIA, MAFB,
NRL

TFs within DMRs

 [1] "JDP2"    "MEF2C"   "MEF2A"   "EBF1"    "MEIS2"   "OLIG2"   "ID2"    
 [8] "NPAS2"   "NFIA"    "MAX"     "NFIX"    "TEF"     "POU6F2"  "MEIS3"  
[15] "MEIS1"   "RFX2"    "RFX3"    "EGR1"    "RFX4"    "TGIF2"   "NRL"    
[22] "MAFB"    "MAFK"    "CREB5"   "EMX2"    "DLX1"    "BCL6"    "FLI1"   
[29] "NFATC1"  "NFATC2"  "TEAD1"   "NFAT5"   "PAX2"    "ATF7"    "ARNTL"  
[36] "CREB3L1" "DLX6"    "LHX2"    "TCFL5"   "LMX1B"   "ARID3B"  "ETV5"   
[43] "PAX6"    "ELK3"    "POU2F3"  "E2F7"    "TCF7L2"  "ETV1"    "GABPA"  
[50] "FOXO1"  


a=subsetByOverlaps(unlist(promoters_by_gene),dmrs_NAvsBA9pos)
a=subsetByOverlaps(unlist(unflattened_features_pc_transcripts$gene),dmrs_NAvsBA9pos_hyper)
a=subsetByOverlaps(unlist(unflattened_features_pc_transcripts$gene),dmrs_NAvsBA9pos_hypo)

a=unique(a$GENEID)
b=strsplit(unlist(a),'[.]')
b=sapply(b,"[[",1)
b=unique(b)


x=getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene","description"),
  values= b,
  mart= grch37)
colnames(x)=c("symbol","external_gene_name" ,"entrezgene","description") 
a=x$external_gene_name
c=match(expressed_Haystack[,2],a)
d=na.omit(expressed_Haystack[c,])
rownames(d)=NULL

d2=unique(d[,2])

7 TFs have DMR in promoter; 2 are DEGS (MEIS3 up and NEUROD2 down)
"NFATC1"-both  "TCF21"   "ELK4"    "MEIS3"-both   "POU2F3"-both  "FOSL2"   "NEUROD2" is downreg



b  # 20 of these 50 are DEGs

subsetByOverlaps(b,dmrs_NAvsBA9pos_hypo)  # 6 are hypo and upregulated
 MEIS2x
 MEIS1
  RFX4x
NFATC2
  DLX6x
 FOXO1x

subsetByOverlaps(b,dmrs_NAvsBA9pos_hyper)  # 16 are hyper and all but three are downregulated


MEIS1 (opp) NFATC2 (opp) has both DMRs
EBF1

### check directionality of DMRs with logFC (even tho not DEGs)
q=rna_atac_meth[, c("expLogFC", "DE", "gene_symbol")]
r=match(expressed_Haystack[,2],q$gene_symbol)
s=as.data.frame(q[r,])
up=s[which(s$expLogFC>0),]
down=s[which(s$expLogFC<0),]
hyper=b
hypo=b

x=getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene","description"),
  values= hypo,
  mart= grch37)
colnames(x)=c("symbol","external_gene_name" ,"entrezgene","description") 
hypo=x$external_gene_name
down_hyper=intersect(down$gene_symbol,hyper)  #27
up_hypo=intersect(up$gene_symbol,hypo)  #11


### how many differentially expressed in total?

diff_TF=intersect(degs$symbol, expressed_Haystack[,2]) #42/158 are DEGs
diff_TF=expressed_Haystack[which(expressed_Haystack[,2] %in% diff_TF),]
diff_TF=diff_TF[,c(2,10:11)]


#### look at list of open chrom regions for FOS from each list and do GO on nearest gene (promoter list)
hypo_FOS=as.data.frame(read.table("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hypo_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/motifs_regions/MA0476.1_motif_region_in_target.bed",stringsAsFactors=FALSE,header=TRUE,sep="\t",quote=""))
colnames(hypo_FOS)=c("chr","start","end","motifHits","numberHits")
hypo_FOS=GRanges(hypo_FOS)


hyper_FOS=as.data.frame(read.table("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hyper_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/motifs_regions/MA0476.1_motif_region_in_target.bed",stringsAsFactors=FALSE,header=TRUE,sep="\t",quote=""))
colnames(hyper_FOS)=c("chr","start","end","motifHits","numberHits")
hyper_FOS=GRanges(hyper_FOS)

hypo_FOS_promoters=unlist(subsetByOverlaps(promoters_by_gene,hypo_FOS))
hypo_FOS_promoters=unique(unlist(hypo_FOS_promoters$GENEID))
names=strsplit(hypo_FOS_promoters,"[.]")
names=sapply(names,"[[",1)

hyper_FOS_promoters=unlist(subsetByOverlaps(promoters_by_gene,hyper_FOS))
hyper_FOS_promoters=unique(unlist(hyper_FOS_promoters$GENEID))
names=strsplit(hyper_FOS_promoters,"[.]")
names=sapply(names,"[[",1)




x=getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene","description"),
  values= names,
  mart= grch37)
colnames(x)=c("symbol","external_gene_name" ,"entrezgene","description") 
hyper=x$external_gene_name

