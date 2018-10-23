#convert bedGraph to bigwig for all pos samples
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
  # manually added lambda size 48502

fetchChromSizes to obtain the actual chrom.sizes information
from UCSC, please do not make up a chrom sizes from your own information.
The input bedGraph file must be sorted, use the unix sort command:
  sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph
options:
   -blockSize=N - Number of items to bundle in r-tree.  Default 256
   -itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024
   -unc - If set, do not use compression.
qrsh -l cegs,mem_free=80G,h_vmem=80G,h_fsize=200G

only need to do BA9 and NA since thats what i have atac for

## SORT BEDGRAPH FILES FIRST

#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq130/Sample_5552_NA_pos/5552_NA_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5552_NA_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq130/Sample_5569_NA_pos/5569_NA_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5569_NA_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5248_NA_pos/5248_NA_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5248_NA_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5284_BA9_pos/5284_BA9_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5284_BA9_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5552_BA9_pos/5552_BA9_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5552_BA9_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5570_BA9_pos/5570_BA9_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5570_BA9_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5248_BA9_pos/5248_BA9_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5248_BA9_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5284_NA_pos/5284_NA_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5284_NA_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5569_BA9_pos/5569_BA9_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5569_BA9_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5570_NA_pos/5570_NA_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5570_NA_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5628_BA9_pos/5628_BA9_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5628_BA9_pos.all.nsorted.bedGraph
#sort -k1,1 -k2,2n /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5628_NA_pos/5628_NA_pos.all.nsorted.bedGraph > /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5628_NA_pos.all.nsorted.bedGraph


#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq130/Sample_5552_NA_pos/5552_NA_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5552_NA_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq130/Sample_5569_NA_pos/5569_NA_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5569_NA_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5248_NA_pos/5248_NA_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5248_NA_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5284_BA9_pos/5284_BA9_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5284_BA9_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5552_BA9_pos/5552_BA9_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5552_BA9_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq137/Sample_5570_BA9_pos/5570_BA9_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5570_BA9_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5248_BA9_pos/5248_BA9_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5248_BA9_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5284_NA_pos/5284_NA_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5284_NA_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5569_BA9_pos/5569_BA9_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5569_BA9_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5570_NA_pos/5570_NA_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5570_NA_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5628_BA9_pos/5628_BA9_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5628_BA9_pos.bw
#bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq140/Sample_5628_NA_pos/5628_NA_pos.all.nsorted.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/5628_NA_pos.bw

qrsh -l mem_free=60G,h_vmem=60G,h_fsize=200G


# ATAC DATA is 
/dcl01/hansen/data/flow-sorted-brain-atac/data/bigWig/


expressed=readRDS("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/rna-seq_expression_with_gene_symbols.rds")
setwd( "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/")
x=expressed$exp
df <- data.frame(matrix(unlist(x), nrow=24161, byrow=T),stringsAsFactors=FALSE)
colnames(df)=colnames(x[[1]])
df$gene=expressed$gene_symbol

df$NAcc_pos=rowMeans(df[,c(2,4,6,9,11)])
df$BA9_pos=rowMeans(df[,c(1,3,5,7,8,10)])

BA9_pos=df[,c(12,14)]
colnames(BA9_pos)=NULL
write.table(BA9_pos,file="BA9_pos_expression.txt",sep="\t",quote=FALSE)



haystack_tf_activity_plane --output_directory /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/TF_activity_plane/ --name NAcc /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/TF_activity_plane/HAYSTACK_on_intersection_dars_dmrs_VS_random_background/ /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/Haystack_pipeline_samples_names.txt NA_pos






NAcc_pos=df[,c(12,13)]
colnames(NAcc_pos)=NULL
write.table(NAcc_pos,file="NAcc_pos_expression.txt",sep="\t",quote=FALSE)


BA9_pos_5343=df[,c(12,1)]
colnames(BA9_pos_5343)=NULL
write.table(BA9_pos_5343,file="5343_BA9_pos_exp.txt",sep="\t")

NAcc_pos_5343=df[,c(12,2)]
colnames(NAcc_pos_5343)=NULL
write.table(NAcc_pos_5343,file="5343_NAcc_pos_exp.txt",sep="\t")

BA9_pos_5347=df[,c(12,3)]
colnames(BA9_pos_5347)=NULL
write.table(BA9_pos_5347,file="5347_BA9_pos_exp.txt",sep="\t")

NAcc_pos_5347=df[,c(12,4)]
colnames(NAcc_pos_5347)=NULL
write.table(NAcc_pos_5347,file="5347_NAcc_pos_exp.txt",sep="\t")

BA9_pos_5358=df[,c(12,5)]
colnames(BA9_pos_5358)=NULL
write.table(BA9_pos_5358,file="5358_BA9_pos_exp.txt",sep="\t")

NAcc_pos_5358=df[,c(12,6)]
colnames(NAcc_pos_5358)=NULL
write.table(NAcc_pos_5358,file="5358_NAcc_pos_exp.txt",sep="\t")

BA9_pos_5404=df[,c(12,7)]
colnames(BA9_pos_5404)=NULL
write.table(BA9_pos_5404,file="5404_BA9_pos_exp.txt",sep="\t")

BA9_pos_5456=df[,c(12,8)]
colnames(BA9_pos_5456)=NULL
write.table(BA9_pos_5456,file="5456_BA9_pos_exp.txt",sep="\t")

NAcc_pos_5456=df[,c(12,9)]
colnames(NAcc_pos_5456)=NULL
write.table(NAcc_pos_5456,file="5456_NAcc_pos_exp.txt",sep="\t")

BA9_pos_5540=df[,c(12,10)]
colnames(BA9_pos_5540)=NULL
write.table(BA9_pos_5540,file="5540_BA9_pos_exp.txt",sep="\t")

NAcc_pos_5540=df[,c(12,11)]
colnames(NAcc_pos_5540)=NULL
write.table(NAcc_pos_5540,file="5540_NAcc_pos_exp.txt",sep="\t")



write.table(test3,file="/dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/Haystack_pipeline_atac_samples_names.txt",sep="\t",row.names=FALSE,quote=FALSE)


### ALL DATA 
samples_names.txt

1_NA_pos 5248_NA_pos.bw 5343_NAcc_pos.rep1.ATAC-seq.cpm.bw 5343_NAcc_pos_exp.txt
2_BA9_pos 5284_BA9_pos.bw 5347_BA9_pos.rep1.ATAC-seq.cpm.bw 5347_BA9_pos_exp.txt
2_NA_pos 5284_NA_pos.bw 5347_NAcc_pos.rep1.ATAC-seq.cpm.bw 5347_NAcc_pos_exp.txt
3_BA9_pos 5552_BA9_pos.bw 5358_BA9_pos.rep1.ATAC-seq.cpm.bw 5358_BA9_pos_exp.txt
3_NA_pos 5552_NA_pos.bw 5358_NAcc_pos.rep1.ATAC-seq.cpm.bw 5358_NAcc_pos_exp.txt
4_BA9_pos 5569_BA9_pos.bw 5456_BA9_pos.rep1.ATAC-seq.cpm.bw 5456_BA9_pos_exp.txt
4_NA_pos 5569_NA_pos.bw 5456_NAcc_pos.rep1.ATAC-seq.cpm.bw 5456_NAcc_pos_exp.txt
5_BA9_pos 5570_BA9_pos.bw 5540_BA9_pos.rep1.ATAC-seq.cpm.bw 5540_BA9_pos_exp.txt
5_NA_pos 5570_NA_pos.bw 5540_NAcc_pos.rep1.ATAC-seq.cpm.bw 5540_NAcc_pos_exp.txt
6_BA9_pos 5628_BA9_pos.bw 5404_BA9_pos.rep1.ATAC-seq.cpm.bw 5404_BA9_pos_exp.txt


To run full pipeline:
qrsh -l cegs,mem_free=100G,h_vmem=100G,h_fsize=200G
cd /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/
haystack_pipeline --input_is_bigwig --temp_directory /dcl01/feinberg/data/gtex/flow_sorted_nuclei/tmp --output_directory /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/Haystack_pipeline_atac_samples_names.txt hg19
# pipeline not working so run in pieces
# started with ATAC data
haystack_hotspots --input_is_bigwig  --output_directory /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/sample_names_hotspots.txt hg19
haystack_motifs myregions.bed hg19
# meth data
haystack_hotspots --input_is_bigwig --output_directory /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/Meth/HAYSTACK_PIPELINE_RESULTS /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/Meth/HAYSTACK_PIPELINE_RESULTS/sample_names_hotspots.txt hg19
haystack_motifs myregions.bed hg19



merge bigwigs together to do haystack_hotspots

bigWigMerge 5248_NA_pos.bw 5552_NA_pos.bw 5570_NA_pos.bw 5284_NA_pos.bw 5569_NA_pos.bw 5628_NA_pos.bw /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/NAcc.bedGraph
bigWigMerge 5248_BA9_pos.bw 5569_BA9_pos.bw 5284_BA9_pos.bw 5570_BA9_pos.bw 5552_BA9_pos.bw 5628_BA9_pos.bw /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/BA9.bedGraph
bigWigMerge 5343_NAcc_pos.rep1.ATAC-seq.cpm.bw 5358_NAcc_pos.rep1.ATAC-seq.cpm.bw 5540_NAcc_pos.rep1.ATAC-seq.cpm.bw 5347_NAcc_pos.rep1.ATAC-seq.cpm.bw 5456_NAcc_pos.rep1.ATAC-seq.cpm.bw /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/NAcc_ATAC.bedGraph
bigWigMerge 5343_BA9_pos.rep1.ATAC-seq.cpm.bw 5358_BA9_pos.rep1.ATAC-seq.cpm.bw 5456_BA9_pos.rep1.ATAC-seq.cpm.bw 5347_BA9_pos.rep1.ATAC-seq.cpm.bw 5404_BA9_pos.rep1.ATAC-seq.cpm.bw 5540_BA9_pos.rep1.ATAC-seq.cpm.bw /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/BA9_ATAC.bedGraph

bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/NAcc.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/NAcc_meth.bw
bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/BA9.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/BA9_meth.bw

bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/NAcc_ATAC.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/NAcc_ATAC.bw
bedGraphToBigWig /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/merged_bw/BA9_ATAC.bedGraph /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/hg19.chrom.sizes /dcl01/feinberg/data/gtex/flow_sorted_nuclei/bigwig/BA9_ATAC.bw


#haystack_motifs --name 1_NA /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_1_NA_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_1_NA_pos.500bp_z_0.25.bed hg19
#haystack_motifs --name 2_BA9 /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_2_BA9_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_2_BA9_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 2_NA /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_2_NA_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_2_NA_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 3_BA9 /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_3_BA9_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_3_BA9_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 3_NA /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_3_NA_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_3_NA_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 4_BA9 /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_4_BA9_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_4_BA9_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 4_NA /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_4_NA_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_4_NA_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 5_BA9 /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_5_BA9_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_5_BA9_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 5_NA /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_5_NA_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_5_NA_pos.500bp_z_0.25.bed hg19
haystack_motifs --name 6_BA9 /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Regions_specific_for_6_BA9_pos.500bp_z_1.50.bed --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/Haystack_pipeline_output/HAYSTACK_PIPELINE_RESULTS/HAYSTACK_HOTSPOTS/SPECIFIC_REGIONS/Background_for_6_BA9_pos.500bp_z_0.25.bed hg19


# orphan data
5248_BA9_pos.bw 
5343_BA9_pos.rep1.ATAC-seq.cpm.bw
5628_NA_pos.bw

hg19.chrom.sizes






What i should really do is just use the dmrs and dars that ive already assessed and just add in expression data








### Haystack with DARs
setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/")

DAR_DMR_hyper=intersect(dars_pos,dmrs_NAvsBA9pos_hyper)






# subset POS DMRs by those that overlap DAR
library(rtracklayer)
x=findOverlaps(dars_pos,dmrs_NAvsBA9pos)
dars_in_dmrs=dars_pos[queryHits(x)] 
dmrs_in_dars=dmrs_NAvsBA9pos[subjectHits(x)]
dars_in_dmrs$DMR_direction=dmrs_in_dars$direction  #181 regions do not agree
dars_in_dmrs=unique(dars_in_dmrs)

hyperMeth_dars_in_dmrs_all=dars_in_dmrs[which(dars_in_dmrs$DMR_direction=="hyper"),]
hypoMeth_dars_in_dmrs_all=dars_in_dmrs[which(dars_in_dmrs$DMR_direction=="hypo"),]
export(hyperMeth_dars_in_dmrs_all,con="hyperMeth_dars_in_dmrs_all.bed",format="bed")
export(hypoMeth_dars_in_dmrs_all,con="hypoMeth_dars_in_dmrs_all.bed",format="bed")

hyperMeth_dars_in_dmrs_promoter=subsetByOverlaps(hyperMeth_dars_in_dmrs_all,promoters) #1071
hypoMeth_dars_in_dmrs_promoter=subsetByOverlaps(hypoMeth_dars_in_dmrs_all,promoters) #328

export(hyperMeth_dars_in_dmrs_promoter,con="hyperMeth_dars_in_dmrs_promoter.bed",format="bed")
export(hypoMeth_dars_in_dmrs_promoter,con="hypoMeth_dars_in_dmrs_promoter.bed",format="bed")

#intersect
y=intersect(dars_pos,dmrs_NAvsBA9pos)
export(y,con="intersection_dars_dmrs.bed",format="bed")

#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/hyperMeth_dars_in_dmrs_all.bed hg19 --disable_ratio  --min_central_enrichment 0 
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/hypoMeth_dars_in_dmrs_all.bed hg19 --disable_ratio  --min_central_enrichment 0 
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/intersection_dars_dmrs.bed hg19 --disable_ratio  --min_central_enrichment 0 
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/hyperMeth_dars_in_dmrs_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/hypoMeth_dars_in_dmrs_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed


### DMRs that don't overlap DARs

dmrs_not_dars=subsetByOverlaps(dmrs_NAvsBA9pos,dars_pos,invert=TRUE)
dmrs_not_dars_hyper=dmrs_not_dars[which(dmrs_not_dars$direction=="hyper"),]
dmrs_not_dars_hypo=dmrs_not_dars[which(dmrs_not_dars$direction=="hypo"),]

export(dmrs_not_dars_hyper,con="dmrs_not_dars_hyper.bed",format="bed")
export(dmrs_not_dars_hypo,con="dmrs_not_dars_hypo.bed",format="bed")

#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/dmrs_not_dars_hyper.bed hg19 --disable_ratio  --min_central_enrichment 0 
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/dmrs_not_dars_hypo.bed hg19 --disable_ratio  --min_central_enrichment 0 

dmrs_not_dars_hyper_promoter=subsetByOverlaps(dmrs_not_dars_hyper,promoters)
dmrs_not_dars_hypo_promoter=subsetByOverlaps(dmrs_not_dars_hypo,promoters)
export(dmrs_not_dars_hyper_promoter,con="dmrs_not_dars_hyper_promoter.bed",format="bed")
export(dmrs_not_dars_hypo_promoter,con="dmrs_not_dars_hypo_promoter.bed",format="bed")
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/dmrs_not_dars_hyper_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/dmrs_not_dars_hypo_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed

dmrs_not_dars_hyper_not_promoter=subsetByOverlaps(dmrs_not_dars_hyper,promoters,invert=TRUE)
dmrs_not_dars_hypo_not_promoter=subsetByOverlaps(dmrs_not_dars_hypo,promoters,invert=TRUE)
export(dmrs_not_dars_hyper_not_promoter,con="dmrs_not_dars_hyper_not_promoter.bed",format="bed")
export(dmrs_not_dars_hypo_not_promoter,con="dmrs_not_dars_hypo_not_promoter.bed",format="bed")
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/dmrs_not_dars_hyper_not_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/not_promoters.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/dmrs_not_dars_hypo_not_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/not_promoters.bed




setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses")
library(GenomicRanges)
library(biomaRt)
library(XML)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gplots)
expressed=readRDS("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/rna-seq_expression_with_gene_symbols.rds")

### Redoing Haystack to get depletion values

# qrsh -l cegs,mem_free=80G,h_vmem=81G,h_fsize=200G
# cd /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/

#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/hyper_DMRs_overlap_promoter_NAvBA9pos.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/hypo_DMRs_overlap_promoter_NAvBA9pos.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/hyper_DAPs_overlap_promoter_NAvBA9pos.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/hypo_DAPs_overlap_promoter_NAvBA9pos.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/promoters_by_gene.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/Hypo_DMRs_NAvBA9pos_not_in_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/not_promoters.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/Hyper_DMRs_NAvBA9pos_not_in_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/not_promoters.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/Hyper_DAPs_NAvBA9pos_not_in_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/not_promoters.bed
#haystack_motifs /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/Hypo_DAPs_NAvBA9pos_not_in_promoter.bed hg19 --disable_ratio  --min_central_enrichment 0 --bed_bg_filename /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/not_promoters.bed

library(XML)
hyper_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_hyper_DMRs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_promoter=hyper_DMRs_promoter$motiftable
hyper_DMRs_promoter=hyper_DMRs_promoter[,1:8]
hyper_DMRs_promoter[,2]=toupper(hyper_DMRs_promoter[,2])
hyper_DMRs_promoter[,8]=as.numeric(hyper_DMRs_promoter[,8])
hyper_DMRs_promoter[,7]=as.numeric(hyper_DMRs_promoter[,7])

hypo_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_hypo_DMRs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_promoter=hypo_DMRs_promoter$motiftable
hypo_DMRs_promoter=hypo_DMRs_promoter[,1:8]
hypo_DMRs_promoter[,2]=toupper(hypo_DMRs_promoter[,2])
hypo_DMRs_promoter[,8]=as.numeric(hypo_DMRs_promoter[,8])
hypo_DMRs_promoter[,7]=as.numeric(hypo_DMRs_promoter[,7])

hypo_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_hypo_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_promoter=hypo_DAPs_promoter$motiftable
hypo_DAPs_promoter=hypo_DAPs_promoter[,1:8]
hypo_DAPs_promoter[,2]=toupper(hypo_DAPs_promoter[,2])
hypo_DAPs_promoter[,8]=as.numeric(hypo_DAPs_promoter[,8])
hypo_DAPs_promoter[,7]=as.numeric(hypo_DAPs_promoter[,7])

hyper_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_hyper_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_promoter=hyper_DAPs_promoter$motiftable
hyper_DAPs_promoter=hyper_DAPs_promoter[,1:8]
hyper_DAPs_promoter[,2]=toupper(hyper_DAPs_promoter[,2])
hyper_DAPs_promoter[,8]=as.numeric(hyper_DAPs_promoter[,8])
hyper_DAPs_promoter[,7]=as.numeric(hyper_DAPs_promoter[,7])

Haystack_promoter=as.data.frame(cbind(hypo_DMRs_promoter,"HYPO_DMR"))
colnames(Haystack_promoter)[9]="Category"
Haystack_promoter2=as.data.frame(cbind(hyper_DMRs_promoter,"HYPER_DMR"))
colnames(Haystack_promoter2)[9]="Category"

Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hyper_DAPs_promoter,"HYPER_DAP"))
colnames(Haystack_promoter2)[9]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hypo_DAPs_promoter,"HYPO_DAP"))
colnames(Haystack_promoter2)[9]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)
Haystack_promoter[,7]=as.numeric(Haystack_promoter[,7])
Haystack_promoter$sig="NO"
Haystack_promoter$sig[which(Haystack_promoter[,7]<=0.01)]="YES"
# with 0.05 cutoff
# NO YES 
# 29 254

# with 0.01 cutoff
# NO YES 
# 117 166

Haystack_promoter[,8]=log10(Haystack_promoter[,8])
Haystack_promoter[,8][Haystack_promoter[,7]>0.01]<-0 # 283


x=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPO_DMR"),]
y=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPER_DAP"),]
z=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPER_DMR"),]
a=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPO_DAP"),]

q=merge(x,y,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPO_DMR","HYPER_DAP"))
r=merge(z,a,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPER_DMR","HYPO_DAP"))

s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)
colnames(s)=c("TF","HYPO_DMR","HYPER_DAP","HYPER_DMR","HYPO_DAP")

rownames(s)=NULL


t=intersect(expressed$gene_symbol,s[,1])
s$Expressed="NO"
s$Expressed[which(s[,1] %in% t)] = "YES"

#Manually change to yes
s$Expressed[c(3,8,9,55,85,88,89,97,130,133,150,152,154,162,170)]="YES"

meth_TF=read.csv("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/meth_effect_TFBinding_Science2017_table.csv",stringsAsFactors=FALSE,header=TRUE)
colnames(meth_TF)=c("TF","MethSens")
meth_TF$MethSensitive="yes"
meth_TF$MethSensitive[which(meth_TF$MethSens %in% c("inconclusive","Little effect"))] = "no"

test=s$TF[which(!s$TF %in% meth_TF$TF)]
 test=as.data.frame(test,"Undetermined")
 test$MethSens="Undetermined"
 test$MethSensitive="no"
 colnames(test)=c("TF","MethSens","MethSensitive")

 test$MethSensitive[which(test$TF=="FOS::JUN")]="yes"
 test$MethSens[which(test$TF=="FOS::JUN")]="MethylMinus"

test$MethSensitive[which(test$TF=="BATF::JUN")]="yes"
 test$MethSens[which(test$TF=="BATF::JUN")]="MethylMinus"

test$MethSensitive[which(test$TF=="MAX::MYC")]="yes"
 test$MethSens[which(test$TF=="MAX::MYC")]="MethylMinus"
   
test$MethSensitive[which(test$TF=="JUND(VAR.2)")]="yes"
 test$MethSens[which(test$TF=="JUND(VAR.2)")]="MethylMinus"

test$MethSensitive[which(test$TF=="JUN(VAR.2)")]="yes"
 test$MethSens[which(test$TF=="JUN(VAR.2)")]="MethylMinus"

test$MethSensitive[which(test$TF=="NR2F6(VAR.2)")]="yes"
 test$MethSens[which(test$TF=="NR2F6(VAR.2)")]="MethylPlus"
test$MethSensitive[which(test$TF=="RARG(VAR.2)")]="yes"
 test$MethSens[which(test$TF=="RARG(VAR.2)")]="MethylPlus"
test$MethSensitive[which(test$TF=="RARA(VAR.2)")]="yes"
 test$MethSens[which(test$TF=="RARA(VAR.2)")]="MethylPlus"
test$MethSensitive[which(test$TF=="RARB(VAR.2)")]="yes"
 test$MethSens[which(test$TF=="RARB(VAR.2)")]="MethylPlus"


meth_TF2=rbind(meth_TF,test)

x=merge(s,meth_TF2,by="TF")


x2=x[which(x$Expressed=="YES"),]
rnames <- x2[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(x2[,2:5])  
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0 

s=s[which(s$Expressed=="YES"),]
rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:5])  
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0 





summary(Haystack_promoter[,8])
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.31876  0.00000  0.00000  0.05706  0.13354  0.51587 

breaks <- seq(-0.5,0.5, length.out = 40)
library(circlize)
library(RColorBrewer)
cols=colorRampPalette(c("#A50026", "#D73027", "#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))(length(breaks) - 1)

mat_data_dmr=mat_data[,c(1,3)]
mat_data_dmr=mat_data_dmr[which(rowSums(mat_data_dmr)!=0),]
x3=x2[which(x2[,1] %in% rownames(mat_data_dmr)),]
#column_tree = hclust(dist(t(mat_data)))

pdf("testing_TFBS.pdf")
Heatmap(mat_data_dmr, name = "log10(central enrichment)", col = rev(cols),
     column_title = "Enrichment", row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 10)) +
Heatmap(x3[,7], column_title = "Methylation Sensitivity", column_title_gp = gpar(fontsize = 10), col = structure(brewer.pal(length(unique(x3[,7])), "Set1"), names = unique(x3[,7])))+
Heatmap(x3[,8], column_title = "Methylation Sensitivity", column_title_gp = gpar(fontsize = 10), col = c("white","red")) 
dev.off()

mat_data2=mat_data[which(rowSums(mat_data)!=0),]
x3=x2[which(x2[,1] %in% rownames(mat_data2)),]
breaks <- seq(-0.3,0.3, length.out = 40)
cols=colorRampPalette(c("#A50026", "#D73027", "#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))(length(breaks) - 1)

pdf("Promoter_TFBS_Haystack.pdf")
Heatmap(mat_data2, name = "log10(central enrichment)", col = rev(cols),
     column_title = "Enrichment", row_names_gp = gpar(fontsize = 4), column_title_gp = gpar(fontsize = 10)) +
Heatmap(x3[,7],  column_title_gp = gpar(fontsize = 10), col = structure(brewer.pal(length(unique(x3[,7])), "Set1"), names = unique(x3[,7])))+
Heatmap(x3[,8], column_title = "Methylation Sensitivity", column_title_gp = gpar(fontsize = 10), col = c("white","red")) 
dev.off()

### two plots one for methylsensitive and one for not methyl sensitive

x2=x[which(x$Expressed=="YES"),]
x2_yes=x2[which(x2$MethSensitive == "yes"),]
rnames <- x2_yes[,1]                            # assign labels in column 1 to "rnames"
mat_data_yes <- data.matrix(x2_yes[,2:5])  
rownames(mat_data_yes) <- rnames 
mat_data_yes[is.na(mat_data_yes)] <- 0 

x2_no=x2[which(x2$MethSensitive == "no"),]
rnames <- x2_no[,1]                            # assign labels in column 1 to "rnames"
mat_data_no <- data.matrix(x2_no[,2:5])  
rownames(mat_data_no) <- rnames 
mat_data_no[is.na(mat_data_no)] <- 0 


mat_data_yes2=mat_data_yes[which(rowSums(mat_data_yes)!=0),]
x3_yes=x2_yes[which(x2_yes[,1] %in% rownames(mat_data_yes2)),]

mat_data_no2=mat_data_no[which(rowSums(mat_data_no)!=0),]
x3_no=x2_no[which(x2_no[,1] %in% rownames(mat_data_no2)),]

pdf("Promoter_TFBS_Haystack_sep.pdf")
Heatmap(mat_data_yes2, name = "log10(central enrichment)", col = rev(cols),
     column_title = "Enrichment", row_names_gp = gpar(fontsize = 4), column_title_gp = gpar(fontsize = 10)) +
Heatmap(x3_yes[,7],  column_title_gp = gpar(fontsize = 10), col = structure(brewer.pal(length(unique(x3_yes[,7])), "Set1"), names = unique(x3_yes[,7])))

Heatmap(mat_data_no2, name = "log10(central enrichment)", col = rev(cols),
     column_title = "Enrichment", row_names_gp = gpar(fontsize = 4), column_title_gp = gpar(fontsize = 10)) +
Heatmap(x3_no[,7],  column_title_gp = gpar(fontsize = 10), col = structure(brewer.pal(length(unique(x3_no[,7])), "Set1"), names = unique(x3_no[,7])))
dev.off()

### just look at enriched not depleted...
mat_data2=mat_data[which(rowSums(mat_data)!=0),]

x3=x2[which(x2[,1] %in% rownames(mat_data2)),]
breaks <- seq(0,0.6, length.out = 40)
cols=colorRampPalette(c("white","green"))(length(breaks) - 1)

pdf("Promoter_TFBS_Haystack2.pdf")
Heatmap(mat_data2, name = "log10(central enrichment)", col = rev(cols),
     column_title = "Enrichment", row_names_gp = gpar(fontsize = 4), column_title_gp = gpar(fontsize = 10)) +
Heatmap(x3[,7],  column_title_gp = gpar(fontsize = 10), col = structure(brewer.pal(length(unique(x3[,7])), "Set1"), names = unique(x3[,7])))+
Heatmap(x3[,8], column_title = "Methylation Sensitivity", column_title_gp = gpar(fontsize = 10), col = c("white","red")) 
dev.off()






mat_data5 = mat_data[ rowSums(mat_data[,c(1,2)])!=0, ] #just a check but all rows are ok
mat_data6 = mat_data[ rowSums(mat_data[,c(3,4)])!=0, ] #just a check but all rows are ok
mat_data7 = mat_data[ rowSums(mat_data[,c(1,3)])!=0, ] #just a check but all rows are ok
mat_data8 = mat_data[ rowSums(mat_data[,c(2,4)])!=0, ] #just a check but all rows are ok

pdf("Promoter_TF_Enrichment_ratio.pdf")
test=heatmap.2(mat_data2,
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
dev.off()


pdf("NEW_Promoter_TF_Enrichment_depletion_sep.pdf")
test=heatmap.2(mat_data5[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data6[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data7[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks)
test=heatmap.2(mat_data8[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks)

dev.off()

################################
## NOT PROMOTER ###
####################

hyper_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_Hyper_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_promoter=hyper_DMRs_promoter$motiftable
hyper_DMRs_promoter=hyper_DMRs_promoter[,1:8]
hyper_DMRs_promoter[,2]=toupper(hyper_DMRs_promoter[,2])
hyper_DMRs_promoter[,8]=as.numeric(hyper_DMRs_promoter[,8])
hyper_DMRs_promoter[,7]=as.numeric(hyper_DMRs_promoter[,7])

hypo_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_Hypo_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_promoter=hypo_DMRs_promoter$motiftable
hypo_DMRs_promoter=hypo_DMRs_promoter[,1:8]
hypo_DMRs_promoter[,2]=toupper(hypo_DMRs_promoter[,2])
hypo_DMRs_promoter[,8]=as.numeric(hypo_DMRs_promoter[,8])
hypo_DMRs_promoter[,7]=as.numeric(hypo_DMRs_promoter[,7])

hypo_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_Hypo_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_promoter=hypo_DAPs_promoter$motiftable
hypo_DAPs_promoter=hypo_DAPs_promoter[,1:8]
hypo_DAPs_promoter[,2]=toupper(hypo_DAPs_promoter[,2])
hypo_DAPs_promoter[,8]=as.numeric(hypo_DAPs_promoter[,8])
hypo_DAPs_promoter[,7]=as.numeric(hypo_DAPs_promoter[,7])

hyper_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/depletion/HAYSTACK_on_Hyper_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_promoter=hyper_DAPs_promoter$motiftable
hyper_DAPs_promoter=hyper_DAPs_promoter[,1:8]
hyper_DAPs_promoter[,2]=toupper(hyper_DAPs_promoter[,2])
hyper_DAPs_promoter[,8]=as.numeric(hyper_DAPs_promoter[,8])
hyper_DAPs_promoter[,7]=as.numeric(hyper_DAPs_promoter[,7])

Haystack_promoter=as.data.frame(cbind(hypo_DMRs_promoter,"HYPO_DMR"))
colnames(Haystack_promoter)[9]="Category"
Haystack_promoter2=as.data.frame(cbind(hyper_DMRs_promoter,"HYPER_DMR"))
colnames(Haystack_promoter2)[9]="Category"

Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hyper_DAPs_promoter,"HYPER_DAP"))
colnames(Haystack_promoter2)[9]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hypo_DAPs_promoter,"HYPO_DAP"))
colnames(Haystack_promoter2)[9]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)
Haystack_promoter[,7]=as.numeric(Haystack_promoter[,7])
Haystack_promoter$sig="NO"
Haystack_promoter$sig[which(Haystack_promoter[,7]<=0.000001)]="YES"

# with 0.000001 cutoff
# NO YES 
# 466 832

Haystack_promoter[,8]=log10(Haystack_promoter[,8])
Haystack_promoter[,8][Haystack_promoter[,7]>0.000001]<-0 


x=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPO_DMR"),]
y=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPER_DAP"),]
z=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPER_DMR"),]
a=Haystack_promoter[,c(2,8)][which(Haystack_promoter$Category=="HYPO_DAP"),]

q=merge(x,y,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPO_DMR","HYPER_DAP"))
r=merge(z,a,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPER_DMR","HYPO_DAP"))

s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)
colnames(s)=c("TF","HYPO_DMR","HYPER_DAP","HYPER_DMR","HYPO_DAP")

rownames(s)=NULL


t=intersect(expressed$gene_symbol,s[,1])
s$Expressed="NO"
s$Expressed[which(s[,1] %in% t)] = "YES"

#Manually change to yes
s$Expressed[c(7,16,21,102,130,201,205,230,262,273,299,350,351,353,355,368,404,406,409,412,438,439,444,445,)]="YES"

s=s[which(s$Expressed=="YES"),]
rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:5])  
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0 
mat_data[is.na(mat_data)] <- 


#mat_data2=mat_data[which(rowSums(mat_data)!=0),]

summary(abs(Haystack_promoter[,8]))
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.31876  0.00000  0.00000  0.05706  0.13354  0.51587 

breaks <- seq(-0.6,0.6, length.out = 40)

library(circlize)
library(RColorBrewer)
cols=colorRampPalette(c("#A50026", "#D73027", "#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))(length(breaks) - 1)


mat_data5 = mat_data[ rowSums(mat_data[,c(1,2)])!=0, ] #just a check but all rows are ok
mat_data6 = mat_data[ rowSums(mat_data[,c(3,4)])!=0, ] #just a check but all rows are ok
mat_data7 = mat_data[ rowSums(mat_data[,c(1,3)])!=0, ] #just a check but all rows are ok
mat_data8 = mat_data[ rowSums(mat_data[,c(2,4)])!=0, ] #just a check but all rows are ok


pdf("NOT_promoter_TF_Enrichment_ratio.pdf")
test=heatmap.2(mat_data,
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
dev.off()


pdf("NEW_NOT_Promoter_TF_Enrichment_depletion_sep.pdf")
test=heatmap.2(mat_data5[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data6[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data7[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks)
test=heatmap.2(mat_data8[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= rev(cols),  
  scale="none", dendrogram="row", 
  breaks=breaks)

dev.off()









#########################################################






test=heatmap.2(mat_data6[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row",
  breaks=breaks) 
test=heatmap.2(mat_data1[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none",dendrogram="row",
  breaks=breaks) 

test=heatmap.2(mat_data2[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none",dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data3[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data4[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
dev.off()













#since already checked for expression in paper......just manually remove rows with no expression for the promoter group at least
AR 
ASCL2
ATOH1
BHLHA15
BHLHE23
EWSR1
FOXA1
FOXA2
HAND1:TCF3
KLF1
MAF:NFE2
MSC
MYF6
MYOD1
NEUROG1 and 2
NFE2
NFIC:TLX1
NKX2-3
OLIG3
TBX15
TBX20
TWIST2
s=s[c(2:3,6:8,10,12:21,23:26,29:33,35:46,48:53,55:64,66:69,72:74,78:82,84:87,89:102,105:110,112:114),]
rownames(s)=NULL
s=s[c(1:62,64:92),]
rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:ncol(s)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0

summary(Haystack_promoter_sig$Ratio)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.060   1.167   1.260   1.381   1.500   2.730

breaks <- seq(1, 3, length.out = 30)
breaks=c(0,breaks)

library(pheatmap)   
library(gplots)

mat_data1 = mat_data[ rowSums(mat_data[,c(1,2)])>=1, ] #just a check but all rows are ok
mat_data2 = mat_data[ rowSums(mat_data[,c(3,4)])>=1, ] #just a check but all rows are ok
mat_data3 = mat_data[ rowSums(mat_data[,c(1,3)])>=1, ] #just a check but all rows are ok
mat_data4 = mat_data[ rowSums(mat_data[,c(2,4)])>=1, ] #just a check but all rows are ok
mat_data5 = mat_data[ rowSums(mat_data[,c(1:4)])>=1, ] #just a check but all rows are ok
mat_data6 = mat_data[ rowSums(mat_data[,c(1:4)])>=2, ] #just a check but all rows are ok

library(viridis)
pdf("Promoter_TF_Enrichment_Ratio_sep.pdf")
test=heatmap.2(mat_data5[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data6[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row",
  breaks=breaks) 
test=heatmap.2(mat_data1[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none",dendrogram="row",
  breaks=breaks) 

test=heatmap.2(mat_data2[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none",dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data3[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data4[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 

dev.off()

##########  Same but not promoter  ###########

hypo_DAPs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hypo_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_not_promoter=hypo_DAPs_not_promoter$motiftable
hypo_DAPs_not_promoter=hypo_DAPs_not_promoter[,1:7]
hypo_DAPs_not_promoter[,2]=toupper(hypo_DAPs_not_promoter[,2])

hypo_DMRs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hypo_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_not_promoter=hypo_DMRs_not_promoter$motiftable
hypo_DMRs_not_promoter=hypo_DMRs_not_promoter[,1:7]
hypo_DMRs_not_promoter[,2]=toupper(hypo_DMRs_not_promoter[,2])

hyper_DAPs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hyper_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_not_promoter=hyper_DAPs_not_promoter$motiftable
hyper_DAPs_not_promoter=hyper_DAPs_not_promoter[,1:7]
hyper_DAPs_not_promoter[,2]=toupper(hyper_DAPs_not_promoter[,2])

hyper_DMRs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hyper_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_not_promoter=hyper_DMRs_not_promoter$motiftable
hyper_DMRs_not_promoter=hyper_DMRs_not_promoter[,1:7]
hyper_DMRs_not_promoter[,2]=toupper(hyper_DMRs_not_promoter[,2])

Haystack_not_promoter=as.data.frame(cbind(hyper_DAPs_not_promoter,"HYPER_DAP"))
colnames(Haystack_not_promoter)[8]="Category"
Haystack_not_promoter2=as.data.frame(cbind(hypo_DAPs_not_promoter,"HYPO_DAP"))
colnames(Haystack_not_promoter2)[8]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)
Haystack_not_promoter2=as.data.frame(cbind(hypo_DMRs_not_promoter,"HYPO_DMR"))
colnames(Haystack_not_promoter2)[8]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)
Haystack_not_promoter2=as.data.frame(cbind(hyper_DMRs_not_promoter,"HYPER_DMR"))
colnames(Haystack_not_promoter2)[8]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)



Haystack_not_promoter[,7]=as.numeric(Haystack_not_promoter[,7])
Haystack_not_promoter_sig=Haystack_not_promoter[which(Haystack_not_promoter[,7]<=0.05),] # none dropped
Haystack_not_promoter[,5]=as.numeric(Haystack_not_promoter[,5])

x=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPO_DMR"),]
y=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPER_DAP"),]
z=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPER_DMR"),]
a=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPO_DAP"),]

q=merge(x,y,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPO_DMR","HYPER_DAP"))
r=merge(z,a,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPER_DMR","HYPO_DAP"))
s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)

colnames(s)=c("TF","HYPO_DMR","HYPER_DAP","HYPER_DMR","HYPO_DAP")
rownames(s)=NULL

t=intersect(expressed$gene_symbol,s[,1])
s$Expressed="NO"
s$Expressed[which(s[,1] %in% t)] = "YES"
#now manually change some to YES

15 20 84 150 154 177 203 257 262 293 296 298 311 317
s$Expressed[c(15,20, 84, 150, 154, 177, 203, 257, 262 ,293, 296, 298, 311, 317)]="YES"

s=s[which(s$Expressed=="YES"),]
rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:5])  
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0

summary(Haystack_not_promoter$Ratio)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.040   1.170   1.250   1.336   1.390   3.460

breaks <- seq(1, 4, length.out = 40)
breaks=c(0,breaks)
mat_data1 = mat_data[ rowSums(mat_data[,c(1,2)])>=1, ] #just a check but all rows are ok
mat_data2 = mat_data[ rowSums(mat_data[,c(3,4)])>=1, ] #just a check but all rows are ok
mat_data3 = mat_data[ rowSums(mat_data[,c(1,3)])>=1, ] #just a check but all rows are ok
mat_data4 = mat_data[ rowSums(mat_data[,c(2,4)])>=1, ] #just a check but all rows are ok
mat_data5 = mat_data[ rowSums(mat_data[,c(1:4)])>=1, ] #just a check but all rows are ok
mat_data6 = mat_data[ rowSums(mat_data[,c(1:4)])>=2, ] #just a check but all rows are ok

pdf("Not_Promoter_TF_Enrichment_Ratio_sep.pdf")
test=heatmap.2(mat_data5[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data6[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row",
  breaks=breaks) 
test=heatmap.2(mat_data1[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none",dendrogram="row",
  breaks=breaks) 

test=heatmap.2(mat_data2[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none",dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data3[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data4[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
dev.off()







#################################################################################

## THESE ARE DARS THAT OVERLAP HYPER OR HYPO DMRS in PROMOTERS
hyper_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/HAYSTACK_on_hyperMeth_dars_in_dmrs_promoter_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_promoter=hyper_DMRs_promoter$motiftable
hyper_DMRs_promoter=hyper_DMRs_promoter[,1:8]
hyper_DMRs_promoter[,2]=toupper(hyper_DMRs_promoter[,2])
hyper_DMRs_promoter[,5]=as.numeric(hyper_DMRs_promoter[,5])
hyper_DMRs_promoter[,6]=as.numeric(hyper_DMRs_promoter[,6])
hyper_DMRs_promoter[,7]=as.numeric(hyper_DMRs_promoter[,7])
hyper_DMRs_promoter[,8]=as.numeric(hyper_DMRs_promoter[,8])

hypo_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/HAYSTACK_on_hypoMeth_dars_in_dmrs_promoter_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_promoter=hypo_DMRs_promoter$motiftable
hypo_DMRs_promoter=hypo_DMRs_promoter[,1:8]
hypo_DMRs_promoter[,2]=toupper(hypo_DMRs_promoter[,2])
hypo_DMRs_promoter[,5]=as.numeric(hypo_DMRs_promoter[,5])
hypo_DMRs_promoter[,6]=as.numeric(hypo_DMRs_promoter[,6])
hypo_DMRs_promoter[,7]=as.numeric(hypo_DMRs_promoter[,7])
hypo_DMRs_promoter[,8]=as.numeric(hypo_DMRs_promoter[,8])





hypo_DMRs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hypo_DMRs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_promoter=hypo_DMRs_promoter$motiftable
hypo_DMRs_promoter=hypo_DMRs_promoter[,1:7]
hypo_DMRs_promoter[,2]=toupper(hypo_DMRs_promoter[,2])
hypo_DMRs_promoter[,5:8]=as.numeric(hypo_DMRs_promoter[,5:8])

hypo_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hypo_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_promoter=hypo_DAPs_promoter$motiftable
hypo_DAPs_promoter=hypo_DAPs_promoter[,1:7]
hypo_DAPs_promoter[,2]=toupper(hypo_DAPs_promoter[,2])
hypo_DAPs_promoter$Ratio=as.numeric(hypo_DAPs_promoter$Ratio)

hyper_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hyper_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_promoter=hyper_DAPs_promoter$motiftable
hyper_DAPs_promoter=hyper_DAPs_promoter[,1:7]
hyper_DAPs_promoter[,2]=toupper(hyper_DAPs_promoter[,2])
hyper_DAPs_promoter$Ratio=as.numeric(hyper_DAPs_promoter$Ratio)

hyper_DMRs_hypo_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hyper_DMRs_hypo_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_hypo_DAPs_promoter=hyper_DMRs_hypo_DAPs_promoter$motiftable
hyper_DMRs_hypo_DAPs_promoter=hyper_DMRs_hypo_DAPs_promoter[,1:7]
hyper_DMRs_hypo_DAPs_promoter[,2]=toupper(hyper_DMRs_hypo_DAPs_promoter[,2])
hyper_DMRs_hypo_DAPs_promoter$Ratio=as.numeric(hyper_DMRs_hypo_DAPs_promoter$Ratio)

hypo_DMRs_hyper_DAPs_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_hypo_DMRs_hyper_DAPs_overlap_promoter_NAvBA9pos_VS_promoters_by_gene/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_hyper_DAPs_promoter=hypo_DMRs_hyper_DAPs_promoter$motiftable
hypo_DMRs_hyper_DAPs_promoter=hypo_DMRs_hyper_DAPs_promoter[,1:7]
hypo_DMRs_hyper_DAPs_promoter[,2]=toupper(hypo_DMRs_hyper_DAPs_promoter[,2])
hypo_DMRs_hyper_DAPs_promoter$Ratio=as.numeric(hypo_DMRs_hyper_DAPs_promoter$Ratio)




Haystack_promoter=as.data.frame(cbind(hypo_DMRs_promoter,"HYPO_DMR"))
colnames(Haystack_promoter)[8]="Category"
Haystack_promoter2=as.data.frame(cbind(hyper_DMRs_promoter,"HYPER_DMR"))
colnames(Haystack_promoter2)[8]="Category"

Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hyper_DAPs_promoter,"HYPER_DAP"))
colnames(Haystack_promoter2)[8]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hypo_DAPs_promoter,"HYPO_DAP"))
colnames(Haystack_promoter2)[8]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hypo_DMRs_hyper_DAPs_promoter,"HYPO_DMR_union"))
colnames(Haystack_promoter2)[8]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter2=as.data.frame(cbind(hyper_DMRs_hypo_DAPs_promoter,"HYPER_DMR_union"))
colnames(Haystack_promoter2)[8]="Category"
Haystack_promoter=rbind(Haystack_promoter,Haystack_promoter2)

Haystack_promoter[,7]=as.numeric(Haystack_promoter[,7])
Haystack_promoter_sig=Haystack_promoter[which(Haystack_promoter[,7]<=0.05),] # only dropped 5 rows

x=Haystack_promoter_sig[,c(2,5)][which(Haystack_promoter_sig$Category=="HYPO_DMR"),]
y=Haystack_promoter_sig[,c(2,5)][which(Haystack_promoter_sig$Category=="HYPER_DAP"),]
z=Haystack_promoter_sig[,c(2,5)][which(Haystack_promoter_sig$Category=="HYPER_DMR"),]
a=Haystack_promoter_sig[,c(2,5)][which(Haystack_promoter_sig$Category=="HYPO_DAP"),]
b=Haystack_promoter_sig[,c(2,5)][which(Haystack_promoter_sig$Category=="HYPO_DMR_union"),]
c=Haystack_promoter_sig[,c(2,5)][which(Haystack_promoter_sig$Category=="HYPER_DMR_union"),]



q=merge(x,y,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPO_DMR","HYPER_DAP"))
r=merge(z,a,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPER_DMR","HYPO_DAP"))
t=merge(b,c,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPO_DMR_union","HYPER_DMR_union"))

p=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)
s=merge(p,t,by="Motif Name",all.x=TRUE,all.y=TRUE)
colnames(s)=c("TF","HYPO_DMR","HYPER_DAP","HYPER_DMR","HYPO_DAP","HYPO_DMR_union","HYPER_DMR_union")

rownames(s)=NULL

#since already checked for expression in paper......just manually remove rows with no expression for the promoter group at least
AR 
ASCL2
ATOH1
BHLHA15
BHLHE23
EWSR1
FOXA1
FOXA2
HAND1:TCF3
KLF1
MAF:NFE2
MSC
MYF6
MYOD1
NEUROG1 and 2
NFE2
NFIC:TLX1
NKX2-3
OLIG3
TBX15
TBX20
TWIST2
s=s[c(2:3,6:8,10,12:21,23:26,29:33,35:46,48:53,55:64,66:69,72:74,78:82,84:87,89:102,105:110,112:114),]
rownames(s)=NULL
s=s[c(1:62,64:92),]
rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:ncol(s)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0

summary(Haystack_promoter_sig$Ratio)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.060   1.167   1.260   1.381   1.500   2.730

breaks <- seq(1, 3, length.out = 30)
breaks=c(0,breaks)

library(pheatmap)   
library(gplots)

mat_data1 = mat_data[ rowSums(mat_data[,c(1,2)])>=1, ] #just a check but all rows are ok
mat_data2 = mat_data[ rowSums(mat_data[,c(3,4)])>=1, ] #just a check but all rows are ok
mat_data3 = mat_data[ rowSums(mat_data[,c(1,3)])>=1, ] #just a check but all rows are ok
mat_data4 = mat_data[ rowSums(mat_data[,c(2,4)])>=1, ] #just a check but all rows are ok
mat_data5 = mat_data[ rowSums(mat_data[,c(1:4)])>=1, ] #just a check but all rows are ok
mat_data6 = mat_data[ rowSums(mat_data[,c(1:4)])>=2, ] #just a check but all rows are ok

library(viridis)
pdf("Promoter_TF_Enrichment_Ratio_sep.pdf")
test=heatmap.2(mat_data5[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data6[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row",
  breaks=breaks) 
test=heatmap.2(mat_data1[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none",dendrogram="row",
  breaks=breaks) 

test=heatmap.2(mat_data2[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none",dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data3[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data4[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.5, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 29)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 

dev.off()

##########  Same but not promoter  ###########

hypo_DAPs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hypo_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DAPs_not_promoter=hypo_DAPs_not_promoter$motiftable
hypo_DAPs_not_promoter=hypo_DAPs_not_promoter[,1:7]
hypo_DAPs_not_promoter[,2]=toupper(hypo_DAPs_not_promoter[,2])

hypo_DMRs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hypo_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hypo_DMRs_not_promoter=hypo_DMRs_not_promoter$motiftable
hypo_DMRs_not_promoter=hypo_DMRs_not_promoter[,1:7]
hypo_DMRs_not_promoter[,2]=toupper(hypo_DMRs_not_promoter[,2])

hyper_DAPs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hyper_DAPs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DAPs_not_promoter=hyper_DAPs_not_promoter$motiftable
hyper_DAPs_not_promoter=hyper_DAPs_not_promoter[,1:7]
hyper_DAPs_not_promoter[,2]=toupper(hyper_DAPs_not_promoter[,2])

hyper_DMRs_not_promoter=readHTMLTable("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Haystack_analyses/Promoter_Analyses/HAYSTACK_on_Hyper_DMRs_NAvBA9pos_not_in_promoter_VS_not_promoters/Haystack_report.html",header=TRUE,skip.rows=4,stringsAsFactors=FALSE)
hyper_DMRs_not_promoter=hyper_DMRs_not_promoter$motiftable
hyper_DMRs_not_promoter=hyper_DMRs_not_promoter[,1:7]
hyper_DMRs_not_promoter[,2]=toupper(hyper_DMRs_not_promoter[,2])

Haystack_not_promoter=as.data.frame(cbind(hyper_DAPs_not_promoter,"HYPER_DAP"))
colnames(Haystack_not_promoter)[8]="Category"
Haystack_not_promoter2=as.data.frame(cbind(hypo_DAPs_not_promoter,"HYPO_DAP"))
colnames(Haystack_not_promoter2)[8]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)
Haystack_not_promoter2=as.data.frame(cbind(hypo_DMRs_not_promoter,"HYPO_DMR"))
colnames(Haystack_not_promoter2)[8]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)
Haystack_not_promoter2=as.data.frame(cbind(hyper_DMRs_not_promoter,"HYPER_DMR"))
colnames(Haystack_not_promoter2)[8]="Category"
Haystack_not_promoter=rbind(Haystack_not_promoter,Haystack_not_promoter2)



Haystack_not_promoter[,7]=as.numeric(Haystack_not_promoter[,7])
Haystack_not_promoter_sig=Haystack_not_promoter[which(Haystack_not_promoter[,7]<=0.05),] # none dropped
Haystack_not_promoter[,5]=as.numeric(Haystack_not_promoter[,5])

x=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPO_DMR"),]
y=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPER_DAP"),]
z=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPER_DMR"),]
a=Haystack_not_promoter[,c(2,5)][which(Haystack_not_promoter$Category=="HYPO_DAP"),]

q=merge(x,y,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPO_DMR","HYPER_DAP"))
r=merge(z,a,by="Motif Name",all.x=TRUE,all.y=TRUE,suffixes=c("HYPER_DMR","HYPO_DAP"))
s=merge(q,r,by="Motif Name",all.x=TRUE,all.y=TRUE)

colnames(s)=c("TF","HYPO_DMR","HYPER_DAP","HYPER_DMR","HYPO_DAP")
rownames(s)=NULL

t=intersect(expressed$gene_symbol,s[,1])
s$Expressed="NO"
s$Expressed[which(s[,1] %in% t)] = "YES"
#now manually change some to YES

15 20 84 150 154 177 203 257 262 293 296 298 311 317
s$Expressed[c(15,20, 84, 150, 154, 177, 203, 257, 262 ,293, 296, 298, 311, 317)]="YES"

s=s[which(s$Expressed=="YES"),]
rnames <- s[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(s[,2:5])  
rownames(mat_data) <- rnames 
mat_data[is.na(mat_data)] <- 0

summary(Haystack_not_promoter$Ratio)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.040   1.170   1.250   1.336   1.390   3.460

breaks <- seq(1, 4, length.out = 40)
breaks=c(0,breaks)
mat_data1 = mat_data[ rowSums(mat_data[,c(1,2)])>=1, ] #just a check but all rows are ok
mat_data2 = mat_data[ rowSums(mat_data[,c(3,4)])>=1, ] #just a check but all rows are ok
mat_data3 = mat_data[ rowSums(mat_data[,c(1,3)])>=1, ] #just a check but all rows are ok
mat_data4 = mat_data[ rowSums(mat_data[,c(2,4)])>=1, ] #just a check but all rows are ok
mat_data5 = mat_data[ rowSums(mat_data[,c(1:4)])>=1, ] #just a check but all rows are ok
mat_data6 = mat_data[ rowSums(mat_data[,c(1:4)])>=2, ] #just a check but all rows are ok

pdf("Not_Promoter_TF_Enrichment_Ratio_sep.pdf")
test=heatmap.2(mat_data5[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data6[,c(1:4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row",
  breaks=breaks) 
test=heatmap.2(mat_data1[,c(1,2)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none",dendrogram="row",
  breaks=breaks) 

test=heatmap.2(mat_data2[,c(3,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none",dendrogram="row", 
  breaks=breaks) 

test=heatmap.2(mat_data3[,c(1,3)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
test=heatmap.2(mat_data4[,c(2,4)],
  main = "Promoter",
  lhei = c(2, 8),
  trace="none", 
  cexRow=.25, cexCol=1, 
  margins =c(9,25),    
  col= c("lightgrey",plasma(n = 39)),  
  scale="none", dendrogram="row", 
  breaks=breaks) 
dev.off()

