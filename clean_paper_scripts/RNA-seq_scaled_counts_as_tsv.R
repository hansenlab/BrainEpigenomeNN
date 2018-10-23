# Export RNA-seq gene-level scaled counts matrix as tab-separated file
# Peter Hickey
# 2017-02-28

txi_gene <- readRDS("../RNA-seq/objects/txi-gene.flow-sorted-brain-rna-seq.rds")

txi_gene_df <- cbind(data.frame(GENE_ID = rownames(txi_gene$counts),
                                stringsAsFactors = FALSE),
                     as.data.frame(txi_gene$counts, row.names = NA))

write.table(x = txi_gene_df,
            file = "flow_sorted_brain.RNA-seq_gene_level_scaled_counts.txt.gz",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
