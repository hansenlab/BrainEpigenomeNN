# Import Salmon output files from flow sorted brain GTEx project to create
# object(s) for DE analysis
# Peter Hickey
# 2016-09-12

# NOTE: Running on Bioc 3.4 (devel)

library(tximport)
library(S4Vectors)
library(readr)
library(GenomicFeatures)
library(limma)
library(edgeR)
library(DESeq2)

extdir <- "../extdata"

#-------------------------------------------------------------------------------
# Load data

cd <- readRDS("../objects/colData-flow-sorted-brain-rna-seq.rds")
quant_sf <- file.path(cd$SALMON_DIR, "quant.sf")
stopifnot(all(file.exists(quant_sf)))
names(quant_sf) <- rownames(cd)

#-------------------------------------------------------------------------------
# Make TxDb object from GENCODE GTF

txdb <- makeTxDbFromGFF(file = file.path(extdir, "flow-sorted-brain-rna-seq",
                                         "data", "GENCODE_v19",
                                         "gencode.v19.annotation.gtf.gz"),
                        format = "auto",
                        organism = "Homo sapiens",
                        dataSource = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz")

#-------------------------------------------------------------------------------
# Import Salmon output and coerce to:
#
# 1. Transcript-level list
# 2. Gene-level  list
#
# for use with limma-voom

# Import quant.sf files
txi_transcript <- tximport(files = quant_sf,
                           type = "salmon",
                           txOut = TRUE,
                           reader = read_tsv)

# Make tx2gene data frame
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

# Summarize transcripts at the gene level
# NOTE: tximport vignette recommended `countsFromAbundance = "lengthScaledTPM`
#       for downstream use with limma-voom
txi_gene <- summarizeToGene(txi = txi_transcript,
                            tx2gene = tx2gene,
                            countsFromAbundance = "lengthScaledTPM")

# Save txi objects
saveRDS(txi_transcript,
        "../objects/txi-transcript.flow-sorted-brain-rna-seq.rds")
saveRDS(txi_gene, "../objects/txi-gene.flow-sorted-brain-rna-seq.rds")

#-------------------------------------------------------------------------------
# Get all TSSs for each GENCODE gene

tss_grl <- transcriptsBy(txdb, by = "gene")
saveRDS(tss_grl, "../objects/GENCODE_v19-TSSs_grl.rds")

