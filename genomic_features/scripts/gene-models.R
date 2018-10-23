# Make gene models and 'flattened' genomic features from GENCODE v19
# annotations
# Peter Hickey
# 2016-10-26

library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)

extdir <- "../extdata"

# NOTE: makeTxDbFromGFF() only imports the "type", "gene_id", "transcript_id",
#       and "exon_id" columns, which is insufficient for our needs. We have to
#       therefore use import() + makeTxDbFromGRanges()
#       (see https://support.bioconductor.org/p/73158/#73288)

###-----------------------------------------------------------------------------
### Import FASTA files (protein-coding transcripts and lncRNAs tranxcripts)
###

fasta_pc <- import(file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                             "GENCODE_v19", "gencode.v19.pc_transcripts.fa.gz"),
                   format = "FASTA")
# > length(fasta_pc)
# [1] 95309

fasta_pc_tx_name <- sapply(strsplit(names(fasta_pc), "\\|"), "[[", 1)
fasta_pc_gene_name <- sapply(strsplit(names(fasta_pc), "\\|"), "[[", 2)

fasta_lnc <- import(file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                              "GENCODE_v19",
                              "gencode.v19.lncRNA_transcripts.fa.gz"),
                    format = "FASTA")
# > length(fasta_lnc)
# [1] 23898

fasta_lnc_tx_name <- sapply(strsplit(names(fasta_lnc), "\\|"), "[[", 1)
fasta_lnc_gene_name <- sapply(strsplit(names(fasta_lnc), "\\|"), "[[", 2)

###-----------------------------------------------------------------------------
### Import GTF file as GRanges objects, only retaining those transcripts in the
### FASTA files
###

gtf <- import(file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                        "GENCODE_v19", "gencode.v19.annotation.gtf.gz"))
# > length(gtf)
# [1] 2619444

fasta_id <- ifelse(gtf$transcript_id %in% sub("\\|.*", "", names(fasta_pc)),
                   "pc_transcripts",
                   ifelse(gtf$transcript_id %in% sub("\\|.*", "",
                                                     names(fasta_lnc)),
                          "lnc_transcripts",
                          NA_character_))
names(fasta_id) <- gtf$transcript_id
gtf$fasta_id <- fasta_id
gtf <- gtf[!is.na(gtf$fasta_id)]
# > length(gtf)
# [1] 2213418

# > any(!sub("\\|.*", "", names(fasta_pc)) %in%
#         +       gtf[gtf$type == "transcript"]$transcript_id)
#
# [1] FALSE
# > any(!sub("\\|.*", "", names(fasta_lnc)) %in%
#         +       gtf[gtf$type == "transcript"]$transcript_id)
# [1] FALSE

###-----------------------------------------------------------------------------
### Add SeqInfo and drop non-autosomal seqlevels and
###

seqinfo(gtf) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
gtf <- keepSeqlevels(gtf, paste0("chr", 1:22))
# > length(gtf)
# [1] 2139571
#
# > table(gtf$type)
#
# gene     transcript           exon            CDS    start_codon
#    0         115461         895868         698286          81156
# stop_codon            UTR Selenocysteine
#      73441         275245            114
#
# > length(unique(gtf$transcript_id))
# [1] 115461
#
# > length(unique(gtf$gene_id))
# [1] 33351

# NOTE: fasta_pc and fasta_lnc may both contain non-autosomal transcripts. We
#       don't want to remove these because it would affect the quasi-mapping.
#       However, we need to bear this in mind when looking for transcripts
#       that are 'missing' GTF annotations (since we only retain autosomal
#       GTF annotations).

###-----------------------------------------------------------------------------
### Convert to TxDb object
###

txdb <- makeTxDbFromGRanges(gr = gtf,
                            drop.stop.codons = FALSE,
                            metadata = NULL,
                            taxonomyId = 9606)
# > txdb
# TxDb object:
#   # Db type: TxDb
#   # Supporting package: GenomicFeatures
#   # Genome: hg19
#   # Organism: Homo sapiens
#   # Taxonomy ID: 9606
#   # transcript_nrow: 115461
#   # exon_nrow: 461387
#   # cds_nrow: 259809
#   # Db created by: GenomicFeatures package from Bioconductor
#   # Creation time: 2016-10-11 18:27:40 -0400 (Tue, 11 Oct 2016)
#   # GenomicFeatures version at creation time: 1.24.5
# # RSQLite version at creation time: 1.0.0
# # DBSCHEMAVERSION: 1.1

###-----------------------------------------------------------------------------
### Extract features of interest
###

transcripts <- transcripts(txdb, columns = columns(txdb))
# > length(transcripts)
# [1] 115461
genes <- genes(txdb, columns = columns(txdb))
# > length(genes)
# [1] 33351
promoters <- promoters(txdb,
                       columns = columns(txdb),
                       upstream = 2000,
                       downstream = 2000)
# > length(promoters)
# [1] 115461
transcripts_by_gene <- transcriptsBy(txdb, by = "gene")
# > length(transcripts_by_gene)
# [1] 33351
five_utrs_by_transcript <- fiveUTRsByTranscript(txdb, use.names = TRUE)
# > length(five_utrs_by_transcript)
# [1] 77532
three_utrs_by_transcript <- threeUTRsByTranscript(txdb, use.names = TRUE)
# > length(three_utrs_by_transcript)
# [1] 69726
exons_by_transcript <- exonsBy(txdb, by = "tx", use.names = TRUE)
# > length(exons_by_transcript)
# [1] 115461
introns_by_transcript <- intronsByTranscript(txdb, use.names = TRUE)
# > length(introns_by_transcript)
# [1] 115461

###-----------------------------------------------------------------------------
### 'Unflattened' features
###

unflattened_features <- list(transcripts = transcripts,
                             genes = genes,
                             promoters = promoters,
                             transcripts_by_gene = transcripts_by_gene,
                             five_utrs_by_transcript = five_utrs_by_transcript,
                             three_utrs_by_transcript = three_utrs_by_transcript,
                             exons_by_transcript = exons_by_transcript,
                             introns_by_transcript = introns_by_transcript)

###-----------------------------------------------------------------------------
### 'Unflattened features' based on protein-coding transcript sequences
###

unflattened_features_pc_transcripts <-
  list(transcripts = transcripts[transcripts$TXNAME %in% fasta_pc_tx_name],
       genes = genes[unlist(genes$GENEID) %in% fasta_pc_gene_name],
       promoters = promoters[promoters$TXNAME %in% fasta_pc_tx_name],
       transcripts_by_gene = transcripts_by_gene[
         na.omit(match(fasta_pc_gene_name,
                       names(transcripts_by_gene)))],
       five_utrs_by_transcript = five_utrs_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(five_utrs_by_transcript)))],
       three_utrs_by_transcript = three_utrs_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(three_utrs_by_transcript)))],
       exons_by_transcript = exons_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(exons_by_transcript)))],
       introns_by_transcript = introns_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(introns_by_transcript)))])

###-----------------------------------------------------------------------------
### 'Unflattened features' based on lncRNA transcript sequences
###

unflattened_features_lnc_transcripts <-
  list(transcripts = transcripts[transcripts$TXNAME %in% fasta_lnc_tx_name],
       genes = genes[unlist(genes$GENEID) %in% fasta_lnc_gene_name],
       promoters = promoters[promoters$TXNAME %in% fasta_lnc_tx_name],
       transcripts_by_gene = transcripts_by_gene[
         na.omit(match(fasta_lnc_gene_name,
                       names(transcripts_by_gene)))],
       five_utrs_by_transcript = five_utrs_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(five_utrs_by_transcript)))],
       three_utrs_by_transcript = three_utrs_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(three_utrs_by_transcript)))],
       exons_by_transcript = exons_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(exons_by_transcript)))],
       introns_by_transcript = introns_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(introns_by_transcript)))])

###-----------------------------------------------------------------------------
### Sanity checking UTRs

# Check that no lncRNA has 5' or 3' UTRs
# > any(names(five_utrs_by_transcript) %in% sub("\\|.*", "", names(fasta_lnc)))
# [1] FALSE
#
# > any(names(three_utrs_by_transcript) %in% sub("\\|.*", "", names(fasta_lnc)))
# [1] FALSE

###-----------------------------------------------------------------------------
### Find the transcripts that lack 5' UTRs
###

no_five_utr <- transcripts[!transcripts$TXNAME %in%
                             names(five_utrs_by_transcript)]

# > table(gtf$fasta_id, gtf$transcript_id %in% no_five_utr$TXNAME)
#
#                   FALSE    TRUE
# lnc_transcripts       0   95235
# pc_transcripts  1786146  258190


# > table(gtf$transcript_type, gtf$transcript_id %in% no_five_utr$TXNAME)
#
#                            FALSE    TRUE
# 3prime_overlapping_ncrna       0      75
# antisense                      0   38089
# IG_C_gene                      0     171
# IG_D_gene                      0      85
# IG_J_gene                     14      50
# IG_V_gene                    861     124
# lincRNA                        0   46186
# miRNA                          0     118
# misc_RNA                       0      32
# nonsense_mediated_decay   212511   64510
# non_stop_decay               860     210
# processed_transcript           0    4843
# protein_coding           1571366  192623
# retained_intron                0    2204
# rRNA                           0       8
# sense_intronic                 0    2354
# sense_overlapping              0    1168
# snoRNA                         0     144
# snRNA                          0      14
# TR_C_gene                      0      53
# TR_D_gene                      0       9
# TR_J_gene                      0     223
# TR_V_gene                    534     132

# > table(is.na(gtf$ccdsid), gtf$transcript_id %in% no_five_utr$TXNAME)
#
#        FALSE   TRUE
# FALSE 872602  49147
# TRUE  913544 304278

###-----------------------------------------------------------------------------
### Find the transcripts that lack 3' UTRs
###

no_three_utr <- transcripts[!transcripts$TXNAME %in%
                              names(three_utrs_by_transcript)]

# > table(gtf$fasta_id, gtf$transcript_id %in% no_three_utr$TXNAME)
#
#                   FALSE    TRUE
# lnc_transcripts       0   95235
# pc_transcripts  1718169  326167

# > table(gtf$transcript_type, gtf$transcript_id %in% no_three_utr$TXNAME)
#
#                            FALSE    TRUE
# 3prime_overlapping_ncrna       0      75
# antisense                      0   38089
# IG_C_gene                    171       0
# IG_D_gene                     10      75
# IG_J_gene                      0      64
# IG_V_gene                     29     956
# lincRNA                        0   46186
# miRNA                          0     118
# misc_RNA                       0      32
# nonsense_mediated_decay   276938      83
# non_stop_decay                 0    1070
# processed_transcript           0    4843
# protein_coding           1440968  323021
# retained_intron                0    2204
# rRNA                           0       8
# sense_intronic                 0    2354
# sense_overlapping              0    1168
# snoRNA                         0     144
# snRNA                          0      14
# TR_C_gene                     53       0
# TR_D_gene                      0       9
# TR_J_gene                      0     223
# TR_V_gene                      0     666

# > table(is.na(gtf$ccdsid), gtf$transcript_id %in% no_three_utr$TXNAME)
#
#        FALSE   TRUE
# FALSE 885160  36589
# TRUE  833009 384813


###-----------------------------------------------------------------------------
### Find the transcripts that have 5' UTRs but lack 3' UTRs
###

have_five_utr_no_three_utr <- transcripts[transcripts$TXNAME %in%
                                            names(five_utrs_by_transcript) &
                                            !transcripts$TXNAME %in%
                                            names(three_utrs_by_transcript)]

# > length(have_five_utr_no_three_utr)
# [1] 18314

# > table(gtf$fasta_id, gtf$transcript_id %in% have_five_utr_no_three_utr$TXNAME)
#
#                   FALSE    TRUE
# lnc_transcripts   95235       0
# pc_transcripts  1783293  261043

# > table(gtf$transcript_type, gtf$transcript_id %in% have_five_utr_no_three_utr$TXNAME)
#
#                            FALSE    TRUE
# 3prime_overlapping_ncrna      75       0
# antisense                  38089       0
# IG_C_gene                    171       0
# IG_D_gene                     85       0
# IG_J_gene                     50      14
# IG_V_gene                    153     832
# lincRNA                    46186       0
# miRNA                        118       0
# misc_RNA                      32       0
# nonsense_mediated_decay   276950      71
# non_stop_decay               210     860
# processed_transcript        4843       0
# protein_coding           1505257  258732
# retained_intron             2204       0
# rRNA                           8       0
# sense_intronic              2354       0
# sense_overlapping           1168       0
# snoRNA                       144       0
# snRNA                         14       0
# TR_C_gene                     53       0
# TR_D_gene                      9       0
# TR_J_gene                    223       0
# TR_V_gene                    132     534

###-----------------------------------------------------------------------------
### 'Flatten' features by unlisting, unstranding, and reducing
### NOTE: This procedure means that a base may belong to multiple categories
###

# NOTE: Intergenic regions are defined as those without a gene on either strand
intergenic <- gaps(reduce(unstrand(genes)))
intergenic <- intergenic[strand(intergenic) == "*"]

flattened_features <- list(genic = reduce(unstrand(genes)),
                           promoter = reduce(unstrand(promoters)),
                           five_utr =
                             reduce(unstrand(unlist(five_utrs_by_transcript))),
                           three_utr =
                             reduce(unstrand(unlist(three_utrs_by_transcript))),
                           exonic =
                             reduce(unstrand(unlist(exons_by_transcript))),
                           intronic =
                             reduce(unstrand(unlist(introns_by_transcript))),
                           intergenic = intergenic)

#-------------------------------------------------------------------------------
# 'Flattened features' based on protein-coding transcript sequences
#

intergenic_pc_transcripts <-
  gaps(reduce(unstrand(genes[names(genes) %in% fasta_pc_gene_name])))
intergenic_pc_transcripts <- intergenic_pc_transcripts[
  strand(intergenic_pc_transcripts) == "*"]

flattened_features_pc_transcripts <-
  list(genic = reduce(unstrand(genes[names(genes) %in% fasta_pc_gene_name])),
       promoter = reduce(unstrand(promoters[promoters$TXNAME %in%
                                              fasta_pc_tx_name])),
       five_utr =
         reduce(unstrand(unlist(five_utrs_by_transcript[
           names(five_utrs_by_transcript) %in% fasta_pc_tx_name]))),
       three_utr =
         reduce(unstrand(unlist(three_utrs_by_transcript[
           names(three_utrs_by_transcript) %in% fasta_pc_tx_name]))),
       exonic =
         reduce(unstrand(unlist(exons_by_transcript[
           names(exons_by_transcript) %in% fasta_pc_tx_name]))),
       intronic =
         reduce(unstrand(unlist(introns_by_transcript[
           names(introns_by_transcript) %in% fasta_pc_tx_name]))),
       intergenic = intergenic_pc_transcripts)

#-------------------------------------------------------------------------------
# 'Flattened features' based on long non-coding RNA transcript sequences
#

intergenic_lnc_transcripts <-
  gaps(reduce(unstrand(genes[names(genes) %in% fasta_lnc_gene_name])))
intergenic_lnc_transcripts <- intergenic_lnc_transcripts[
  strand(intergenic_lnc_transcripts) == "*"]

flattened_features_lnc_transcripts <-
  list(genic = reduce(unstrand(genes[names(genes) %in% fasta_lnc_gene_name])),
       promoter = reduce(unstrand(promoters[promoters$TXNAME %in%
                                              fasta_lnc_tx_name])),
       five_utr =
         reduce(unstrand(unlist(five_utrs_by_transcript[
           names(five_utrs_by_transcript) %in% fasta_lnc_tx_name]))),
       three_utr =
         reduce(unstrand(unlist(three_utrs_by_transcript[
           names(three_utrs_by_transcript) %in% fasta_lnc_tx_name]))),
       exonic =
         reduce(unstrand(unlist(exons_by_transcript[
           names(exons_by_transcript) %in% fasta_lnc_tx_name]))),
       intronic =
         reduce(unstrand(unlist(introns_by_transcript[
           names(introns_by_transcript) %in% fasta_lnc_tx_name]))),
       intergenic = intergenic_lnc_transcripts)

#-------------------------------------------------------------------------------
# Save unflattened features objects
#

save(unflattened_features,
     unflattened_features_pc_transcripts,
     unflattened_features_lnc_transcripts,
     file = "../objects/unflattened-GENCODE-v19-features.rda",
     compress = "xz")

#-------------------------------------------------------------------------------
# Save flattened features objects
#

save(flattened_features,
     flattened_features_lnc_transcripts,
     flattened_features_pc_transcripts,
     file = "../objects/flattened-GENCODE-v19-features.rda")

#-------------------------------------------------------------------------------
# Make object of genes/transcripts for use in plotting routines
#

saveRDS(transcripts_by_gene,
        file = "../objects/GENCODE-v19-transcripts-by-gene.rds")
saveRDS(sort(granges(genes)),
        file = "../objects/GENCODE-v19-genes.rds")

###-----------------------------------------------------------------------------
### NOTES/QUESTIONS
###

# 1. length(five_utrs_by_transcript) != length(three_utrs_by_transcript)
#   - A transcript can have a 5' UTR and no 3' UTR (see 2 below)
# 2. Don't really understand why some transcripts have no 5' and/or 3' UTRs
#   - E.g., A transcript with transcript_type == "non_stop_decay" will not have
#     a 3' UTR (https://www.gencodegenes.org/gencode_biotypes.html)
# 3. How do fiveUTRsByTranscript() and threeUTRsByTranscript() work?
#   - Each first calls GenomicFeatures:::.getSplicingsForTranscriptsWithCDSs()
#     As the name suggests, this extracts the splicings **for transcripts with
#     CDSs** and not all transcripts have CDSs. Each then extracts the
#     requested UTRs

