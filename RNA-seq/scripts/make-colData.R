# Make colData object for flow sorted brain GTEx project bulk RNA-seq data to
# store metadata of these samples
# Peter Hickey
# 2016-09-10

library(S4Vectors)
library(readr)
library(dplyr)
library(stringr)

extdir <- "../extdata"

#-------------------------------------------------------------------------------
# Read in sample sheet
sample_sheet <- read_csv(file.path(extdir, "RNA-seq_160905_FlowSortedPaper.csv"),
                         skip = 1)
# Sanitise colnames
colnames(sample_sheet) <- c("LIBRARY",
                            "LIBRARY_CONSTRUCTED_BY",
                            "SAMPLE",
                            "AVERAGE_LENGTH_LIBRARY_IN_BP",
                            "qPCR_IN_nM_FROM_BioA_TRACE",
                            "SIZING_LIBRARY_QC_ASSAY",
                            "qPCR_ADJUSTMENT_FACTOR",
                            "POST_qPCR_IN_nM",
                            "ADDITIONAL_LIBRARY_QUANTIFICATION",
                            "PRIMER_DETAILS_REMARKS",
                            "INDEX_SEQUENCE_GENERATED_DURING_SEQUENCING",
                            "LIBRARY_TYPE",
                            "LIBRARY_KIT",
                            "CT_CONVERSION_KIT",
                            "HISEQ_OPERATOR",
                            "HISEQ_CHEMISTRY",
                            "HISEQ_NUMBER",
                            "FC_LANES",
                            "TOTAL_LIB_CONC_IN_pM_CLUSTERED_PER_LANE",
                            "PhiX_SPIKE_IN_PER_LANE",
                            "RUN_PARAMETERS",
                            "qPCR_CLUSTERING_BY",
                            "CLUSTERING_DATE",
                            "ASSAY_TO_QUANT_INPUT_DNA_IN_ng",
                            "INPUT_DNA_AMOUNT_IN_ng_FOR_NGS_LIB_PREP",
                            "PCR_CYCLES_USED_DURING_LIBPREP",
                            "ELUTION_VOLUME_OF_FINAL_LIBRARY_IN_ul",
                            "LIBRARY_CONSTRUCTION_COMPLETION_DATE",
                            "NOTES",
                            "LIBRARY_STORED_IN_PLATE",
                            "GRID_REFERENCE",
                            "NOT_USED_1",
                            "NOT_USED_2",
                            "NOT_USED_3",
                            "NOT_USED_4",
                            "NOT_USED_5")

# Drop unused columns
sample_sheet <- sample_sheet %>%
  select(-NOT_USED_1, -NOT_USED_2, -NOT_USED_3, -NOT_USED_4, -NOT_USED_5)

#-------------------------------------------------------------------------------
# Find all FASTQ files

hiseq_run_dirs <- c(file.path(extdir, "HiSeq233", "RNAseq"))
stopifnot(all(dir.exists(hiseq_run_dirs)))

fastq_dirs <- list.dirs(hiseq_run_dirs, recursive = FALSE)
names(fastq_dirs) <- basename(fastq_dirs)


#-------------------------------------------------------------------------------
# Construct data frame and write to disk

r1 <- lapply(fastq_dirs, function(dir) {
  list.files(path = dir, pattern = "^.*R1.*\\.fastq\\.gz$", full.names = TRUE)
})
r2 <- lapply(fastq_dirs, function(dir) {
  list.files(path = dir, pattern = "^.*R2.*\\.fastq\\.gz$", full.names = TRUE)
})

sample <- str_extract(basename(unlist(lapply(r1, "[[", 1))),
                      "[0-9]+-[:alnum:]+-[a-z]+")
replicate <- ifelse(grepl("REP", basename(unlist(lapply(r1, "[[", 1)))),
                    "rep2",
                    "rep1")
prefix <- paste(sample, replicate, sep = ".")

tbl <- data_frame(PREFIX = prefix,
                  R1 = unlist(lapply(r1, paste0, collapse = " ")),
                  R2 = unlist(lapply(r2, paste0, collapse = " ")),
                  SAMPLE = sample,
                  REPLICATE = replicate)

write_tsv(tbl, file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                         "RNA-seq-samples.txt"))

#-------------------------------------------------------------------------------
# Construct colData and save to disk

# NOTE: This includes full path (on JHPCE) to FASTQ and BAM files so that
#       object isn't as dependent on the current working directory

# TODO: Need to fix DONOR_COLOR

cd <- tbl %>%
  mutate(SAMPLE = gsub("-", "_", SAMPLE),
         R1_FULL_PATH = gsub("../extdata",
                             "/dcs01/feinberglab/core/sequencing/hiseq/",
                             R1),
         R2_FULL_PATH = gsub("../extdata",
                             "/dcs01/feinberglab/core/sequencing/hiseq/",
                             R2),
         BAM = file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                         "bam", paste0(PREFIX, ".sorted.bam")),
         BAM_FULL_PATH = gsub("../extdata",
                              "/dcl01/hansen/data",
                              BAM),
         SALMON_DIR = file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                                "salmon",
                                paste0(PREFIX, ".transcripts_quant")),
         SALMON_DIR_FULL_PATH = gsub("../extdata",
                                     "/dcl01/hansen/data",
                                     SALMON_DIR),
         DONOR = substr(SAMPLE, 1, 4),
         TISSUE = vapply(strsplit(tbl$SAMPLE, "-"), "[[", character(1), 2),
         NEUN = vapply(strsplit(tbl$SAMPLE, "-"), "[[", character(1), 3),
         NEUN_COLOR = ifelse(NEUN == "pos", "deepskyblue", "firebrick3"),
         TISSUE_COLOR = ifelse(TISSUE == "BA9", "deepskyblue", "chocolate1")) %>%
  inner_join(data_frame(DONOR = c("5343", "5347", "5358", "5404", "5456",
                                  "5540", "5546"),
                        DONOR_COLOR = c("aquamarine", "chocolate1",
                                        "firebrick3", "deepskyblue",
                                        "purple", "darkgrey", "black")),
             c("DONOR" = "DONOR")) %>%
  inner_join(sample_sheet, c("SAMPLE" = "SAMPLE")) %>%
  inner_join(data_frame(DONOR = c("5347", "5404", "5358", "5343", "5456",
                                   "5540"),
                         FLOW_DATE = factor(c(1, 1, 1, 2, 2, 2))),
             c("DONOR" = "DONOR")) %>%
  mutate(SAMPLE = gsub("_", "-", SAMPLE)) %>%
  DataFrame()

rownames(cd) <- cd$PREFIX

saveRDS(cd, "../objects/colData-flow-sorted-brain-rna-seq.rds")

