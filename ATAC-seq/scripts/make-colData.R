# Make colData object for flow sorted brain GTEx project bulk ATAC-seq data to
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
sample_sheet <- read_csv(file.path(extdir, "LR_Edit_08_25_and_08_31_2016_sorted_nuclei_index_table.csv"))

# Sanitise colnames
colnames(sample_sheet) <- c("SAMPLE_ID", "INDEX",
                            "PREFIX", "NUCLEI_COUNT", "EXP_BIOANALYZER_DATE",
                            "POOLS", "NOT_USED1", "NOT_USED2", "NOT_USED3")

# Drop unused columns
sample_sheet <- sample_sheet %>%
  select(-NOT_USED1, -NOT_USED2, -NOT_USED3)

# Fix up PREFIX
sample_sheet$PREFIX <- gsub("NA-", "NA-neg", sample_sheet$PREFIX)
sample_sheet$PREFIX <- gsub("NA\\+", "NA-pos", sample_sheet$PREFIX)
sample_sheet$PREFIX <- gsub("BA9-", "BA9-neg", sample_sheet$PREFIX)
sample_sheet$PREFIX <- gsub("BA9\\+", "BA9-pos", sample_sheet$PREFIX)
sample_sheet$PREFIX <- gsub("REP", ".rep2", sample_sheet$PREFIX)
sample_sheet$PREFIX <- ifelse(grepl("rep2", sample_sheet$PREFIX),
                              sample_sheet$PREFIX,
                              paste0(sample_sheet$PREFIX, ".rep1"))

#-------------------------------------------------------------------------------
# Find all FASTQ files

hiseq_run_dirs <- c(file.path(extdir, "HiSeq234", "Bulk-ATAC-seq"),
                    file.path(extdir, "HiSeq233", "Bulk-ATAC-seq"))
stopifnot(all(dir.exists(hiseq_run_dirs)))

fastq_dirs <- list.dirs(hiseq_run_dirs, recursive = FALSE)
names(fastq_dirs) <- basename(fastq_dirs)

r1 <- lapply(fastq_dirs, function(dir) {
  list.files(path = dir, pattern = "^.*R1.*\\.fastq\\.gz$", full.names = TRUE)
})
r2 <- lapply(fastq_dirs, function(dir) {
  list.files(path = dir, pattern = "^.*R2.*\\.fastq\\.gz$", full.names = TRUE)
})

sample <- str_extract(basename(unlist(lapply(r1, "[[", 1))),
                      "[0-9]+-[:alnum:]+-[a-z]+")
# WARNING: Sample with prefix '5546-NA-pos' is a typo. Should be '5456-NA-pos'
#         I haven't changed the FASTQ filenames but all subsequent files use the
#         corrected '5456-NA-pos' prefix
sample <- gsub("5546", "5456", sample)

replicate <- ifelse(grepl("REP", basename(unlist(lapply(r1, "[[", 1)))),
                    "rep2",
                    "rep1")
prefix <- paste(sample, replicate, sep = ".")

tbl <- data_frame(PREFIX = prefix,
                  R1 = unlist(lapply(r1, paste0, collapse = " ")),
                  R2 = unlist(lapply(r2, paste0, collapse = " ")),
                  SAMPLE = sample,
                  REPLICATE = replicate)

# Need to combine across runs at the prefix level
tbl <- bind_rows(lapply(split(tbl, tbl$PREFIX), function(x) {
  data_frame(PREFIX = unique(x$PREFIX),
             R1 = paste0(x$R1, collapse = " "),
             R2 = paste0(x$R2, collapse = " "),
             SAMPLE = unique(x$SAMPLE),
             REPLICATE = unique(x$REPLICATE))
}))

write_tsv(tbl, file.path(extdir, "flow-sorted-brain-atac", "data",
                         "ATAC-seq-samples.txt"))

#-------------------------------------------------------------------------------
# Construct colData and save to disk

# NOTE: This includes full path (on JHPCE) to FASTQ and BAM files so that
#       object isn't as dependent on the current working directory

cd <- tbl %>%
  mutate(R1_FULL_PATH = gsub("../extdata",
                             "/dcs01/feinberglab/core/sequencing/hiseq/",
                             R1),
         R2_FULL_PATH = gsub("../extdata",
                             "/dcs01/feinberglab/core/sequencing/hiseq/",
                             R2),
         BAM = file.path(extdir, "flow-sorted-brain-atac", "data",
                         "bam", paste0(PREFIX, ".markdup.bam")),
         BAM_FULL_PATH = gsub("../extdata",
                              "/dcl01/hansen/data",
                              BAM),
         DONOR = substr(SAMPLE, 1, 4),
         TISSUE = vapply(strsplit(tbl$SAMPLE, "-"), "[[", character(1), 2),
         NEUN = vapply(strsplit(tbl$SAMPLE, "-"), "[[", character(1), 3),
         TISSUE_NEUN = paste0(TISSUE, "_", NEUN),
         NEUN_COLOR = ifelse(NEUN == "pos", "deepskyblue", "firebrick3"),
         TISSUE_COLOR = ifelse(TISSUE == "BA9", "deepskyblue", "chocolate1")) %>%
  inner_join(data_frame(DONOR = c("5343", "5347", "5358", "5404", "5456",
                                  "5540"),
                        DONOR_COLOR = c("aquamarine", "chocolate1",
                                        "firebrick3", "deepskyblue",
                                        "purple", "darkgrey")),
             c("DONOR" = "DONOR")) %>%
  inner_join(data_frame(DONOR = c("5347", "5404", "5358", "5343", "5456",
                                  "5540"),
                        FLOW_DATE = factor(c(1, 1, 1, 2, 2, 2))),
             c("DONOR" = "DONOR")) %>%
  inner_join(sample_sheet, c("PREFIX" = "PREFIX")) %>%
  DataFrame()

rownames(cd) <- cd$PREFIX

saveRDS(cd, "../objects/colData-flow-sorted-brain-atac-seq.rds")

