# Make a BSseq object for all CH methylation data in the flow sorting WGBS
# project. This uses the HDF5Array branch of bsseq to created a HDF5-backed
# BSseq object
# Peter Hickey
# 2016-09-07

library(bsseq)
extdir <- "../extdata"
options("mc.cores" = 45)

x <- readRDS(file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                       "BS.fit.small.sorted.somatic.all.hdf5.rds"))
cd <- colData(x)
rm(x)

#-------------------------------------------------------------------------------
# Load data and construct a HDF5-backed BSseq object for CHG loci in each
# seqlevel

context <- "CHG"
reports <- list.files(file.path(extdir, "fornonCG"),
                      pattern = paste0("^.*", context, "_report\\.txt$"),
                      full.names = TRUE)
# Drop the chromosome-specfic files
reports <- grep("chr", reports, value = TRUE, invert = TRUE)
seqlevels <- paste0("chr", c(1:22, "X", "Y"))

dir.create(paste0("/scratch/temp/phickey/", context))
setHDF5DumpDir(paste0("/scratch/temp/phickey/", context))
list_of_chg_bsseq <- lapply(seqlevels, function(seqlevel, context) {
  message(seqlevel)
  # Extract just the loci on the seqlevel of interest and in the context of
  # interest
  parsed <- unlist(mclapply(reports, function(x) {
    tf <- tempfile(pattern = paste0(basename(x), "-", context, "-",
                                    seqlevel, "."),
                   tmpdir = "/scratch/temp/phickey")
    cmd <- paste0("awk '{ if ($1 == \"", seqlevel, "\" && $6 == \"", context,
                  "\") print $0 }' ", x, " > ", tf)
    system(cmd)
    tf
  }))

  # Create the HDF5-backed BSseq object
  bsseq <- read.bismark(files = parsed,
                        sampleNames = gsub(paste0("_", context, "_report.txt"),
                                           "", basename(reports)),
                        rmZeroCov = FALSE,
                        strandCollapse = FALSE,
                        fileType = "cytosineReport",
                        verbose = FALSE,
                        mc.cores = getOption("mc.cores"),
                        hdf5 = TRUE)

  # Remove the temporary files
  unlink(parsed)

  # Return the HDF5-backed BSseq object
  bsseq
}, context = context)

# Consolidate assays into a single .h5 file per object in
# /dcl01/hansen/data/flow-sorted-brain-wgbs/hdf5
list_of_chg_bsseq <- mclapply(list_of_chg_bsseq, function(bsseq, context) {
  # NOTE: Use full path and not relative path for better portability (on JHPCE)
  setHDF5DumpFile(paste0("/dcl01/hansen/data/flow-sorted-brain-wgbs/hdf5/",
                         context, "-", unique(seqnames(bsseq)),
                         "-flow-sorted-brain-wgbs.h5"))
  assays <- endoapply(assays(bsseq), function(a) {
    HDF5Array(a)
  })
  assays(bsseq) <- assays
  message("Done ", unique(seqnames(bsseq)))
  bsseq
}, context = context)

# Add colData
list_of_chg_bsseq <- mclapply(list_of_chg_bsseq, function(bsseq) {
  colData(bsseq) <- cd[match(colnames(bsseq), rownames(cd)), ]
  bsseq
})

# Save seqlevel objects
mclapply(list_of_chg_bsseq, function(bsseq, context) {
  saveRDS(bsseq,
          file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(context, "-", unique(seqnames(bsseq)),
                           "-flow-sorted-brain-wgbs.rds")))
}, context = context)

#-------------------------------------------------------------------------------
# Load data and construct a HDF5-backed BSseq object for CHH loci in each
# seqlevel
# NOTE: Ended up parsing the all files and then creating all the BSseq objects
#       for CHH loci (did this on a seqlevel-basis for CHG loci). This approach
#       is slightly more efficient but requires more disk space.

context <- "CHH"
reports <- list.files(file.path(extdir, "fornonCG"),
                      pattern = paste0("^.*", context, "_report\\.txt$"),
                      full.names = TRUE)
# Drop the chromosome-specfic files
reports <- grep("chr", reports, value = TRUE, invert = TRUE)
seqlevels <- paste0("chr", c(1:22, "X", "Y"))

dir.create(paste0("/scratch/temp/phickey/", context))
setHDF5DumpDir(paste0("/scratch/temp/phickey/", context))

# Create seqlevel-specific CHH report files for each sample
df <- data.frame(report = rep(reports, length(seqlevels)),
                 seqlevel = rep(seqlevels, each = length(reports)),
                 stringsAsFactors = FALSE)
parsed <- unlist(mclapply(seq_len(nrow(df)), function(i) {
  report <- df$report[i]
  seqlevel <- df$seqlevel[i]
  message(seqlevel, ", ", report)
  tf <- tempfile(pattern = paste0(basename(report), "-", context, "-",
                                  seqlevel, "."),
                 tmpdir = "/scratch/temp/phickey")
  cmd <- paste0("awk '{ if ($1 == \"", seqlevel, "\" && $6 == \"", context,
                "\") print $0 }' ", report, " > ", tf)
  system(cmd)
  tf
}))

# Create the HDF5-backed BSseq objects
df <- cbind(df, parsed, stringsAsFactors = FALSE)
df$seqlevel <- factor(df$seqlevel, levels = paste0("chr", c(1:22, "X", "Y")),
                      ordered = TRUE)
split_df <- split(df, df$seqlevel)
list_of_chh_bsseq <- lapply(split_df, function(df_) {
  message(unique(df_$seqlevel))

  read.bismark(files = df_$parsed,
               sampleNames = gsub(paste0("_", context, "_report.txt"), "",
                                  basename(df_$report)),
               rmZeroCov = FALSE,
               strandCollapse = FALSE,
               fileType = "cytosineReport",
               verbose = FALSE,
               mc.cores = getOption("mc.cores"),
               hdf5 = TRUE)
})

# Consolidate assays into a single .h5 file per object in
# /dcl01/hansen/data/flow-sorted-brain-wgbs/hdf5
list_of_chh_bsseq <- mclapply(list_of_chh_bsseq, function(bsseq, context) {
  # NOTE: Use full path and not relative path for better portability (on JHPCE)
  setHDF5DumpFile(paste0("/dcl01/hansen/data/flow-sorted-brain-wgbs/hdf5/",
                         context, "-", unique(seqnames(bsseq)),
                         "-flow-sorted-brain-wgbs.h5"))
  assays <- endoapply(assays(bsseq), function(a) {
    HDF5Array(a)
  })
  assays(bsseq) <- assays
  bsseq
}, context = context)

# Add colData
list_of_chh_bsseq <- mclapply(list_of_chh_bsseq, function(bsseq) {
  colData(bsseq) <- cd[match(colnames(bsseq), rownames(cd)), ]
  bsseq
})

# Save seqlevel objects
mclapply(list_of_chh_bsseq, function(bsseq, context) {
  saveRDS(bsseq,
          file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(context, "-", unique(seqnames(bsseq)),
                           "-flow-sorted-brain-wgbs.rds")))
}, context = context)

#-------------------------------------------------------------------------------
# Combine into a single CH BSseq object (one per seqlevel)

list_of_ch_bsseq <- mcmapply(function(chg, chh) {
  sort(rbind(chg, chh))
}, chg = list_of_chg_bsseq, chh = list_of_chh_bsseq)

# Consolidate assays into a single .h5 file per object in
# /dcl01/hansen/data/flow-sorted-brain-wgbs/hdf5
list_of_ch_bsseq <- mclapply(list_of_ch_bsseq, function(bsseq, context) {
  # NOTE: Use full path and not relative path for better portability (on JHPCE)
  setHDF5DumpFile(paste0("/dcl01/hansen/data/flow-sorted-brain-wgbs/hdf5/",
                         "CH-", unique(seqnames(bsseq)),
                         "-flow-sorted-brain-wgbs.h5"))
  assays <- endoapply(assays(bsseq), function(a) {
    HDF5Array(a)
  })
  assays(bsseq) <- assays
  bsseq
})

# Save seqlevel objects
mclapply(list_of_ch_bsseq, function(bsseq, context) {
  saveRDS(bsseq,
          file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0("CH-", unique(seqnames(bsseq)),
                           "-flow-sorted-brain-wgbs.rds")))
}, context = context)
