# Download data for use with LDSC (https://github.com/bulik/ldsc)
# Peter Hickey
# 2018-06-11

# ------------------------------------------------------------------------------
# Shell variables
#

# NOTE: Need to manually specify path to ldsc directory
LDSC="/users/phickey/software/ldsc"

PHASE1="/dcl01/hansen/data/ldsc/Phase1"
mkdir -p ${PHASE1}
GWASSS="/dcl01/hansen/data/ldsc/GWAS_summary_stats"
mkdir -p ${GWASSS}
MUNGEDSS="/dcl01/hansen/data/ldsc/munge_sumstats"
mkdir -p ${MUNGEDSS}

### ============================================================================
### Download data for use with LDSC
###

# Baseline model LD scores; this contains the `baseline.*` files
wget -O ${PHASE1}/1000G_Phase1_baseline_ldscores.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_baseline_ldscores.tgz

# Regression weights; this contains the `weights.*` files
wget -O ${PHASE1}/weights_hm3_no_hla.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz

# Allele frequencies; this contains the `1000G.mac5eur.*` files
wget -O ${PHASE1}/1000G_Phase1_frq.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_frq.tgz

# Cell type files
wget -O ${PHASE1}/1000G_Phase1_cell_type_groups.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_cell_type_groups.tgz

# HapMap3 SNPs
wget -O ${PHASE1}/w_hm3.snplist.bz2 \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
wget -O ${PHASE1}/hapmap3_snps.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz

# PLINK files
wget -O ${PHASE1}/1000G_Phase1_plinkfiles.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_plinkfiles.tgz

### ============================================================================
### Download GWAS summary stats
###

# ==============================================================================
# Trying to get all data in 'Data Sources' table of
# https://doi.org/10.1101/048991
#

# ------------------------------------------------------------------------------
# Psychiatric disorders
#

# ADHD - PGC-ADD2 (June 2017)
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: package.zip
# MD5 Checksum: c81ea9633e88b3513ccf634b5ada40d6

# Anorexia nervosa - PGC-ED
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: pgc.ed.freeze1.July2017.zip
# MD5 Checksum: a761bdacc52f2543fc2cdce078096d4b

# Anxiety disorder - ANGST
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: angst.study.results.zip
# MD5 Checksum: 6bc88d2c9d3385a6acdd5a1bc209d888

# Autism spectrum disorders(65) - PGC-AUT (2017)
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: daner_AUT_meta14.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tar.gz
# MD5 Checksum: 08d87fa42c7e926533f6a749be677705

# Bipolar disorder - PGC-BIP2
# TODO: Unavailable from http://www.med.unc.edu/pgc/results-and-downloads

# Bipolar disorder (older version)
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: pgc.bip.2012-04.zip
# MD5 Checksum: a675a65c7ff3eb824f90910e8ffc0fdd

# Major depressive disorder - PGC-MDD2
# TODO: Unavailable from http://www.med.unc.edu/pgc/results-and-downloads

# Major depressive disorder (older)
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: pgc.mdd.2012-04.zip
# MD5 Checksum: 895b82d43e11c5e915c34cab5a5e3335

# OCD - PGC-OCDTS
# TODO: Unavailable from http://www.med.unc.edu/pgc/results-and-downloads

# PTSD - PGC-PTSD
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: All.zip
# MD5 Checksum: 430915c1079cc33fb82e30f9d5b159f6

# Schizophrenia - PGC-SCZ2
# Downloaded manually from http://www.med.unc.edu/pgc/results-and-downloads and
# copied to ${GWASSS}
# Filename: ckqny.scz2snpres.gz
# MD5 Checksum: af7b9b521a196ce711d99060426fe01e

# Tourette Syndrome - TSAIGC
# TODO: Unavailable from http://www.med.unc.edu/pgc/results-and-downloads

# ------------------------------------------------------------------------------
# Neurological disorders
#

# Alzheimer's disease - IGAP
# Downloaded manually from
# http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php and
# copied to ${GWASSS}
# Filename: IGAP_summary_statistics.zip
# MD5 Checksum: 99459f68fd909a596d58ec4ed5a4713f

# Epilepsy and subtypes, focal and generalized - ILAE
# Downloaded manually from http://www.epigad.org/gwas_ilae2014/ and copied to
# ${GWASSS}
# There are 3 files, one for each trait ('overall', 'focal', 'GGE')
# Filename: ILAE_All_Epi_11.8.14.txt.gz
# MD5 Checksum: None supplied
# Filename: ILAE_Focal_5.8.14.txt.gz
# MD5 Checksum: None supplied
# Filename: ILAE_GGE_5.8.14.txt.gz
# MD5 Checksum: None supplied

# Intracerebral hemorrhage - ISGC
# Couldn't find data on http://www.strokegenetics.com/
# But did find what seems to be the correct file at
# http://www.cerebrovascularportal.org/informational/downloads and copied to
# ${GWASSS}
# Filename: 3980413.Woo.2014.zip
# MD5 Checksum: None suplied

# Ischemic stroke and subtypes (cardioembolic, early-onset, small-vessel and
# large-vessel)
# Couldn't find data on http://www.strokegenetics.com/
# But did find what seems to be the correct file at
# http://www.cerebrovascularportal.org/informational/downloads and copied to
# ${GWASSS}
# Filename: 4818561.Malik.2016.zip
# MD5 Checksum: None supplied

# Migraine and subtypes, migraine with and without aura - IHGC
# TODO: Data unavailable from
# http://www.headachegenetics.org/content/datasets-and-cohorts

# Multiple sclerosis - IMSGC
# TODO: Data unavailable from http://eaglep.case.edu/imsgc_web

# Parkinson's disease - IPDGC
# TODO: Data unavailable from http://www.pdgene.org/

# ------------------------------------------------------------------------------
# Behavioural-cognitive disorders
#

# College attainment, years of education - SSGAC
# Data manually downloaded from https://www.thessgac.org/data
# (http://ssgac.org/documents/SSGAC_Rietveld2013.zip) and copied to ${GWASSS}
# Filename: SSGAC_Rietveld2013.zip
# MD5 Checksum: None supplied

# Childhood cognitive performance - SSGAC
# Data manually downloaded from https://www.thessgac.org/data
# (http://ssgac.org/documents/CHIC_Summary_Benyamin2014.txt.gz) and copied to
# ${GWASSS}
# Filename: CHIC_Summary_Benyamin2014.txt.gz
# MD5 Checksum: None supplied

# Extraversion, agreeableness, conscientiousness and openness - GPC
# Data manually donwloaded from http://www.tweelingenregister.org/GPC/ and
# copied to ${GWASSS}
# Filename: GPC-1.BigFiveNEO.zip
# MD5 Checksum: None supplied

# IQ - CTG
# Data manually downloaded from http://ctg.cncr.nl/software/summary_statistics
# (http://ctg.cncr.nl/documents/p1651/sumstats.txt.gz) and copied to ${GWASSS}
# Filename: sumstats.txt.gz
# MD5 Checksum: None supplied

# Neuroticism, depressive symptoms and subjective well-being - SSGAC
# Data manually downloaded from https://www.thessgac.org/data and copied to
# ${GWASSS}
# There are 3 files, one for each trait
# (http://ssgac.org/documents/SWB_Full.txt.gz,
# http://ssgac.org/documents/Neuroticism_Full.txt.gz,
# http://ssgac.org/documents/DS_Full.txt.gz)
# Filename: SWB_Full.txt.gz
# MD5 Checksum: None supplied
# Filename: Neuroticism_Full.txt.gz
# MD5 Checksum: None supplied
# Filename: DS_Full.txt.gz
# MD5 Checksum: None supplied

# Never/ever smoked, cigarettes per day - TAG
# Data manually downloaded from
# http://www.med.unc.edu/pgc/results-and-downloads and copied to ${GWASSS}
# There are 2 files, one for each trait
# Filename: tag.evrsmk.tbl.gz
# MD5 Checksum: 9216b2742390819b4d5d0e4673926885
# Filename: tag.cpd.tbl.gz
# MD5 Checksum: ee43737cec4f1b184074ce07ea600b85

# ------------------------------------------------------------------------------
# Additional phenotypes
#

# BMI - GIANT
# Data manually downloaded from
# https://portals.broadinstitute.org/collaboration/giant/index.php/Main_Page
# (https://portals.broadinstitute.org/collaboration/giant/images/f/f0/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz)
# and copied to ${GWASSS}
# Filename: All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz
# MD5 Checksum: None supplied

# Height - GIANT
# Data manually downloaded from
# https://portals.broadinstitute.org/collaboration/giant/index.php/Main_Page
# (https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz)
# and copied to ${GWASSS}
# Filename: GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
# MD5 Checksum: None supplied

# Crohn’s disease - IIBDGC
# Data manually downloaded from https://www.ibdgenetics.org/downloads.html
# (ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/iibdgc-trans-ancestry-filtered-summary-stats.tgz)
# and copied to ${GWASSS}
# Filename: iibdgc-trans-ancestry-filtered-summary-stats.tgz
# MD5 Checksum: None supplied

# Coronary artery disease - Cardiogram
# Data manually downloaded from http://www.cardiogramplusc4d.org/
# (http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cardiogram_gwas_results.zip)
# and copied to ${GWASSS}
# Filename: cardiogram_gwas_results.zip
# MD5 Checksum: None supplied

### ============================================================================
### Download baseline BED files

wget -O ${Phase1}/baseline_bedfiles.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/baseline_bedfiles.tgz

### ============================================================================
### Inflate and extract downloaded data
###

# ------------------------------------------------------------------------------
# Phase1
#

tar xvfz ${PHASE1}/1000G_Phase1_baseline_ldscores.tgz -C ${PHASE1}
tar xvfz ${PHASE1}/weights_hm3_no_hla.tgz -C ${PHASE1}
tar xvfz ${PHASE1}/1000G_Phase1_frq.tgz -C ${PHASE1}
tar xvfz ${PHASE1}/1000G_Phase1_cell_type_groups.tgz -C ${PHASE1}
bzip2 -d ${PHASE1}/w_hm3.snplist.bz2
tar xvfz ${PHASE1}/hapmap3_snps.tgz -C ${PHASE1}
tar xvfz ${PHASE1}/1000G_Phase1_plinkfiles.tgz -C ${PHASE1}

# ------------------------------------------------------------------------------
# GWAS summary stats
#

unzip -d ${GWASSS} ${GWASSS}/package.zip
gzip -d ${GWASSS}/adhd_jul2017.gz
unzip -d ${GWASSS} ${GWASSS}/pgc.ed.freeze1.July2017.zip
unzip -d ${GWASSS} ${GWASSS}/angst.study.results.zip
gzip -d ${GWASSS}/anxiety.meta.full.fs.tbl.gz
gzip -d ${GWASSS}/anxiety.meta.full.cc.tbl.gz
tar xvfz ${GWASSS}/daner_AUT_meta14.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tar.gz -C ${GWASSS}
gzip -d ${GWASSS}/daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv.gz
gzip -d ${GWASSS}/daner_AUT_meta14_WW_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv.gz
gzip -d ${GWASSS}/daner_AUT_meta14_WW_noswm3.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv.gz
unzip -d ${GWASSS} ${GWASSS}/pgc.bip.2012-04.zip
unzip -d ${GWASSS} ${GWASSS}/pgc.mdd.2012-04.zip
unzip -d ${GWASSS} ${GWASSS}/All.zip
gzip -d ${GWASSS}/ckqny.scz2snpres.gz
unzip -d ${GWASSS} ${GWASSS}/IGAP_summary_statistics.zip
gzip -d ${GWASSS}/ILAE_All_Epi_11.8.14.txt.gz
gzip -d ${GWASSS}/ILAE_Focal_5.8.14.txt.gz
gzip -d ${GWASSS}/ILAE_GGE_5.8.14.txt.gz
unzip -d ${GWASSS} ${GWASSS}/3980413.Woo.2014.zip
unzip -d ${GWASSS} ${GWASSS}/4818561.Malik.2016.zip
unzip -d ${GWASSS} ${GWASSS}/SSGAC_Rietveld2013.zip
gzip -d ${GWASSS}/CHIC_Summary_Benyamin2014.txt.gz
unzip -d ${GWASSS} ${GWASSS}/GPC-1.BigFiveNEO.zip
gzip -d ${GWASSS}/sumstats.txt.gz
gzip -d ${GWASSS}/SWB_Full.txt.gz
gzip -d ${GWASSS}/Neuroticism_Full.txt.gz
gzip -d ${GWASSS}/DS_Full.txt.gz
gzip -d ${GWASSS}/tag.evrsmk.tbl.gz
gzip -d ${GWASSS}/tag.cpd.tbl.gz
gzip -d ${GWASSS}/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz
gzip -d ${GWASSS}/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
tar xvfz ${GWASSS}/iibdgc-trans-ancestry-filtered-summary-stats.tgz -C ${GWASSS}
unzip -d ${GWASSS} ${GWASSS}/cardiogram_gwas_results.zip

# ------------------------------------------------------------------------------
# Baseline BED files
#

tar xvfz ${PHASE1}/baseline_bedfiles.tgz -C ${PHASE1}

### ============================================================================
### Convert GWAS summary stats to `.sumstats` format
###

# ------------------------------------------------------------------------------
# Psychiatric disorders
#

# ADHD - PGC-ADD2 (June 2017)
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/adhd_jul2017 \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/ADHD \
  --N-cas 12645 \
  --N-con 84435 \
  --a1-inc

# Anorexia nervosa - PGC-ED
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/pgc.ed.freeze1.summarystatistics.July2017.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Anorexia_nervosa \
  --N-cas 3495 \
  --N-con 11105

# Anxiety disorder - ANGST
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/anxiety.meta.full.cc.tbl \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Anxiety_disorder \
  --N-cas 5761  \
  --N-con 11765

# Autism spectrum disorders(65) - PGC-AUT (2017)
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
# NOTE: Inferred column names from data-release_Jun2017.readme
echo -e "chr\tbp_hg19\tsnp\ta1\ta2\tor\tlb95\tub95\teffect\tse\tp\tfrq_a1\tinfo\tN\tdirection" | \
    cat - ${GWASSS}/daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv > \
    tmp.tsv
python ${LDSC}/munge_sumstats.py \
  --sumstats tmp.tsv \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Autism_spectrum_disorder \
  --N-cas 6197  \
  --N-con 7377 \
  --signed-sumstats or,1
rm tmp.tsv

# Bipolar disorder (older version)
# NOTE: --N 16731 based on Supplementary Table 3 of Finucane et al.
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/pgc.bip.full.2012-04.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Bipolar_disorder \
  --N 16731

# Major depressive disorder (older)
# NOTE: --N-cas/--N-con based on abstract of doi: 10.1038/mp.2012.21
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/pgc.mdd.full.2012-04.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Major_depressive_disorder \
  --N-cas 9240 \
  --N-con 9519

# PTSD - PGC-PTSD
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/SORTED_PTSD_EA9_AA7_LA1_SA2_ALL_study_specific_PCs1.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/PTSD \
  --N-cas 2424 \
  --N-con 7113

# Schizophrenia - PGC-SCZ2
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/ckqny.scz2snpres \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Schizophrenia \
  --N-cas 33640 \
  --N-con 43456

# ------------------------------------------------------------------------------
# Neurological disorders
#

# Alzheimer's disease
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/IGAP_stage_1.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Alzheimers_disease \
  --N-cas 17008 \
  --N-con 37154

# Epilepsy and subtypes, focal and generalized - ILAE
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/ILAE_All_Epi_11.8.14.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Epilepsy \
  --N-cas 7779 \
  --N-con 20439
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/ILAE_Focal_5.8.14.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Focal_epilepsy \
  --N-cas 4601 \
  --N-con 17985
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/ILAE_GGE_5.8.14.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Generalized_epilepsy \
  --N-cas 2525 \
  --N-con 16244

# Intracerebral hemorrhage
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/ICH_GWAS_phase1_finalresults_allICH \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Intracarebral_hemorrhage \
  --N-cas 1545 \
  --N-con 1481

# Ischemic stroke and subtypes (cardioembolic, early-onset, small-vessel and large-vessel)
# NOTE: --N-cas/--N-con based on Table 1 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/metastroke.all.chr.bp \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Ischemic_stroke \
  --N-cas 10307 \
  --N-con 19326
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/metastroke.ce.chr.bp \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Cardioembolic_stroke \
  --N-cas 1859 \
  --N-con 17708
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/metastroke.lvd.chr.bp \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Large-vessel_disease \
  --N-cas 1817 \
  --N-con 17708
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/metastroke.svd.chr.bp \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Small-vessel_disease \
  --N-cas 1349 \
  --N-con 19326

# ------------------------------------------------------------------------------
# Behavioural-cognitive disorders
#

# College attainment, years of education - SSGAC
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/SSGAC_College_Rietveld2013_publicrelease.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/College_attainment \
  --N 120917
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/SSGAC_EduYears_Rietveld2013_publicrelease.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Years_of_education \
  --N 293723

# Childhood cognitive performance - SSGAC
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/CHIC_Summary_Benyamin2014.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Childhood_cognitive_performance \
  --N 17989 \
  --signed-sumstats EFFECT_A1,0 \
  --frq FREQ_A1

# Extraversion, agreeableness, conscientiousness and openness - GPC
# NOTE: --N based on Table 2 of Anttila et al. (2017)
# NOTE: Inferred column names from ReadmeGPC-1.pdf
echo -e "SNPID\tCHR\tBP\tA1\tA2\tBETA\tSE\tPVALUE\tINFO\tNCOH\tMAF" | \
    cat - ${GWASSS}/GPC-1.NEO-EXTRAVERSION.full.txt > \
    tmp.txt
python ${LDSC}/munge_sumstats.py \
  --sumstats tmp.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Extraversion \
  --N 63030
rm tmp.txt
echo -e "SNPID\tCHR\tBP\tA1\tA2\tBETA\tSE\tPVALUE\tINFO\tNCOH\tMAF" | \
    cat - ${GWASSS}/GPC-1.NEO-AGREEABLENESS.full.txt > \
    tmp.txt
python ${LDSC}/munge_sumstats.py \
  --sumstats tmp.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Agreeableness \
  --N 17375
rm tmp.txt
echo -e "SNPID\tCHR\tBP\tA1\tA2\tBETA\tSE\tPVALUE\tINFO\tNCOH\tMAF" | \
    cat - ${GWASSS}/GPC-1.NEO-CONSCIENTIOUSNESS.full.txt > \
    tmp.txt
python ${LDSC}/munge_sumstats.py \
  --sumstats tmp.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Conscientiousness \
  --N 17375
rm tmp.txt
echo -e "SNPID\tCHR\tBP\tA1\tA2\tBETA\tSE\tPVALUE\tINFO\tNCOH\tMAF" | \
    cat - ${GWASSS}/GPC-1.NEO-OPENNESS.full.txt > \
    tmp.txt
python ${LDSC}/munge_sumstats.py \
  --sumstats tmp.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Openness \
  --N 17375
rm tmp.txt

# IQ - CTG
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/sumstats.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/IQ \
  --N 78308 \
  --signed-sumstats Beta,0 \
  --a1 ref \
  --a2 alt

# Neuroticism, depressive symptoms and subjective well-being - SSGAC
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/Neuroticism_Full.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Neuroticism \
  --N 170911
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/DS_Full.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Depressive_symptoms \
  --N 161460
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/SWB_Full.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Subjective_well-being \
  --N 298420

# Never/ever smoked, cigarettes per day - TAG
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/tag.evrsmk.tbl \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Ever_smoked \
  --N 74035 \
  --signed-sumstats OR,0
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/tag.cpd.tbl \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Cigarettes_per_day \
  --N 38617 \
  --signed-sumstats OR,0

# ------------------------------------------------------------------------------
# Additional phenotypes
#

# BMI - GIANT
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/BMI \
  --N 339224

# Height - GIANT
# NOTE: --N based on Table 2 of Anttila et al. (2017)
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Height \
  --N 253288

# Crohn’s disease - IIBDGC
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/EUR.CD.gwas_info03_filtered.assoc \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Crohns_disease \
  --N 20883

# Coronary artery disease - Cardiogram
python ${LDSC}/munge_sumstats.py \
  --sumstats ${GWASSS}/CARDIoGRAM_GWAS_RESULTS.txt \
  --merge-alleles ${PHASE1}/w_hm3.snplist \
  --out ${MUNGEDSS}/Coronary_artery_disease \
  --N 86995
