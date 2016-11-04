############################################################
## Get data for cycle day models
############################################################

library(dplyr)
library(stringr)

data_dir <- "/Users/Jessica/Work/2016_mivf/r_projects/data/"

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# LOAD DATA

# Expression data
load(paste0(data_dir, "combined_datasets/bg_rsn_combat_2batch_exprs.rda"))
load(paste0(data_dir, "combined_datasets/combat_2batch_cycle_exprs.rda"))

# Phenotype info
phenotype <- read.table(
    paste0(data_dir, "combined_datasets/phenotype_v2.tsv"), 
    header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  mutate(endo=factor(endo))
rownames(phenotype) <- phenotype$sample_id

# Probe info
probe_df <- read.table(
      paste0(data_dir, "combined_datasets/",
             "illumina_v4_annotation_with_detection_stats.tsv"), 
      header=TRUE, sep="\t", stringsAsFactors=FALSE) %>% 
  S4Vectors::rename(X10.="10%", X50.="50%", X90.="90%") %>%
  dplyr::select(-OverlappingSNP, -ArrayAddress) %>%
  S4Vectors::rename(SymbolReannotated="Symbol")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SUBSET DATA

# Only get samples with a 28 day cycle value
phenotype_data <- phenotype %>% 
  filter(sample_id %in% colnames(combat_2batch_cycle_exprs),
         ! is.na(day_cycle_v2))

# Check all samples in phenotype dataframe have valid 28 day cycle values
stopifnot(phenotype_data$day_cycle_v2 %in% seq(0.5,28,0.5))
stopifnot(phenotype_data$svm_predicted_day %in% seq(0,28,0.5))

# Subset expression data
combat_2batch_exprs <- bg_rsn_combat_2batch_exprs[,phenotype_data$sample_id]
combat_2batch_cycle_exprs <- combat_2batch_cycle_exprs[,phenotype_data$sample_id]

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PROBES OF INTEREST

# Probes of interest from Jane
probes_raw <- "ILMN_1764096
ILMN_3194087
ILMN_3238259
ILMN_3272768
ILMN_1734552
ILMN_1652431
ILMN_1696183
ILMN_1743290
ILMN_1786015
ILMN_1807529
ILMN_2091454
ILMN_2169736
ILMN_2367126
ILMN_1740706
ILMN_2060719
ILMN_1784217
ILMN_1729033
ILMN_1782743"

# Make list of probes of interest for drop down menu
probes_of_interest <- probes_raw %>% str_split("\n") %>% unlist
probe_list <- list()
for (p in probes_of_interest) {
  probe_list[[p]] <- p
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SAVE RDATA

phenotype <- phenotype_data
save(combat_2batch_exprs, combat_2batch_cycle_exprs, phenotype, probe_df, 
     probe_list, file="data/cycle_data.RData")
