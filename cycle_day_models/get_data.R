############################################################
## Get data for cycle day models
############################################################

library(dplyr)
library(stringr)

data_dir <- "/Users/Jessica/Work/2016_mivf/r_projects/data/"

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# LOAD DATA

# Expression data
load(paste0(data_dir, "combined_datasets/combat_exprs.rda"))

# Phenotype info
phenotype <- read.table(
    paste0(data_dir, "combined_datasets/phenotype_2016-09-27.tsv"), 
    header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  mutate(endo=factor(endo))
rownames(phenotype) <- phenotype$sample_id

# Check sample order
stopifnot(colnames(combat_exprs) == rownames(phenotype))

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
phenotype_data <- phenotype %>% filter(! is.na(svm_predicted_day))

# Check all samples in phenotype dataframe have valid 28 day cycle values
stopifnot(phenotype_data$day_cycle %in% seq(0.5,28,0.5))
stopifnot(phenotype_data$svm_predicted_day %in% seq(0.5,28,0.5))

# Subset expression data
combat_exprs <- combat_exprs[,phenotype_data$sample_id]

# Scale expression around zero and make std dev one
scaled_exprs <- combat_exprs %>% t %>% scale %>% t

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
save(combat_exprs, scaled_exprs, phenotype, probe_df, probe_list,
     file="data/cycle_data.RData")
