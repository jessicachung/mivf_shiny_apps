# Create RData object for bmi_models

library(readxl)
library(dplyr)
library(stringr)

cc_exprs_filename <- "/Users/jess/Work/2018_mivf/mivf_endometriosis/data/array_data/cc_combat_exprs.rds"
phenotype_filename <- "/Users/jess/Work/2018_mivf/mivf_endometriosis/data/array_data/combat_phenotype.rds"
bmi_data_filename <- "/Users/jess/Work/2018_mivf/mivf_analysis/data/sarah/gene array list for Jessica - endometrium and BMI analysis (July 2017).xlsx"

cc_exprs <- readRDS(cc_exprs_filename)
combat_phenotype <- readRDS(phenotype_filename)
bmi_data <- read_excel(bmi_data_filename)
prior_genes <- readRDS("data/prior_genes.rds")
probe_info <- readRDS("data/probe_info.rds")

# Parse BMI xlsx data
colnames(bmi_data) <- c("sample_id", "bmi", "group", "subgroup")
bmi_data <- bmi_data %>%
  mutate(sample_id=paste0("X", sample_id),
         group=group %>% tolower %>% str_replace_all("-", ""),
         subgroup=subgroup %>% tolower %>% str_replace_all(" ", "_")) %>%
  mutate(subgroup=ifelse(group == "obese", paste(group, subgroup, sep="_"), subgroup),
         subgroup=ifelse(is.na(subgroup), group, subgroup))

# Merge BMI data with other phenotype data
phenotype <- merge(bmi_data, combat_phenotype, by="sample_id") %>%
  filter(! is.na(cycle_stage)) %>%
  filter(sample_id %in% colnames(cc_exprs))

# Rename group column to class
phenotype <- phenotype %>% rename(class=group)

# Add text column
text <- with(phenotype,
             sprintf("%s [endo: %d, afs: %d]", sample_id, endo, afs_score))
phenotype$text <- text

save(cc_exprs, phenotype, prior_genes, probe_info, file="data/bmi_models.RData")
