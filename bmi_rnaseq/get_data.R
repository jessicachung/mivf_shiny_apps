# Create RData object for bmi_models

library(readxl)
library(dplyr)
library(stringr)
library(limma)
library(edgeR)
library(sva)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load and parse data

qld_counts_filename <- "/Users/jess/Work/2019_mivf/mivf_rna_seq/data/tidy_data/qld_counts_list.rds"
phenotype_filename <- "/Users/jess/Work/2019_mivf/mivf_rna_seq/data/tidy_data/sample_info_2019-02-08.rds"
#phenotype_filename <- "/Users/jess/Work/2018_mivf/mivf_endometriosis/data/array_data/combat_phenotype.rds"
bmi_data_filename <- "/Users/jess/Work/2018_mivf/mivf_analysis/data/sarah/gene array list for Jessica - endometrium and BMI analysis (July 2017).xlsx"

qld_counts <- readRDS(qld_counts_filename)
phenotype <- readRDS(phenotype_filename)
bmi_data <- read_excel(bmi_data_filename)

# Parse BMI xlsx data
colnames(bmi_data) <- c("sample_id", "bmi", "group", "subgroup")
bmi_data <- bmi_data %>%
  mutate(sample_id=paste0("X", sample_id),
         group=group %>% tolower %>% str_replace_all("-", ""),
         subgroup=subgroup %>% tolower %>% str_replace_all(" ", "_")) %>%
  mutate(subgroup=ifelse(group == "obese", paste(group, subgroup, sep="_"), subgroup),
         subgroup=ifelse(is.na(subgroup), group, subgroup))

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Normalise counts and correct for cycle stage

# Filter low counts
ok <- filterByExpr(qld_counts$counts)
counts <- qld_counts$counts[ok,]

# Normalise counts
y <- DGEList(counts=counts)
y <- calcNormFactors(y, method="TMM")
v <- voom(y)
nc <- cpm(y,normalized.lib.sizes = TRUE, log=TRUE)

# Check v$E and nc are similar
# for (i in seq_len(ncol(nc))) {
#   print(cor(v$E[,i], nc[,i]))
# }

# Get cycle stage
cycle_stage <- phenotype$cycle_stage %>% setNames(phenotype$sample_id)
cycle_stage <- cycle_stage[cycle_stage %in% 1:7]

# Check samples without a valid cycle_stage:
# colnames(v$E)[! colnames(v$E) %in% names(cycle_stage)]

ok_samples <- colnames(v$E)[colnames(v$E) %in% names(cycle_stage)]
cycle_stage <- cycle_stage[ok_samples] %>% factor
exprs <- v$E[,ok_samples]
cc_exprs <- ComBat(exprs, batch=cycle_stage, mod=NULL)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Tidy phenotype data

# Samples with seq data but no BMI data
colnames(cc_exprs)[! colnames(cc_exprs) %in% bmi_data$sample_id]

# Samples with BMI data but no seq data
bmi_data$sample_id[! bmi_data$sample_id %in% colnames(cc_exprs)]

# Merge BMI data with other phenotype data
phenotype <- merge(bmi_data, phenotype, by="sample_id") %>%
  filter(! is.na(cycle_stage)) %>%
  filter(sample_id %in% colnames(cc_exprs))

# Rename group column to class
phenotype <- phenotype %>% rename(class=group)

# Add endo status
phenotype <- phenotype %>%
  mutate(endo=ifelse(afs_score > 0, 1, 0))

# Add text column
text <- with(phenotype,
             sprintf("%s [endo: %d, afs: %s]", sample_id, endo, afs_score))
phenotype$text <- text

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Write to output

save(cc_exprs, phenotype, qld_counts, file="data/bmi_rnaseq.RData")
