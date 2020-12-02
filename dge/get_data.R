library(limma)
library(edgeR)
library(dplyr)
library(tibble)
library(stringr)
library(reshape2)
library(sva)

############################################################
## RNA-Seq Data

# Load data
rna_phenotype <- readRDS("../../endo_molecular_model/data/tidy_data/combat_phenotype_2020-09-08.rds")
rna_bc <- readRDS("../../endo_molecular_model/data/tidy_data/batch_normalised_exprs.rds")
rna_cc <- readRDS("../../endo_molecular_model/data/tidy_data/cycle_normalised_exprs_2020-09-08_k30_g4.rds")
qld_counts_list <- readRDS("../../endo_molecular_model/data/tidy_data/qld_counts_list.rds")

rna_gene_info <- qld_counts_list$gene_info %>%
  filter(ensembl_id %in% rownames(rna_bc))

# Remove unnecessary columns
rna_phenotype <- rna_phenotype %>% 
  select(sample_id, study, age, endo, cycle_stage, transformed_time_2) %>%
  rename(model_time="transformed_time_2")


############################################################
## Microarray data

# Load data
# array_phenotype <- readRDS("../../mivf_microarray/data/phenotype_data/combat_phenotype.rds")
array_phenotype <- readRDS("../../endo_molecular_model/data/tidy_data/array_combat_phenotype_2020-12-02.rds")
probe_info <- readRDS("../../mivf_endometriosis/data/array_data/illumina_v4_annotation.rds")
array_bc <- readRDS("../../endo_molecular_model/data/array_data/combat_exprs.rds")
array_cc <- readRDS("../../endo_molecular_model/data/tidy_data/array_cycle_normalised_exprs_2020-12-02_k30_g4.rds")
pvals <- readRDS("../../mivf_endometriosis/data/array_data/combined_pval.rds")

# Filter probes with poor detection p-values (liberally) and poor annotation quality
stopifnot(rownames(pvals) == probe_info$IlluminaID)
probe_info <- probe_info %>% 
  mutate(p_detected=rowMeans(pvals < 0.05),
         pass_quality=str_detect(ProbeQuality, "Perfect|Good"),
         pass_detect=p_detected > 0.2)

# Manually add genes of interest
probe_info <- probe_info %>%
  mutate(special=SymbolReannotated == "CAMK4")

# Get probes
ok <- probe_info %>%
  filter(special | (pass_quality & pass_detect)) %>%
  pull(IlluminaID)
# length(ok)

# Get data
array_bc <- array_bc[ok,array_phenotype$sample_id]
array_cc <- array_cc[ok,array_phenotype$sample_id]

# Remove unnecessary columns
# array_phenotype <- phenotype %>% 
#   select(sample_id, model_day, batch, study, day_of_cycle)
array_phenotype <- array_phenotype %>% 
  select(sample_id, study, age, endo, cycle_stage, transformed_time_2) %>%
  rename(model_time="transformed_time_2")
array_probe_info <- probe_info %>%
  select(IlluminaID, EntrezReannotated, SymbolReannotated) %>%
  rename(illumina_id="IlluminaID",
         entrez="EntrezReannotated",
         symbol="SymbolReannotated") %>%
  filter(illumina_id %in% rownames(array_bc))




############################################################

# Save data
save(rna_phenotype, rna_gene_info, rna_bc, rna_cc, 
     array_phenotype, array_probe_info, array_bc, array_cc,
     version=2,
     file="data/data.RData")

