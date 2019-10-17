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
rna_phenotype <- readRDS("../../mivf_rna_seq/data/tidy_data/combat_phenotype.rds")
rna_bc <- readRDS("../../mivf_rna_seq/data/tidy_data/batch_normalised_exprs.rds")
rna_cc <- readRDS("../../mivf_rna_seq/data/tidy_data/cycle_normalised_exprs.rds")
qld_counts_list <- readRDS("../../mivf_rna_seq/data/tidy_data/qld_counts_list.rds")

rna_gene_info <- qld_counts_list$gene_info %>%
  filter(ensembl_id %in% rownames(rna_bc))

############################################################
## Microarray data

# Load data
phenotype <- readRDS("../../mivf_microarray/data/phenotype_data/combat_phenotype.rds")
probe_info <- readRDS("../../mivf_endometriosis/data/array_data/illumina_v4_annotation.rds")
array_bc <- readRDS("../../mivf_microarray/data/array_data/combined_combat_exprs_corrected.rds")
array_cc <- readRDS("../../mivf_microarray/data/array_data/day_normalised_exprs.rds")
pvals <- readRDS("../../mivf_microarray/data/array_data/combined_detection_pvals.rds")

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
array_bc <- array_bc[ok,]
array_cc <- array_cc[ok,]

# Remove unnecessary columns
array_phenotype <- phenotype %>% 
  select(sample_id, model_day, batch, study, day_of_cycle)
array_probe_info <- probe_info %>%
  select(IlluminaID, EntrezReannotated, SymbolReannotated) %>%
  rename(illumina_id="IlluminaID",
         entrez="EntrezReannotated",
         symbol="SymbolReannotated") %>%
  filter(illumina_id %in% rownames(array_bc))

# Not sure if some of the replicate samples should be removed, but for now
# keep them

# Wrangle phenotype info...
# .... Can add more phenotype info later

############################################################

# Save data
save(rna_phenotype, rna_gene_info, rna_bc, rna_cc, 
     array_phenotype, array_probe_info, array_bc, array_cc,
     version=2,
     file="data/data.RData")

