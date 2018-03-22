# Data for Shiny: Secretory pathology samples

library(dplyr)

# Array data
batch_corrected_filename <- "data/array_data/batch_corrected.rds"

# Phenotype data
endo_pheno_filename <- "data/patient_data/endo_phenotype.rds"
urs_pod_filename <- "data/patient_data/urs_pod.rds"
endo_pod_filename <- "data/patient_data/endo_pod.rds"
path_data_filename <- "data/patient_data/path_data.rds"
urs_pheno_filename <- "data/patient_data/urs_phenotype.rds"
batch_corrected_phenotype_filename <- "data/array_data/phenotype.rds"

# Probe data
probe_info_filename <- "data/array_data/illumina_v4_annotation.rds"

# Protein atlas files
pa_filename <- "/Users/jess/Work/2018_mivf/mivf_woi/cache/protein_atlas_top.rds"

# Model k=4
results_filename <- "../cache/model_gam_k4_reml.rds"

######################################################################

# Load expression data
exprs <- readRDS(batch_corrected_filename)
probes <- rownames(exprs)

# Get R^2 from k=4 model fitting
model_results <- readRDS(results_filename)
r2 <- sapply(probes, function(p) model_results[[p]]$pseudo_r2)

# Load phenotype dataframes
phenotype <- readRDS(batch_corrected_phenotype_filename)
endo_phenotype <- readRDS(endo_pheno_filename)
urs_phenotype <- readRDS(urs_pheno_filename)
probe_info <- readRDS(probe_info_filename) %>%
  select(-OverlappingSNP, -ArrayAddress) %>%
  rename(Symbol=SymbolReannotated)

# Combine phenotype dataframes
phenotype <- phenotype %>% 
  merge(endo_phenotype %>% select(sample_id:endo_group), by="sample_id", all.x=TRUE) %>% 
  merge(urs_phenotype %>% select(sample_id, group), by="sample_id", all.x=TRUE) %>% 
  mutate(endo=factor(endo)) %>%
  rename(urs=group)

# Parse Protein Atlas files
pa_df <- readRDS(pa_filename)
pa_list <- list("tissue_enriched"=list(),
                "group_enriched"=list(),
                "tissue_enhanced"=list())

for (i in seq_len(nrow(pa_df))) {
  x <- pa_df[i,]
  name <- sprintf("%s - %s (%.3f)", x$IlluminaID, x$Gene, x$r)
  # print(name)
  pa_list[[x$group]][[name]] <- x$IlluminaID
}

save(exprs, probes, r2, phenotype, probe_info, pa_list, file="data/shiny.RData")
