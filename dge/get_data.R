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
phenotype_all <- readRDS("../../mivf_rna_seq/data/tidy_data/sample_info.rds")
qld_counts_list <- readRDS("../../mivf_rna_seq/data/tidy_data/qld_counts_list.rds")

# Get samples
phenotype <- phenotype_all %>% 
  filter(! str_detect(sample_id, "_2$")) %>%
  filter(cycle_stage %in% 1:7)

# Get gene info
gene_info <- qld_counts_list$gene_info

# Filter counts (liberally)
y <- qld_counts_list$counts[,phenotype$sample_id]
y <- y[rowMeans(cpm(y) > 0.5) > 0.2,]
dim(y)

# Remove batch and cycle effects
y <- DGEList(counts=y)
y <- calcNormFactors(y, method="TMM")
v <- voom(y)
cycle <- factor(phenotype$cycle_stage)
batch <- factor(phenotype$flowcell_batch)
mod <- model.matrix(~cycle)
rna_bc <- ComBat(dat=v$E, batch=batch, mod=mod)
rna_cc <- ComBat(dat=rna_bc, batch=cycle, mod=NULL)
rna_phenotype <- phenotype
rna_gene_info <- gene_info

############################################################
## Microarray data

# Load data
phenotype <- readRDS("../../mivf_microarray/data/phenotype_data/combat_phenotype.rds")
probe_info <- readRDS("../../mivf_endometriosis/data/array_data/illumina_v4_annotation.rds")
array_bc <- readRDS("../../mivf_microarray/data/array_data/combined_combat_exprs_corrected.rds")
array_cc <- readRDS("../../mivf_microarray/data/array_data/day_normalised_exprs.rds")
pvals <- readRDS("../../mivf_microarray/data/array_data/combined_detection_pvals.rds")

# Filter probes with poor detection p-values (liberally)
ok <- rownames(pvals)[rowMeans(pvals < 0.05) > 0.2]
# length(ok)
array_bc <- array_bc[ok,]
array_cc <- array_cc[ok,]

# Remove unnecessary columns
array_phenotype <- phenotype %>% 
  select(sample_id, model_day, batch, study, day_of_cycle)
array_probe_info <- probe_info %>%
  select(IlluminaID, SymbolReannotated, EntrezReannotated) %>%
  rename(illumina_id="IlluminaID",
         entrez="EntrezReannotated",
         symbol="SymbolReannotated")

# Not sure if some of the replicate samples should be removed, but for now
# keep them

# Wrangle phenotype info...
# .... Can add more phenotype info later

############################################################

# Save data
save(rna_phenotype, rna_gene_info, rna_bc, rna_cc, 
     array_phenotype, array_probe_info, array_bc, array_cc,
     file="data/data.RData")

