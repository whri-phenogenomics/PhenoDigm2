### IMPC post processing - models and PhenoDigm scores
### It requires:
### 1) tables from PhenoDigm sqlite database in Apocrita
### 2) different IMPC files from the ftp repository
### 3) Auxiliary files and scripts: orthologue mapping and symbol checker
### https://www.gentar.org/orthology-api/api/ortholog/one_to_one/impc/write_to_tsv_file
### http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-22.0/results/


# Set library path for the interpreter
.libPaths('/data/WHRI-Phenogenomics/projects/PhenoDigm2/post_processing_dependencies/r_lib_paths/R/x86_64-pc-linux-gnu-library/4.4.1')

# Intake CLI arguments Here we pass the absolute path to the phenodigm DB

args <- commandArgs(trailingOnly = TRUE)
post_processing_path <- args[1]


# Install/Load packages as needed 
packages <- c("R.utils", "dplyr", "readr", "data.table", "arrow", "ggplot2", "fst", "logr")

for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

library("dplyr")
library("readr")
library(data.table)
library(arrow)
library(ggplot2)
library(fst)
library(logr)
library(jsonlite)



print(nslots)
# Set the root directory dynamically
setwd(post_processing_path)

# Execute external Rscript
source("./scripts/auxiliary/hgnc_symbol_checker.R")

# Initialise a logger for output ------------------------------------------

# Create a temporary file location
log_open(
  "./data/output/post_processing_results.log"
  )

# Remove notes and blank lines
# options("logr.notes" = FALSE)
options("logr.notes" = FALSE, "logr.compact" = TRUE)

# Function to log events
log_event <- function(...) {
  record <- list(...)
  msg <- toJSON(record, auto_unbox = TRUE)
  logr::log_print(msg)
}


# import data -------------------------------------------------------------

model_genotype <- fread("./data/phenodigm/model_genotype.tsv.gz") %>% 
  # model_genotype <-read_parquet("./data/phenodigm/model_genotype.parquet") %>%
  select(id,gene_id) %>%
  rename(mgi_id = gene_id)


model <- fread("./data/phenodigm/model.tsv.gz") %>%
  inner_join(model_genotype)

hgnc <- protein_coding_genes %>%
  select(hgnc_id,symbol) %>%
  dplyr::rename(gene_symbol = symbol)


hm_ortho_symbol <- read_delim("./data/auxiliary/one_to_one_orthologs.tsv") %>%
  select(`Human Gene Symbol`,`Hgnc Acc Id`,`Mgi Gene Acc Id`) %>%
  rename(gene_symbol = `Human Gene Symbol`,
         hgnc_id = `Hgnc Acc Id`,
         mgi_id = `Mgi Gene Acc Id`) %>%
  as.data.frame()

disease <- fread("./data/phenodigm/disease.tsv.gz") %>%
  select(id,term) %>%
  dplyr::rename(disorder_id = id, disorder_name = term) 


gene_disease <- read_delim("./data/phenodigm/disease_gene_mapping.tsv.gz", 
                           delim = "\t", col_names = TRUE) 

# Log number of gene disease entries
# put(paste0("Gene disease table: ",length(unique(gene_disease$match))))

log_event(
  gene_disease_table = length(unique(gene_disease$match))
  )

gene_disease <- read_delim("./data/phenodigm/disease_gene_mapping.tsv.gz", 
                           delim = "\t", col_names = TRUE) %>%
  filter(source !="MGI") %>%
  select(query,match) %>%
  right_join(disease,by =c("query" = "disorder_id")) %>%
  dplyr::rename(disorder_id = query, hgnc_id = match) %>%
  distinct() %>%
  left_join(hgnc) %>%
  select(disorder_id, disorder_name, hgnc_id, gene_symbol) %>%
  replace(is.na(.),"-") %>%
  inner_join(hm_ortho_symbol) %>%
  mutate(gene_disease = paste0(hgnc_id, "_", disorder_id)) %>%
  distinct() %>%
  mutate(id = paste(gene_symbol, hgnc_id, mgi_id, sep = "_"))

# Log number of disease genes
log_event(disease_genes = length(unique(gene_disease$gene_disease)))

# Log total number of human disease genes with mouse orthologs
log_event(disease_genes_with_mouse_orthologs = length(unique(gene_disease$hgnc_id)))

gene_disease_orthologs <- gene_disease %>%
  distinct(mgi_id, .keep_all=TRUE) %>%
  filter(!is.na(mgi_id)) 

# # Check how many non-disease genes
# genes <- read_delim("./data/phenodigm/gene.tsv.gz", 
#                     delim = "\t", col_names = TRUE) %>%
#   filter(organism == "human") %>%
#   filter(!id %in% gene_disease$hgnc_id)
# 
# print(length(unique(genes$id)))


# import impc data: phenotype associations, viability, all -----------------
# "Cache" assertions
geno_pheno_assertions <- fread("./data/impc/genotype-phenotype-assertions-ALL.csv.gz") 

# Log total number of genes with geno-pheno
log_event(genes_with_geno_pheno_assertions = length(unique(geno_pheno_assertions$marker_accession_id)))
# Genes with a human ortholog and a geno-pheno annotation
geno_pheno_with_orthologs <- geno_pheno_assertions %>%
  filter(marker_accession_id %in% hm_ortho_symbol$mgi_id)

# Log number of human orthologs and an geno-pheno annotation
log_event(orthologs_with_geno_pheno_assertions = length(unique(geno_pheno_with_orthologs$marker_accession_id)) )

# How many are not in disease genes
geno_pheno_with_orthologs_no_disease_genes <- geno_pheno_with_orthologs %>%
  rename(mgi_id = marker_accession_id) %>%
  inner_join(hm_ortho_symbol) %>%
  filter(!hgnc_id %in% gene_disease$hgnc_id)

# Log number of genes with phenotype hits with no disease associated
log_event(genes_with_phenotype_hits_no_disease_associated = length(unique(geno_pheno_with_orthologs_no_disease_genes$mgi_id)))


assertions_not_in_disease_genes <- geno_pheno_assertions %>%
  filter(!marker_accession_id %in% gene_disease$mgi_id)
# Log number of assertions not in disease genes
log_event(assertions_not_in_disease_genes = length(unique(assertions_not_in_disease_genes$marker_accession_id)))



phenotypes <- geno_pheno_assertions %>%
  dplyr::rename(mgi_id = marker_accession_id) %>%
  select(mgi_id, top_level_mp_term_id, top_level_mp_term_name,
         mp_term_id, mp_term_name) %>%
  filter(!is.na(top_level_mp_term_id)) %>%
  separate_rows(top_level_mp_term_id, top_level_mp_term_name, sep = "\\,") %>%
  distinct() %>%
  select(mgi_id,top_level_mp_term_name) %>%
  arrange(mgi_id,top_level_mp_term_name) %>%
  distinct() %>%
  dplyr::group_by(mgi_id) %>%
  dplyr::mutate(top_level_mp = paste0(top_level_mp_term_name, collapse = "|")) %>%
  select(mgi_id, top_level_mp) %>%
  distinct()


phenotypes_all <- geno_pheno_assertions %>%
  dplyr::rename(mgi_id = marker_accession_id) %>%
  select(mgi_id, top_level_mp_term_id, top_level_mp_term_name,
         mp_term_id, mp_term_name) %>%
  filter(!is.na(top_level_mp_term_id)) %>%
  separate_rows(top_level_mp_term_id, top_level_mp_term_name, sep = "\\,") %>%
  distinct()

viability <- read_delim("./data/impc/viability.csv.gz", 
                        delim = ",", col_names = TRUE) %>%
  filter(is.na(Comment)) %>%
  dplyr::rename(mgi_id = 'Gene Accession Id',
                viability = 'Viability Phenotype HOMs/HEMIs') %>%
  select(mgi_id, viability) %>%
  arrange(mgi_id, viability) %>%
  distinct()

# "Cache" stats results
stats_all <- fread("./data/impc/statistical-results-ALL.csv.gz")

# TODO: Log this in the future, not useful right now.
length(unique(stats_all$marker_accession_id))

# How many orthologs are in the stats
orthologs_in_stats <- stats_all %>%
  filter(marker_accession_id %in% unique(hm_ortho_symbol$mgi_id))

# TODO: Log this in the future, not useful right now.
length(unique(orthologs_in_stats$marker_accession_id))

all_filter <- stats_all %>%
  dplyr::rename(mgi_id = marker_accession_id) %>%
  select(mgi_id) %>%
  distinct() %>%
  left_join(phenotypes) %>%
  left_join(viability) %>%
  replace(is.na(.),"-") %>%
  filter(mgi_id !="-") %>%
  filter(!is.na(mgi_id))


mouse_symbol = stats_all %>%
  dplyr::rename(mgi_id = marker_accession_id,
                mm_gene_symbol = marker_symbol) %>%
  dplyr::select(mgi_id, mm_gene_symbol) %>%
  distinct()

all_disease <- all_filter %>%
  filter(mgi_id %in% unique(gene_disease$mgi_id)) %>%
  mutate(phenotype = ifelse(top_level_mp == "-" & viability =="-","n",
                            ifelse(top_level_mp == "-" & viability =="viable","n","y"))) %>%
  rename(top_level_mp_term_name = top_level_mp)

# TODO: To log these table values appropriately
table(all_disease$phenotype)
table(all_disease$viability)


all_disease_viability <- all_disease %>%
  mutate(type_phenotype = ifelse(top_level_mp_term_name %in% c("mortality/aging","-") &
                                   viability %in% c("lethal","subviable"), 
                                 "lethal_only_phenotype",
                                 ifelse(!top_level_mp_term_name %in% c("mortality/aging","-") &
                                          viability %in% c("lethal","subviable"),"lethal_and_other",
                                        ifelse(top_level_mp_term_name %in% c("mortality/aging","-") &
                                                 viability == "viable", "viable_no_other_phenotype",
                                               ifelse(!top_level_mp_term_name %in% c("mortality/aging","-") &
                                                        viability == "viable","viable_and_phenotype",
                                                      ifelse(!top_level_mp_term_name %in% c("mortality/aging","-") &
                                                               viability =="-","phenotype_no_viability_data",
                                                             ifelse(top_level_mp_term_name =="-" &
                                                                      viability =="-","nada","-")))))))
# TODO: To log these table values appropriately
table(all_disease_viability$type_phenotype)


all_disease_by_gene <- all_disease %>%
  select(mgi_id, viability, phenotype) %>%
  distinct()

all_disease_viability_group <- all_disease_viability %>%
  arrange(mgi_id, type_phenotype) %>%
  group_by(mgi_id) %>%
  summarise(types = paste0(unique(type_phenotype), collapse = "|")) %>%
  mutate(phenotype_type = ifelse(types == "lethal_and_other|lethal_only_phenotype",
                                 "lethal_and_other_phenotype",
                                 ifelse(types == "viable_and_phenotype|viable_no_other_phenotype",
                                        "viable_and_phenotype",
                                        ifelse(types == "-|phenotype_no_viability_data",
                                               "phenotype_no_viability_data",
                                               ifelse(types %in% c("-","nada"),
                                                      "no_phenotype",types))))) %>%
  select(mgi_id,phenotype_type) %>%
  distinct() %>%
  inner_join(all_disease_by_gene)

all_disease_viability_group %>%
  count(phenotype_type)

nomp <- all_disease_viability_group %>%
  filter(phenotype == "y") %>%
  select(mgi_id) %>%
  pull(mgi_id)


hpo <- fread("./data/hpo/genes_to_phenotype.txt")

hpo_genes <- hgnc.checker(unique(hpo$gene_symbol), protein_coding_genes) %>%
  filter(hgnc_id != "-") %>%
  select(hgnc_id) %>%
  distinct() %>%
  mutate(hpo_annotations = "y") %>%
  inner_join(hm_ortho_symbol)

# one-to-one orthologs with IMPC models
disease_genes_impc_categories <- all_disease_viability_group %>%
  left_join(hpo_genes) %>%
  replace(is.na(.),"-")

# Log one-to-one orthologs with IMPC models
log_event(disease_genes_with_orthologs_with_impc_models = dim(disease_genes_impc_categories)[1])

# TODO: Find what information this provides
disease_genes_impc_pairs <- gene_disease %>%
  filter(mgi_id %in% disease_genes_impc_categories$mgi_id)

# Disease genes with orthologs, impc models and mp annot - IMPC mouse phenotypes encoded as MP terms
genes_mouse_phenotypes <- disease_genes_impc_categories %>%
  filter(phenotype == "y")
log_event(disease_genes_orthologs_impc_models_with_mp = dim(genes_mouse_phenotypes)[1])


# Disease genes with orthologs, impc models and hp annot - with human phenotypes encoded as HPO terms
genes_human_phenotypes <- disease_genes_impc_categories %>%
  filter(hpo_annotations == "y")
log_event(disease_genes_orthologs_impc_models_with_hp = dim(genes_human_phenotypes)[1])

# Disease genes with orthologs, impc models and hp and mp annot
genes_mouse_human_phenotypes <- genes_mouse_phenotypes %>%
  filter(hpo_annotations == "y")
log_event(disease_genes_orthologs_impc_models_with_hp_and_mp = dim(genes_mouse_human_phenotypes)[1])

# phenodigm scores --------------------------------------------------------

### genes with phenodigm match
model_disease_omim_score_no0 <- fread("./data/phenodigm/disease_model_association_omim_impc.tsv.gz") %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  inner_join(model, by = c("match" = "id")) %>%
  filter(match !="MGI") %>%
  filter(mgi_id %in% all_disease$mgi_id)%>%
  filter(mgi_id %in% genes_mouse_human_phenotypes$mgi_id)%>%  
  inner_join(hm_ortho_symbol) %>%
  mutate(gene_disease = paste0(hgnc_id,"_", query)) %>%
  inner_join(disease_genes_impc_pairs) %>%
  distinct() #%>%
#filter(!grepl("NOT-RELEASED", match))



model_disease_orphanet_score_no0 <- read_delim("./data/phenodigm/disease_model_association_orphanet_impc.tsv.gz", 
                                               delim = "\t", col_names = TRUE) %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  inner_join(model, by = c("match" = "id")) %>%
  filter(match !="MGI") %>%
  filter(mgi_id %in% all_disease$mgi_id)%>%
  filter(mgi_id %in% genes_mouse_human_phenotypes$mgi_id)%>%  
  inner_join(hm_ortho_symbol) %>%
  mutate(gene_disease = paste0(hgnc_id,"_", query)) %>%
  inner_join(disease_genes_impc_pairs) %>%
  distinct() #%>%
#filter(!grepl("NOT-RELEASED", match))


phenodigm_matches <- model_disease_omim_score_no0 %>%
  bind_rows(model_disease_orphanet_score_no0) %>%
  inner_join(mouse_symbol)

# Genes that recapitulate clinical features of human disease
phenodigm_matches_impc = phenodigm_matches

matches_tidy <- phenodigm_matches_impc %>%
  select(disorder_id, disorder_name, hgnc_id,
         gene_symbol, description,
         score, query_phenotype, match_phenotype)

head(matches_tidy)

# Log number of genes associated to disease that recapitulate human disease in some way
log_event(genes_associated_recapitulate_human_disease = length(unique(phenodigm_matches$hgnc_id)))

# Distribution plots 
base <-  ggplot(matches_tidy, aes(y=score))
histogram <- base +
  geom_histogram(fill = "grey") +
  scale_y_continuous(position = "right")

box <- base + 
  geom_boxplot(fill = "lightblue")

density <- base +
  geom_density(fill = "lightblue", adjust = 1/4)


violin_score_associated <- ggplot(matches_tidy, aes(x = 0, y = score)) +
  geom_violin(fill = "#00AFBB", trim = FALSE, width = 0.3, adjust = 0.3) +
  geom_boxplot(width = 0.02, fill = "white", outlier.shape = NA) +
  stat_summary(
    aes(shape = "Mean ± SD"),
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", color = "black"
  ) +
  geom_hline(
    data = data.frame(yintercept = 40),
    aes(yintercept = yintercept, , linetype = "Score threshold (40)"),
    color = "purple", linewidth = 0.8
  ) +
  scale_shape_manual(values = c("Mean ± SD" = 16)) +
  scale_linetype_manual(values = c("Score threshold (40)" = "dashed")) +
  guides(
    shape = guide_legend(title = NULL),
    linetype = guide_legend(title = NULL, override.aes = list(color = "purple"))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 0),
    legend.position = c(0.9, 0.9)
  ) +
  labs(
    x = "PhenoDigm match hits",
    y = "PhenoDigm Score"
  )


histogram
box
density
violin_score_associated

# Disease gene(associated) summary stats:
summary_stats_associated <- matches_tidy %>%
        summarise(
          min = min(score, na.rm = TRUE),
          q1 = quantile(score, 0.25, na.rm = TRUE),
          median = median(score, na.rm = TRUE),
          mean = mean(score, na.rm = TRUE),
          q3 = quantile(score, 0.75, na.rm = TRUE),
          max = max(score, na.rm = TRUE),
          iqr = IQR(score, na.rm = TRUE)
        )

# Log disease-associated gene stats
log_event(score_disease_associated_min = summary_stats_associated$min)
log_event(score_disease_associated_q1 = summary_stats_associated$q1)
log_event(score_disease_associated_median = summary_stats_associated$median)
log_event(score_disease_associated_mean = summary_stats_associated$mean)
log_event(score_disease_associated_q3 = summary_stats_associated$q3)
log_event(score_disease_associated_max = summary_stats_associated$max)
log_event(score_disease_associated_iqr = summary_stats_associated$iqr)

# export files for shiny app ----------------------------------------------

write_parquet(phenodigm_matches_impc,
              "./data/output/phenodigm_matches_full.parquet",
              compression="zstd", compression_level=5)

write_parquet(matches_tidy,
              "./data/output/phenodigm_matches.parquet",
              compression="zstd", compression_level=5)

# gene summary file -------------------------------------------------------

# Do not keep unique disease names, this helps separating them in the DMP home tab
# NOTE:  Unique id's are kept, names have to be duplicated to get the right amount of rows in the portal. (Two different id's can have the same name.)
gene_disease_by_gene <- gene_disease %>%
  group_by(hgnc_id) %>%
  summarise(disorder_id = paste0(unique(disorder_id), collapse = "|"),
            disorder_name = paste0(disorder_name, collapse = "|"))

max_score <- phenodigm_matches %>%
  group_by(hgnc_id) %>%
  summarise(max_score = max(score))

gene_summary_df <- hm_ortho_symbol %>%
  left_join(gene_disease_by_gene) %>%
  mutate(IMPC_pipeline = ifelse(mgi_id %in% all_filter$mgi_id,"yes","no")) %>%
  mutate(IMPC_phenotypes = ifelse(mgi_id %in% genes_mouse_phenotypes$mgi_id,"yes","no") )%>%
  mutate(HPO_phenotypes = ifelse(mgi_id %in% genes_human_phenotypes$mgi_id,"yes","no") )%>%
  mutate(PhenoDigm_match = ifelse (mgi_id %in% unique(phenodigm_matches$mgi_id),
                                   "yes","no")) %>%
  replace(is.na(.),"-")  %>%
  left_join(max_score) 

filter_matches <- gene_summary_df %>%
  filter(PhenoDigm_match == "yes")


# export gene summary file ------------------------------------------------

write_parquet(gene_summary_df,
              "./data/output/gene_summary.parquet",
              compression="zstd", compression_level=5)

write.fst(gene_summary_df, "./data/output/gene_summary.fst")


# home page gene summary --------------------------------------------------

# Summarise phenodigm match phenotypes by gene symbol and disorder id
gene_disease_home_page <- matches_tidy  %>%
   select(
          disorder_id,
          gene_symbol,
          hgnc_id,
          query_phenotype,
          match_phenotype) %>%
group_by(gene_symbol,disorder_id) %>%
   summarise(
     query_phenotype=paste(query_phenotype, collapse=","),
     match_phenotype=paste(match_phenotype, collapse=",")
   )

# Explode the gene info, to separate diseases by rows
gene_info <- gene_summary_df%>%
  mutate(disease_gene = ifelse(disorder_id != "-", "yes", "no")) %>%
  select(
    gene_symbol,
    disease_gene,
    IMPC_pipeline,
    IMPC_phenotypes,
    PhenoDigm_match,
    max_score,
    disorder_id,
    disorder_name
  ) %>%
  separate_longer_delim(c("disorder_id","disorder_name"), delim ="|") 
 
  # Join by disorder name and gene name
 
 gene_info_with_phenotypes<- gene_info %>%
  left_join(gene_disease_home_page, by=c("gene_symbol","disorder_id")) %>%
  distinct() %>%
  mutate(query_phenotype = ifelse(!is.na(query_phenotype),query_phenotype,"-")) %>%
  mutate(match_phenotype = ifelse(!is.na(match_phenotype),match_phenotype,"-"))

 # Summarise as you wish 
 gene_info_with_phenotypes_reframe <- gene_info_with_phenotypes %>%
  group_by(gene_symbol) %>%
  reframe(
    disease_gene,
    IMPC_pipeline,
    IMPC_phenotypes,
    PhenoDigm_match,
    max_score,
    disorder_id=paste(disorder_id, collapse="|"),
    disorder_name=paste(disorder_name, collapse="|"),
    query_phenotype=paste(query_phenotype, collapse="|"),
    match_phenotype=paste(match_phenotype, collapse="|")
  ) %>%
  distinct()

# export home page gene summary --------------------------------------------------

write_parquet(gene_info_with_phenotypes_reframe, "./data/output/home_gene_summary.parquet")

# match vs no match -------------------------------------------------------
# Genes with PhenoDigm score
genes_score = unique(c(model_disease_orphanet_score_no0$mgi_id,
                       model_disease_omim_score_no0$mgi_id))

# Number of genes with a match
# length(genes_score)
# length(unique(phenodigm_matches$mgi_id))

genes_pheno_hpo_nomatch <- disease_genes_impc_categories %>%
  filter(!mgi_id %in% genes_score) %>%
  filter(phenotype !="n")%>%
  filter(hpo_annotations == "y") %>%
  inner_join(hm_ortho_symbol)


genes_pheno_hpo_nomatch_match <- disease_genes_impc_categories %>%
  mutate(phenodigm_match = ifelse(!mgi_id %in% genes_score,"n","y"))

# Number of no match genes
log_event(genes_pheno_hp_no_match = length(unique(genes_pheno_hpo_nomatch$mgi_id)))

# lethal curation ---------------------------------------------------------
# Cache omim_curation.tsv
omim_curation <- read_delim("./data/auxiliary/omim_curation.tsv",
                            delim = "\t",col_names = TRUE)

lethal_humans_nomatch <- omim_curation %>%
  select(hgnc_id, earliest_lethality_category, omim_id) %>%
  distinct() %>%
  inner_join(genes_pheno_hpo_nomatch_match) %>%
  filter(hpo_annotations == "y") %>%
  filter(phenotype == "y") %>%
  filter(phenodigm_match == "n") %>%
  filter(viability %in% c("subviable","lethal"))

# TODO: Inlcude these tables in the logs in a readable format
table(lethal_humans_nomatch$earliest_lethality_category)

lethal_humans <- omim_curation %>%
  select(hgnc_id, earliest_lethality_category) %>%
  distinct() %>%
  inner_join(genes_pheno_hpo_nomatch)

# TODO: Inlcude these tables in the logs in a readable format
table(lethal_humans$earliest_lethality_category)

lethal_humans_check = lethal_humans %>%
  filter(earliest_lethality_category %in% c("L1","L2","L3"))

# TODO: Silenced for now since it was only used for the DM paper.
# panel <- read_delim("./data/auxiliary/panelapp_by_gene_panel_aus_data_accessed_050424.txt",
#                     delim = "\t",col_names = TRUE)
# 
# green_disease_category  <- panel %>%
#   filter(confidence_level_color == "green") %>%
#   filter(panel_disease_group !="") %>%
#   select(hgnc_id, panel_disease_group) %>%
#   distinct()
# 
# total_disease_category <-  lethal_humans_check %>%
#   inner_join(green_disease_category) %>%
#   group_by(panel_disease_group) %>%
#   tally() %>%
#   rename(n_disease_category = n )

# matches_disease_category is called, but cannot find reference to it.
# SKIPPED: since it's not needed for DM portal or benchmarking. 
# total_disease_category_match <-  matches_disease_category %>%
#   group_by(panel_disease_group, match) %>%
#   tally() %>%
#   inner_join(total_disease_category) %>%
#   mutate(percentage = (n/n_disease_category)*100) %>%
#   filter(match == "y") %>%
#   arrange(-percentage)


# PhenoDigm scores IMPC ---------------------------------------------------

# Get how many have MP annotations
# Get OMIM diseases with Hp annotations
# Check how manny genes linked to these 
# How many have a score over 40
####################################################################
model_no0 <- fread("./data/phenodigm/disease_model_association_omim_impc.tsv.gz") %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  inner_join(model, by = c("match" = "id")) %>%
  inner_join(hm_ortho_symbol) %>%
  inner_join(disease, by = c("query" = "disorder_id")) %>%
  mutate(gene_disease_pair = paste0(hgnc_id,"_", query)) 

model_no0_nodisease <- model_no0 %>%
  filter(!mgi_id %in% unique(gene_disease$mgi_id)) %>%
  arrange(-score) %>%
  select(query, match, mgi_id, hgnc_id, gene_symbol, disorder_name, score, 
         life_stage,description,
         query_phenotype, match_phenotype) %>%
  filter(score > 40)



# Genes not associated with a disease that recapitulate clinical features of disease by phenotypic sim.
model_no0_nodisease_dist <- model_no0 %>%
  filter(!mgi_id %in% unique(gene_disease$mgi_id))

# Log number of genes-predicted/other that recapitulate human disease with a PD score above 40
log_event("genes_predicted_recapitulate_human_disease_phenodigm_score>40" = length(unique(model_no0_nodisease_dist$mgi_id)))


# arrange(-score) %>%
# select(query, match, mgi_id, hgnc_id, disorder_name, score, 
#        life_stage,description,
#        query_phenotype, match_phenotype) %>%
# filter(score > 40)


# Non-disease gene(predicted/other) summary stats:
summary_stats_predicted <- model_no0_nodisease_dist %>%
        summarise(
          min = min(score, na.rm = TRUE),
          q1 = quantile(score, 0.25, na.rm = TRUE),
          median = median(score, na.rm = TRUE),
          mean = mean(score, na.rm = TRUE),
          q3 = quantile(score, 0.75, na.rm = TRUE),
          max = max(score, na.rm = TRUE),
          iqr = IQR(score, na.rm = TRUE)
        )

# Log disease-predicted gene stats
log_event(score_disease_predicted_min = summary_stats_predicted$min)
log_event(score_disease_predicted_q1 = summary_stats_predicted$q1)
log_event(score_disease_predicted_median = summary_stats_predicted$median)
log_event(score_disease_predicted_mean = summary_stats_predicted$mean)
log_event(score_disease_predicted_q3 = summary_stats_predicted$q3)
log_event(score_disease_predicted_max = summary_stats_predicted$max)
log_event(score_disease_predicted_iqr = summary_stats_predicted$iqr)

# Violin plot for visualisation and cutoff
violin_score_predicted <- ggplot(model_no0_nodisease_dist, aes(x = 0, y = score)) +
  geom_violin(fill = "#00AFBB", trim = FALSE, width=0.3, adjust = 0.3) +
  geom_boxplot(width = 0.02, fill = "white", outlier.shape = NA) +
  stat_summary(
    aes(shape = "Mean ± SD"),
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", color = "black"
  ) +
  geom_hline(
    data = data.frame(yintercept = 40),
    aes(yintercept = yintercept, , linetype = "Score threshold (40)"),
    color = "purple", linewidth = 0.8
  )  +
  scale_shape_manual(values = c("Mean ± SD" = 16)) +
  scale_linetype_manual(values = c("Score threshold (40)" = "dashed")) +
  guides(
    shape = guide_legend(title = NULL),
    linetype = guide_legend(title = NULL, override.aes = list(color = "purple"))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 0),
    legend.position = c(0.9,0.9)
  ) +
  labs(y = "PhenoDigm Score") +
  labs(x = "PhenoDigm other hits")

violin_score_predicted

write_parquet(model_no0_nodisease,
              "./data/output/phenodigm_other.parquet", compression="zstd", compression_level=5)

# PhenoDigm other is getting quite large, we can split in multiple files of 35MB to comply with GitHub enterprise's limits.
# 1. Calculate the the number of rows per chunk needed to fit the required file size.
target_mb      <- 90
file_size_mb   <- as.numeric(object.size(model_no0_nodisease)) / (1024^2)
total_rows     <- nrow(model_no0_nodisease)
rows_per_mb    <- total_rows / file_size_mb
rows_per_chunk <- floor(rows_per_mb * target_mb)

# 2. Write to a new dir, in parquet format. Overwrites if existing. 
write_dataset(
  dataset = model_no0_nodisease,
  path = "./data/output/phenodigm_other",
  format = "parquet",
  max_rows_per_file = rows_per_chunk,
  existing_data_behavior = "overwrite",
  basename_template = "part-{i}-phenodigm_other.parquet"
)

###########################################################################
#### impc vs non impc match
############################################################################

impc_match = genes_pheno_hpo_nomatch_match %>%
  filter(phenodigm_match == "y")


model_nonimpc_disease_omim_score_no0 <- open_dataset("./data/phenodigm/disease_model_association_omim_nonimpc.tsv.gz", 
                                                   format = "tsv", col_names = TRUE) %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  inner_join(model, by = c("match" = "id"),relationship = "many-to-many") %>%
  #filter(match =="MGI") %>%
  filter(mgi_id %in% all_disease$mgi_id)%>%
  filter(mgi_id %in% genes_mouse_human_phenotypes$mgi_id)%>%
  inner_join(hm_ortho_symbol) %>%
  mutate(gene_disease = paste0(hgnc_id,"_", query)) %>%
  inner_join(disease_genes_impc_pairs) %>%
  distinct() %>%
  collect()


model_nonimpc_disease_orphanet_score_no0 <- open_dataset("./data/phenodigm/disease_model_association_orphanet_nonimpc.tsv.gz", 
                                                       format= "tsv", 
                                                       col_names = TRUE) %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  inner_join(model, by = c("match" = "id"),relationship = "many-to-many") %>%
  #filter(match =="MGI") %>%
  filter(mgi_id %in% all_disease$mgi_id)%>%
  filter(mgi_id %in% genes_mouse_human_phenotypes$mgi_id)%>%
  inner_join(hm_ortho_symbol) %>%
  mutate(gene_disease = paste0(hgnc_id,"_", query)) %>%
  inner_join(disease_genes_impc_pairs) %>%
  distinct() %>%
  collect()


phenodigm_matches_mgi <- model_nonimpc_disease_omim_score_no0  %>%
  bind_rows(model_nonimpc_disease_orphanet_score_no0)

b = head(phenodigm_matches_mgi)

# comparison --------------------------------------------------------------

phenodigm_matches_impc_summary <- phenodigm_matches_impc  %>%
  group_by(gene_disease) %>%
  mutate(max_score_impc = max(score)) %>%
  select(gene_disease, max_score_impc,gene_symbol,mgi_id) %>%
  distinct()

phenodigm_matches_mgi_summary <- phenodigm_matches_mgi  %>%
  group_by(gene_disease) %>%
  mutate(max_score_mgi = max(score)) %>%
  select(gene_disease, max_score_mgi) %>%
  distinct()

max_score_comparison <- phenodigm_matches_impc_summary %>%
  full_join(phenodigm_matches_mgi_summary,
            by=c("gene_disease" = "gene_disease")) %>%
  filter(!is.na(max_score_impc)) %>%
  mutate(comparison = ifelse(is.na(max_score_mgi),"impc_novel",
                             ifelse(max_score_impc > max_score_mgi,"impc_higher","-"))) 

## ones that are not higher by gene but by gene-disease


phenodigm_matches_impc_summary_gene <- phenodigm_matches_impc  %>%
  group_by(gene_symbol) %>%
  mutate(max_score_gene_impc = max(score)) %>%
  select(gene_symbol, max_score_gene_impc) %>%
  distinct()

phenodigm_matches_mgi_summary_gene <- phenodigm_matches_mgi  %>%
  group_by(gene_symbol) %>%
  mutate(max_score_gene_mgi = max(score)) %>%
  select(gene_symbol,max_score_gene_mgi) %>%
  distinct()

max_score_gene_comparison <- phenodigm_matches_impc_summary_gene %>%
  full_join(phenodigm_matches_mgi_summary_gene) %>%
  filter(!is.na(max_score_gene_impc)) %>%
  mutate(comparison_gene = ifelse(is.na(max_score_gene_mgi),"impc_gene_novel",
                                  ifelse(max_score_gene_impc > max_score_gene_mgi,"impc_gene_higher","-")))



### merge both

max_score_comparison %>%
  inner_join(max_score_gene_comparison,by = c("gene_symbol")) %>%
  mutate(novel_assoc_notgene = ifelse(
    comparison_gene !="impc_gene_novel" &
      comparison == "impc_novel", "impc_gene_disease_novel","-")) -> 
  max_score_comparison_final

## Write table to check max score and new candidate models
# write.table(max_score_comparison_final,
#             "./data/output/max_score_comparison_final.tsv",
#             quote = F, sep = "\t", row.names = F)


# These tables are not entirely necessary but they are a nice to have for summary stats. Not sure what the last table is. 
table(max_score_comparison$comparison)
table(max_score_gene_comparison$comparison_gene)
table(max_score_comparison_final$novel_assoc_notgene)

#  Number of genes in the categories of columns: comparison, comparison_gene, novel_association_notgene
# Generate counts for all categories in each column
counts_comparison <- table(max_score_comparison_final$comparison)
counts_comparison_gene <- table(max_score_comparison_final$comparison_gene)
counts_novel_assoc_notgene <- table(max_score_comparison_final$novel_assoc_notgene)

# Access pattern for those counts
# Comparison 
counts_comparison["impc_higher"]
counts_comparison["impc_novel"]
counts_comparison["-"]

# Comparison gene
counts_comparison_gene["-"]
counts_comparison_gene["impc_gene_higher"]
counts_comparison_gene["impc_gene_novel"]

# Novel association not gene
counts_novel_assoc_notgene["-"]
counts_novel_assoc_notgene["impc_gene_disease_novel"]

# Log stats of MGI comparison
  # Comparison 
log_event(mgi_comparison_null = counts_comparison[["-"]])
log_event(mgi_comparison_impc_higher = counts_comparison[["impc_higher"]])
log_event(mgi_comparison_impc_novel = counts_comparison[["impc_novel"]])
  
  # Comparison gene
log_event(mgi_comparison_gene_null = counts_comparison_gene[["-"]])
log_event(mgi_comparison_gene_impc_gene_higher = counts_comparison_gene[["impc_gene_higher"]])
log_event(mgi_comparison_gene_impc_gene_novel = counts_comparison_gene[["impc_gene_novel"]])
  
  # Novel association not gene
log_event(mgi_novel_assoc_not_gene_null = counts_novel_assoc_notgene[["-"]])
log_event(mgi_novel_assoc_not_gene_impc_gene_disease_novel = counts_novel_assoc_notgene[["impc_gene_disease_novel"]])

# table(max_score_comparison_final)

# Close logger
log_close()


# Write the json info in the logger as a df for easy manipulation
# Read file
lines <- readLines("./data/output/log/post_processing_results.log", warn = FALSE)

# Keep only lines starting with {
json_lines <- lines[startsWith(trimws(lines), "{")]

# Convert to long format with keys as rows
stats_df <- lapply(seq_along(json_lines), function(i) {
  json_obj <- fromJSON(json_lines[i])
  data.frame(
    line_num = i,
    metric = names(json_obj),
    value = as.character(unlist(json_obj)),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write_parquet(stats_df,"./data/output/log/post_processing_results.parquet")

##########################  Pheval-impc Benchmarking ####################################################
# Finding all OMIM genes with an IMPC gene
# Filtered version of stats_all file 
tested_impc_genes <- stats_all %>%
  dplyr::rename(mgi_id = marker_accession_id) %>%
  select(mgi_id) %>%
  distinct() %>%
  replace(is.na(.),"-") %>%
  filter(mgi_id !="-") %>%
  filter(!is.na(mgi_id))

# TODO: is this the same as disease_genes_impc_categories?
# This gets all omim ids with impc genes
omim_diseases_with_impc_genes <- gene_disease %>%
  filter(grepl("OMIM",disorder_id)) %>%
  filter(mgi_id %in% tested_impc_genes$mgi_id) %>%
  select(disorder_id) %>%
  distinct()

#### Produce benchmarking files ####
# Should be all disease ids with curated disease genes, keeping the highest score per association.
# For example:
# omim_id_1 - hgnc:1 - score
# omim_id_1 - hgnc:2 - score
# omim_id_2 - hgnc:5 - score


# Read disease_model_association files
# 1. omim-all models
# 2. omim-impc models-only
dma_omim <- fread("./data/phenodigm/disease_model_association_omim.tsv.gz")
dma_impc_omim <- fread("./data/phenodigm/disease_model_association_omim_impc.tsv.gz")


# 1  Select scores for all models
# Set phenodigm filter thresholds used in the pipeline
phenodigm_min_perc <- 60
phenodigm_min_raw <- 1.3

# TODO: this should be a function. 
# Group by genes, and keep the highest scoring row for each gene
top_genes_all_models <- dma_omim %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  filter(score_avg_norm > phenodigm_min_perc | score_max_norm > phenodigm_min_perc | score_avg_raw > phenodigm_min_raw | score_max_raw > phenodigm_min_raw) %>%
  inner_join(model, by = c("match" = "id"),relationship = "many-to-many") %>%
  inner_join(hm_ortho_symbol) %>%
  filter(query %in% omim_diseases_with_impc_genes$disorder_id) %>%
  group_by(query) %>%
  arrange(desc(score)) %>%
  distinct(query, hgnc_id, .keep_all = TRUE)


top_genes_all_models_tidy <- top_genes_all_models %>%
  dplyr::rename(disorder_id = query) %>%
  select(disorder_id,
         gene_symbol, hgnc_id,
         score, source) %>%
  distinct()


# 2 Same but for IMPC models only
top_genes_impc_models <- dma_impc_omim %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  filter(score_avg_norm > phenodigm_min_perc | score_max_norm > phenodigm_min_perc | score_avg_raw > phenodigm_min_raw | score_max_raw > phenodigm_min_raw) %>%
  inner_join(model, by = c("match" = "id"),relationship = "many-to-many") %>%
  inner_join(hm_ortho_symbol) %>%
  filter(query %in% omim_diseases_with_impc_genes$disorder_id) %>%
  group_by(query) %>%
  arrange(desc(score)) %>%
  distinct(query, hgnc_id, .keep_all = TRUE)

top_genes_impc_models_tidy <- top_genes_impc_models %>%
  dplyr::rename(disorder_id = query) %>%
  select(disorder_id,
         gene_symbol, hgnc_id,
         score, source) %>%
  distinct()


#########################################################################################################

# Write benchmarking files to be read by pheval_impc_phenodigm

# Write all models file
write_parquet(top_genes_all_models_tidy,
              "./data/output/phenodigm_scores_benchmarking_all_models.parquet",
              compression="zstd", compression_level=5)

# Write impc models only file
write_parquet(top_genes_impc_models_tidy,
              "./data/output/phenodigm_scores_benchmarking_impc_models.parquet",
              compression="zstd", compression_level=5)
