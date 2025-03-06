### IMPC DR21 models and PhenoDigm scores
### It requires:
### 1) tables from PhenoDigm 22 sqlite database in Apocrita
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
packages <- c("R.utils", "dplyr", "readr", "data.table", "arrow")

for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

library("dplyr")
library("readr")
library(data.table)
library(arrow)

# Set the root directory dynamically
setwd(post_processing_path)

# Execute external Rscript
source("./scripts/auxiliary/hgnc_symbol_checker.R")


# import data -------------------------------------------------------------

model_genotype <- fread("./data/phenodigm/model_genotype.tsv.gz") %>%
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


length(unique(gene_disease$gene_disease))
length(unique(gene_disease$hgnc_id))

# import impc data: phenotype associations, viability, all -----------------
# "Cache" assertions
geno_pheno_assertions <- fread("./data/impc/genotype-phenotype-assertions-ALL.csv.gz")

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

all <- stats_all %>%
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

all_disease <- all %>%
  filter(mgi_id %in% unique(gene_disease$mgi_id)) %>%
  mutate(phenotype = ifelse(top_level_mp == "-" & viability =="-","n",
                            ifelse(top_level_mp == "-" & viability =="viable","n","y"))) %>%
  rename(top_level_mp_term_name = top_level_mp)

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


disease_genes_impc_categories <- all_disease_viability_group %>%
  left_join(hpo_genes) %>%
  replace(is.na(.),"-")

head(disease_genes_impc_categories)
dim(disease_genes_impc_categories)

## 2848 with a human orthologue and disease association

disease_genes_impc_pairs <- gene_disease %>%
  filter(mgi_id %in% disease_genes_impc_categories$mgi_id)

genes_mouse_phenotypes <- disease_genes_impc_categories %>%
  filter(phenotype == "y")

dim(genes_mouse_phenotypes)

## 2621 with IMPC mouse phenotypes encoded as MP terms

genes_human_phenotypes <- disease_genes_impc_categories %>%
  filter(hpo_annotations == "y")

dim(genes_human_phenotypes)

genes_mouse_human_phenotypes <- genes_mouse_phenotypes %>%
  filter(hpo_annotations == "y")

dim(genes_mouse_human_phenotypes)

## 2504  with human phenotypes encoded as HPO terms

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

phenodigm_matches_impc = phenodigm_matches

matches_tidy <- phenodigm_matches_impc %>%
  select(disorder_id, disorder_name,
         gene_symbol, description,
         score, query_phenotype, match_phenotype)


# export files for shiny app ----------------------------------------------


write.table(phenodigm_matches_impc,
            "./data/output/phenodigm_matches_DR22.txt",
            quote = F, sep = "\t", row.names = F)

write.table(matches_tidy,
            "./data/output/phenodigm_matches_tidy_DR22.txt",
            quote = F, sep = "\t", row.names = F)

# gene summary file -------------------------------------------------------

gene_disease_by_gene <- gene_disease %>%
  group_by(hgnc_id) %>%
  summarise(disorder_id = paste0(unique(disorder_id), collapse = "|"),
            disorder_name = paste0(unique(disorder_name), collapse = "|"))

max_score <- phenodigm_matches %>%
  group_by(hgnc_id) %>%
  summarise(max_score = max(score))

gene_summary_df <- hm_ortho_symbol %>%
  left_join(gene_disease_by_gene) %>%
  mutate(IMPC_pipeline = ifelse(mgi_id %in% all$mgi_id,"yes","no")) %>%
  mutate(IMPC_phenotypes = ifelse(mgi_id %in% genes_mouse_phenotypes$mgi_id,"yes","no") )%>%
  mutate(HPO_phenotypes = ifelse(mgi_id %in% genes_human_phenotypes$mgi_id,"yes","no") )%>%
  mutate(PhenoDigm_match = ifelse (mgi_id %in% unique(phenodigm_matches$mgi_id),
                                   "yes","no")) %>%
  replace(is.na(.),"-")  %>%
  left_join(max_score) 

filter_matches <- gene_summary_df %>%
  filter(PhenoDigm_match == "yes")


# export gene summary file ------------------------------------------------


write.table(gene_summary_df ,
            "./data/output/gene_summary_DR22.txt",
            quote = F, sep = "\t", row.names = F)


# match vs no match -------------------------------------------------------

genes_score = unique(c(model_disease_orphanet_score_no0$mgi_id,
                       model_disease_omim_score_no0$mgi_id))

length(genes_score)
length(unique(phenodigm_matches$mgi_id))

genes_pheno_hpo_nomatch <- disease_genes_impc_categories %>%
  filter(!mgi_id %in% genes_score) %>%
  filter(phenotype !="n")%>%
  filter(hpo_annotations == "y") %>%
  inner_join(hm_ortho_symbol)


genes_pheno_hpo_nomatch_match <- disease_genes_impc_categories %>%
  mutate(phenodigm_match = ifelse(!mgi_id %in% genes_score,"n","y"))

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

table(lethal_humans_nomatch$earliest_lethality_category)

lethal_humans <- omim_curation %>%
  select(hgnc_id, earliest_lethality_category) %>%
  distinct() %>%
  inner_join(genes_pheno_hpo_nomatch)

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
  select(query, match, mgi_id, hgnc_id, disorder_name, score, 
         life_stage,description,
         query_phenotype, match_phenotype) %>%
  filter(score > 40)

write.table(model_no0_nodisease,
            "./data/output/model_no0_nodisease_DR22.txt",
            quote = F, sep = "\t", row.names = F)

###########################################################################
#### impc vs non impc match
############################################################################

impc_match = genes_pheno_hpo_nomatch_match %>%
  filter(phenodigm_match == "y")

model_nonimpc_disease_omim_score_no0 <- read_delim("./data/phenodigm/disease_model_association_omim_nonimpc.tsv.gz", 
                                                   delim = "\t", col_names = TRUE) %>%
  mutate(score =(score_avg_norm + score_max_norm)/2) %>%
  filter(score > 0) %>%
  inner_join(model, by = c("match" = "id"),relationship = "many-to-many") %>%
  #filter(match =="MGI") %>%
  filter(mgi_id %in% all_disease$mgi_id)%>%
  filter(mgi_id %in% genes_mouse_human_phenotypes$mgi_id)%>%
  inner_join(hm_ortho_symbol) %>%
  mutate(gene_disease = paste0(hgnc_id,"_", query)) %>%
  inner_join(disease_genes_impc_pairs) %>%
  distinct()


model_nonimpc_disease_orphanet_score_no0 <- read_delim("./data/phenodigm/disease_model_association_orphanet_nonimpc.tsv.gz", 
                                                       delim = "\t", 
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
  distinct()


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


# These tables are not entirely necessary but they are a nice to have for summary stats. Not sure what the last table is. 
table(max_score_comparison$comparison)
table(max_score_gene_comparison$comparison_gene)
# table(all_score_comparison$novel_assoc_notgene)


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
write.table(top_genes_all_models_tidy,
            "./data/output/phenodigm_scores_benchmarking_DR22_all_models.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Write impc models only file
write.table(top_genes_impc_models_tidy,
            "./data/output/phenodigm_scores_benchmarking_DR22_impc_models.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

