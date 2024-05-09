# Notes on changes 

i.e. discrepancies/changes from original PhenoDigm build to this
implementation (PhenoDigm2)

## XML files

Changed catalog-v001.xml into phenodigm2_catalog.xml. 
Changed uri mappings to fit new downloads folder structure
Changed these two lines 

```
<uri id="User Entered Import Resolution" 
name="http://purl.obolibrary.org/obo/upheno/zp/zp-equivalence-axioms-subq-ubr.owl" 
uri="obo/upheno/zp.owl"/>
<uri id="User Entered Import Resolution" 
name="http://purl.obolibrary.org/obo/wbphenotype/wbphenotype-equivalence-axioms-subq-ubr.owl" 
uri="obo/wbphenotype/wbphenotype-equivalence-axioms-subq-ubr.owl"/>
```

In the first line, changed target uri. That line then becomes a repeat and thus I removed it. 
Removed the second line because it was a repeat of an earlier line.



## Databases

MySQL -> SQLite

Table definitions omit primary key columns `id`.

Disease table omits column `disease_classes`. This does not seem to be
used in the original repo.

Tables with ontology (mp, hp) just use `id` rather than `mp_id` 
and `hp_id`. This allows using a single interface to manipulate both 
tables. 

Tables with ontologies (mp, hp) combined into a single table `ontology`.
Same with other tables like `hp_synonyms`.

Tables with mappings (e.g. hp_hp_mapping) just use `id` rather than
`hp_id` and just `id_hit` rather than `hp_id_hit`. This allows using 
a single interface to manipulate several tables with mappings. 
Tables with several types of mappings combined into a single 
`ontology_ontology_mapping` table.

Tables with mappings omit `term` columns. These can be added computationally
from tables like `hp_term`.  

In `mouse_model`, changed column `hom_het` to `zygosity`

Table `disease_hp` changed to `disease_phenotype`. In this table, column
`hp_id` changed to `phenotype_id`.

Table `mouse_model_mp` changed to `mouse_model_phenotype` and inner column
`mp_id` changed to `phenotype_id`.


Changed multiple table names to remove the word `mouse`. 
Eg. `mouse_model_gene_orthology` changed to `model_gene_ortholog`;
 `mouse_gene_summary` changed to `model_gene_summary`;
 `mouse_mouse_association` changed to `model_model_association`.

Removed `best_impc` tables. 

In `model_model_association` changed to `score`.

Merged tables `mouse_mouse_association` and `mouse_mouse_association_details`
into table `model_model_association`.

Changed column names in all association tables to have the same
names, i.e. `id`, `match` and `score`.



## Downloads

Zebrafish downloads are omitted

Added download from orphadata.org for `en_produce1.xml`. This replaces 
a missing data source `orphanet_disease_mapping` 


## Scripts

downloads - OK (fish skipped)
run_biomart_queries - OK
create db - OK
produce_orthologs - OK (using biomart tables only)
load_phenotype_ontologies -OK (hp, mp) 
processOrphanetOmim - partial (?) missing mouse_disease_gene_summary
parse_mgi_data - skipped? (impc download incorporates this, no?)
parse_impc_data - OK?
import_omim_disease_2_MGI_models - 
produce_fish-tables - skipped
flagLocusCandidates - skipped (missing data? new tables?)
loadOntologyMappings - 


produce_monarch_files.pl - incorporated into load_data.py



## Data

Zebrafish data are omitted

Ortholog mapping used to come from `HOM_MouseHumanSequence.rpt`, but 
this table has many-to-many mappings for each ortholog_id. The 
alternative file `human_mouse_mapping.txt` has similar mappings. Here
just use the second data source.

While loading phenotype data, original script had a procedure to check 
that mouse genes are included into the `mouse_gene_ortholog` table. 
The procedure then added new genes without links to orthologs. I am 
skipping this step - i.e. I do not modify `mouse_gene_ortholog` after
the initial load.

Skipped population of `mouse_model_gene_ortholog` - is this necessary?
The original script had INSERT statements in parse_impc_data.pl.

When populating disease_hp, original script skipped OMIM disease genes?
I am just transfering all phenotype_annotation.tab into the db. 

For model tables, new tables use model_id for the form `MODEL:xxx`.


## Questions

For section: flagLocusCandidates.pl
Missing input file: human_band_mapping.txt


