#!/bin/bash

# Function to display usage help
usage() {
    echo "Usage: $0 --env <virtuan_env_path> --db <database_path> --code <path_to_phenodigm_repo> --output_dir <output_directory>"
    exit 1  
}

# Parse flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
    --env) VENV_PATH="$2"; shift ;;
    --db) DB_PATH="$2"; shift ;;
    --code) CODE_PATH="$2"; shift ;;
    --output_dir) OUTPUT_DIR="$2"; shift ;;
    *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Check for mandatory arguments
if [[ -z "$VENV_PATH" || -z "$DB_PATH" || -z "$CODE_PATH" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# TODO: after this nothgin will work because the rest is designed for decompressed variable. 
# Activate pigz for compression
module load pigz/2.6 
# 1. Decompress the db if codeits compressed
# Conditional to check 
if [[ "$DB_PATH" == *.tar.gz ]]; then 
    echo "Decompressing files..."
    tar --use-compress-program=pigz -xvf "$DB_PATH"
    DB_PATH="${DB_PATH%.tar.gz}"
else
    echo "No decompression needed, the file does not have a .tar.gz extension."
fi

# 2. Activate the virutal environment
echo "Activating virtual environment..."
module load python/3.8.5
source "$VENV_PATH/bin/activate" || { echo "Failed to activate virtual environment"; exit 1; }

#3. Create output directory:
mkdir -p "$OUTPUT_DIR"

# 4. a. Extract all of the required tables
echo "Extracting tables..."
tables=("model" "model_genotype" "disease" "disease_gene_mapping" "gene_gene_mapping")
for i in "${!tables[@]}"; do
    table="${tables[$i]}"
    if python3 "$CODE_PATH/phenodigm2.py" export --db "$DB_PATH" --table "$table" | pigz -k -p4 > "$OUTPUT_DIR/${table}.tsv.gz"; then
        echo "Processed table: $table, output saved to: $OUTPUT_DIR/${table}.tsv.gz"
    else
        echo "Error processing table: $table" >&2
    fi
done


# 4 b. Extract the conditional disease_model_associations
echo "Extracting conditional tables..."
conditions=("query LIKE '%OMIM%' AND match LIKE '%#%'" "query LIKE '%ORPHA%' AND match LIKE '%#%'" "query LIKE '%DECIPHER%' AND match LIKE '%#%'" "query LIKE '%OMIM%' AND match NOT LIKE '%#%'" "query LIKE '%ORPHA%' AND match NOT LIKE '%#%'" "query LIKE '%DECIPHER%' AND match NOT LIKE '%#%'")
filenames=("disease_model_association_omim_impc" "disease_model_association_orphanet_impc" "disease_model_association_decipher_impc" "disease_model_association_omim_nonimpc" "disease_model_association_orphanet_nonimpc" "disease_model_association_decipher_nonimpc")

for i in "${!conditions[@]}"; do
    condition="${conditions[$i]}"
    file="${filenames[$i]}"

    if python3 "$CODE_PATH/phenodigm2.py" export --db "$DB_PATH" --table disease_model_association --where "${condition}" | pigz -k -p4 > "$OUTPUT_DIR"/${file}.tsv.gz; then
        echo "Processed file: $file, output saved to: $OUTPUT_DIR/${file}.tsv.gz"
    else
        echo "Error processing table: $file" >&2
    fi
done
