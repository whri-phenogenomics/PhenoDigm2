""" Runs post processing of Phenodigm release.
    This involves prodcing necessary files for data analysis of each DR and the Disease Models ShinyApp
    
    By: Diego Pava"""


from . import tools as pd2tools
import luigi
from luigi import Task, LocalTarget
import os
from pathlib import Path
import requests
from typing import Dict

# Contants with paths to download data
def post_process_paths(config) -> Dict[str, Path]:

    post_proc_dir = _get_post_proc_dir(config)
    return {
        "data": Path(post_proc_dir, "data"),
        "data_aux": Path(post_proc_dir, "data", "auxiliary"),
        "hpo": Path(post_proc_dir, "data", "hpo"),
        "impc": Path(post_proc_dir ,"data", "impc"),
        "intermediate": Path(post_proc_dir ,"data", "intermediate"),
        "output": Path(post_proc_dir,"data", "output"),
        "phenodigm": Path(post_proc_dir,"data", "phenodigm"),
        "scripts": Path(post_proc_dir,"scripts"),
        "script_aux": Path(post_proc_dir,"scripts", "auxiliary")
    }

# helper functions here
def _get_post_proc_dir(config) -> Path:
    rootdir = Path(config.db).resolve()
    return Path(rootdir, "post_processing")

def _get_pipeline_status_dir(config) -> None:
    post_proc_dir = _get_post_proc_dir(config)
    return Path(post_proc_dir, "pipeline_status")

def run_post_process(config):
    pd2tools.log("Running post processing pipeline")
    
    #LOGIC HERE
    # post_proc_dir = _get_post_proc_dir(config)
    # Or use default (current directory)
    luigi.build([CreatePostProcessDirs(config=config)], local_scheduler=True)

def download_data(url, filename):
    """Generic function to download data passing a URL.

    Args:
        url (str): URL of the ontology/gwas catalogue
        filename (str): Filename to save the data. Can pass a path.
    """
    print(f"Fetching {filename}...")
    response = requests.get(url, timeout=10)
    with open(filename, "wb") as f:
        f.write(response.content)
    print("Done")


# luigi tasks here
class CreatePostProcessDirs(Task):

    config = luigi.Parameter()

    def output(self):
        post_proc_dir = _get_post_proc_dir(self.config)
        return LocalTarget(Path(post_proc_dir, "pipeline_status",".directories_created"))

    def run(self):

        paths = post_process_paths(self.config)
        directories = [v for _, v in paths.items()]

        try:
            # Create all project directories
            for dir_path in directories:
                os.makedirs(dir_path, exist_ok=True)
                # print(f"making_dir_path:{dir_path}")
                pd2tools.log(f"Created directory: {dir_path}")

            # Create the marker file to indicate successful completion
            with self.output().open("w") as f:
                f.write("Directories created successfully.")
            pd2tools.log("Created post-processing directories successfully.")

        except Exception as e:
            pd2tools.log(f"Error creating directories: {e}")
            raise e
        
        # Create a pipeline_status dir
        os.makedirs(_get_pipeline_status_dir(self.config), exist_ok=True)


class DownloadResources(Task):
    def requires(self):
        return CreatePostProcessDirs()
    
    def output(self):
        post_proc_dir = _get_post_proc_dir(self.config)
        return LocalTarget(Path(post_proc_dir, "pipeline_status",".files_downloaded"))

    def run(self):
        pass
        # TODO: Files to download
        # 1. GENTAR Orthologues: ### https://www.gentar.org/orthology-api/api/ortholog/one_to_one/impc/write_to_tsv_file
        # 2. IMPC DR geno-pheno assertions all http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/genotype-phenotype-assertions-ALL.csv.gz
        # 3. IMPC viability: http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz
        # 4. IMPC statistical results: http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/statistical-results-ALL.csv.gz
        # 5. HPO genes_to_phenotype: https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-08-13/genes_to_phenotype.txt
        # 6. 

        # Auxiliary static
        # 1. "./data/auxiliary/omim_curation.tsv"
        
        # NOTE: For now we will download everything here because the versions of the files are to remain flexible and cannot be determined at the begining of the build
