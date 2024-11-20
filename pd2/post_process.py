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
        "impc": Path(post_proc_dir, "data", "impc"),
        "intermediate": Path(post_proc_dir, "data", "intermediate"),
        "output": Path(post_proc_dir, "data", "output"),
        "phenodigm": Path(post_proc_dir, "data", "phenodigm"),
        "scripts": Path(post_proc_dir, "scripts"),
        "script_aux": Path(post_proc_dir, "scripts", "auxiliary")
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
    luigi.build([DownloadResources(config=config)], local_scheduler=True)

def download_data(url, filename):
    """Generic function to download data passing a URL.

    Args:
        url (str): URL of the ontology/gwas catalogue
        filename (str): Filename to save the data. Can pass a path.
    """
    pd2tools.log(f"Fetching {filename}")
    response = requests.get(url, timeout=10)
    with open(filename, "wb") as f:
        f.write(response.content)

def downloads_dict(impc_data_release: str = 'latest'):
    """Generic dictionary containing urls, filenames and paths to download post_processing files

    Args:
        impc_data_release (str, optional): IMPC data release to fetch. For older DRs type in the format: 'release-21.1' Defaults to 'latest'.

    Returns:
        Dict[str]: Returns a dictionary with the urls, filenames and paths to download post_processing files
    """
    downloads = {
            "orthology": {
                "url": "https://www.gentar.org/orthology-api/api/ortholog/one_to_one/impc/write_to_tsv_file",
                "filename": "one_to_one_orthologs.tsv",
                "targetdir": "data_aux"
            },
            "impc_geno_pheno_assertions":{
                "url": f"http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/{impc_data_release}/results/genotype-phenotype-assertions-ALL.csv.gz",
                "filename": "genotype-phenotype-assertions-ALL.csv.gz",
                "targetdir": "impc"
            },
            "impc_viability": {
                "url": f"http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/{impc_data_release}/results/viability.csv.gz",
                "filename": "viability.csv.gz",
                "targetdir": "impc",
            },
            "impc_stat_result": {
                "url": f"http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/{impc_data_release}/results/statistical-results-ALL.csv.gz",
                "filename": "statistical-results-ALL.csv.gz",
                "targetdir": "impc"
            },
            "hpo_gene_to_pheno": {
                "url": "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2022-02-14/genes_to_phenotype.txt",
                "filename": "genes_to_phenotype.txt",
                "targetdir": "hpo"
            }
        }
    return downloads


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

    config = luigi.Parameter()

    def requires(self):
        return CreatePostProcessDirs(config=self.config)
    
    def output(self):
        post_proc_dir = _get_post_proc_dir(self.config)
        return LocalTarget(Path(post_proc_dir, "pipeline_status",".files_downloaded"))

    def run(self):
        # Download data
        pd2tools.log("Downloading post_process resources...")
        downloads = downloads_dict(impc_data_release='release-14.0')
        download_paths = post_process_paths(self.config)
        for resource, values in downloads.items():
            url = values["url"]
            filename = values["filename"]
            targetdir = values["targetdir"]
            target_path = download_paths[targetdir]
            file_path = Path(target_path, filename)
            
            download_data(url, file_path)
        pd2tools.log(f"Now you should add the omim_curation.tsv file to {download_paths['data_aux']}")

        # Create the marker file to indicate successful completion
        with self.output().open("w") as f:
            f.write("Files downloaded successfully.")
        pd2tools.log("Downloaded post-processing files successfully.")

        # NOTE: For now we will download everything here because the versions of the files are to remain flexible and cannot be determined at the begining of the build
