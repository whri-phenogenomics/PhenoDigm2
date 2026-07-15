"""
Helper module to load and parse a config.yaml file to locate external resources.
"""

import yaml
from pathlib import Path
from typing import Optional
from pydantic import BaseModel


# pydantic model to optionally override the bundled post-processing inputs.
# The omim curation file and both R scripts ship inside the package
# (pd2/resources, pd2/rscripts) and are used by default. Each optional path
# overrides a bundled copy, e.g. to supply a freshly curated omim file or run a
# locally edited R script for a different R version.
class PostProcessConfig(BaseModel):
    omim_curation_path: Optional[Path] = None
    main_r_script_path: Optional[Path] = None
    hgnc_symbol_checker_script_path: Optional[Path] = None


# helper function
def load_config(config_path: str) -> PostProcessConfig:
    with open(config_path, "r") as f:
        # An empty or fully-commented yaml parses to None -> use all defaults.
        config_dict = yaml.safe_load(f) or {}
    return PostProcessConfig(**config_dict)
