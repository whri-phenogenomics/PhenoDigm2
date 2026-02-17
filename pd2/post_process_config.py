"""
Helper module to load and parse a config.yaml file to locate external resources.
"""

import yaml
from pathlib import Path
from pydantic import BaseModel


# pydantic model to take external paths to locate external resources
class PostProcessConfig(BaseModel):
    omim_curation_path: Path
    main_r_script_path: Path
    hgnc_symbol_checker_script_path: Path


# helper function
def load_config(config_path: str) -> PostProcessConfig:
    with open(config_path, "r") as f:
        config_dict = yaml.safe_load(f)
    return PostProcessConfig(**config_dict)
