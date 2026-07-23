"""phenodigm2 package.

Public API:
    from phenodigm2 import PhenoDigm
"""

from phenodigm2.api import PhenoDigm
from phenodigm2.cli import main

__all__ = ["PhenoDigm", "main"]
