"""Back-compat entry point for PhenoDigm2.

The implementation now lives in ``pd2.cli``. Prefer the installed console
command ``phenodigm2`` (or ``phenodigm``); this shim keeps
``python3 phenodigm2.py ...`` working.

@author: Tomasz Konopka
"""

from pd2.cli import main


if __name__ == "__main__":
    main()
