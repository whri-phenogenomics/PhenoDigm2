"""Example Airflow DAG driving the PhenoDigm2 pipeline via the Python API.

Each pipeline stage is a ``PythonOperator`` whose callable invokes a
:class:`phenodigm2.PhenoDigm` method. ``PhenoDigm(DB)`` is constructed *inside*
each task callable on purpose: Airflow runs tasks in separate worker processes,
and the build directory on disk is the shared state, so the object itself is
stateless and cheap to recreate per task.

Airflow is intentionally NOT a runtime dependency of ``phenodigm2``. Install it
in your orchestration environment alongside the package, e.g.::

    pip install apache-airflow phenodigm2

Configure the build directory with the ``PHENODIGM2_DB`` environment variable
(or edit ``DB`` below). If ``download`` needs the OMIM feed, also set
``OMIM_API_KEY`` in the worker environment.
"""

from __future__ import annotations

import os
from datetime import datetime

from airflow import DAG
from airflow.operators.python import PythonOperator

from phenodigm2 import PhenoDigm

DB = os.environ.get("PHENODIGM2_DB", "/data/phenodigm2_build")

# Stages to run, in dependency order.
STAGES = [
    "download",
    "build",
    "ontology_mapping",
    "score",
    "index",
    "solr_prepare",
    "solr",
]


def _stage(name: str, **overrides):
    """Return a zero-arg callable that runs one PhenoDigm stage."""

    def _run():
        getattr(PhenoDigm(DB), name)(**overrides)

    return _run


with DAG(
    dag_id="phenodigm2_pipeline",
    description="Build the PhenoDigm2 database and Solr core",
    start_date=datetime(2024, 1, 1),
    schedule=None,
    catchup=False,
    tags=["phenodigm2"],
) as dag:
    tasks = {
        name: PythonOperator(task_id=name, python_callable=_stage(name))
        for name in STAGES
    }

    # Wire a linear pipeline: download >> build >> ... >> solr
    previous = None
    for name in STAGES:
        if previous is not None:
            previous >> tasks[name]
        previous = tasks[name]
