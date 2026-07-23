"""Programmatic Python API for PhenoDigm2.

Wraps the same pipeline stages the ``phenodigm2`` CLI runs (see :mod:`cli`) so
they can be driven as ordinary Python methods -- e.g. from Airflow
``PythonOperator`` callables::

    from phenodigm2 import PhenoDigm

    pd = PhenoDigm("/data/pd2_build")
    pd.download()
    pd.build()
    pd.score()

Each method mirrors one CLI action. Configuration is produced by reusing
:func:`cli.build_parser`, so the API and CLI share a single source of default
values and never drift. Unlike the CLI, methods never call ``exit()``: on
failure the underlying exception propagates so an orchestrator (Airflow) can
mark the task failed.
"""

from .cli import build_parser
from . import dbbuild
from . import export
from . import load
from . import ontology_mapping
from . import parquet as parquet_export
from . import post_process
from . import prep
from . import query
from . import score
from . import solr
from . import status
from . import tools


class PhenoDigm:
    """Facade over the PhenoDigm2 pipeline stages, bound to one build dir.

    Arguments:
    db -- build directory (the same value as the CLI ``--db`` flag).
    options -- global overrides applied to every stage's config, e.g.
        ``PhenoDigm(db, cores=8, verbose=True)``. Per-call ``**overrides``
        passed to a method take precedence over these.
    """

    def __init__(self, db, **options):
        self.db = str(db)
        self.options = options

    # -- config construction -------------------------------------------------

    def _config(self, action, **overrides):
        """Build the argparse ``Namespace`` the stage functions expect.

        ``action`` is the CLI action token (hyphenated where relevant, e.g.
        ``"solr-prepare"``); ``"init"`` is special-cased onto the ``--init``
        flag path. Reusing the CLI parser guarantees the same defaults as the
        command line.
        """

        argv = ["--init", self.db] if action == "init" else [action]
        config = build_parser().parse_args(argv)
        config.db = self.db
        for key, value in {**self.options, **overrides}.items():
            setattr(config, key, value)
        if action != "init":
            config.dbfile = tools.getDBfile(config)
        return config

    # -- pipeline stages (mirror cli.main dispatch) --------------------------

    def init(self, **overrides):
        """Create a new build directory and seed its bundled resources."""
        config = self._config("init", **overrides)
        tools.log("Starting PhenoDigm2 [init]")
        prep.runDirInit(config)
        tools.log("Done")
        return self

    def download(self, **overrides):
        """Download all declared data inputs into the build directory."""
        config = self._config("download", **overrides)
        tools.log("Starting PhenoDigm2 [download]")
        prep.runDirPrep(config)
        prep.runDownloads(config)
        tools.log("Done")
        return self

    def build(self, **overrides):
        """Build the sqlite db schema and load all data tables."""
        config = self._config("build", **overrides)
        tools.log("Starting PhenoDigm2 [build]")
        prep.runDirPrep(config)
        dbbuild.runDBBuild(config)
        load.runLoad(config)
        tools.log("Done")
        return self

    def ontology_mapping(self, **overrides):
        """Process Phenio mappings into ontology term-term similarities."""
        config = self._config("ontology_mapping", **overrides)
        tools.log("Starting PhenoDigm2 [ontology_mapping]")
        prep.runDirPrep(config)
        ontology_mapping.run_ontology_mapping_processing(config)
        tools.log("Done")
        return self

    def score(self, **overrides):
        """Compute PhenoDigm2 disease-model association scores."""
        config = self._config("score", **overrides)
        tools.log("Starting PhenoDigm2 [score]")
        prep.runDirPrep(config)
        score.runScoring(config)
        tools.log("Done")
        return self

    def index(self, **overrides):
        """Create sqlite indexes on the built db."""
        config = self._config("index", **overrides)
        tools.log("Starting PhenoDigm2 [index]")
        prep.runDirPrep(config)
        dbbuild.runDBIndexing(config)
        tools.log("Done")
        return self

    def solr_prepare(self, **overrides):
        """Prepare the Solr bundle (CLI action ``solr-prepare``)."""
        config = self._config("solr-prepare", **overrides)
        tools.log("Starting PhenoDigm2 [solr-prepare]")
        prep.runDirPrep(config)
        solr.prepareSolrBundle(config)
        tools.log("Done")
        return self

    def solr(self, **overrides):
        """Build the Solr core."""
        config = self._config("solr", **overrides)
        tools.log("Starting PhenoDigm2 [solr]")
        prep.runDirPrep(config)
        solr.runSolrCoreBuild(config)
        tools.log("Done")
        return self

    def post_process(self, **overrides):
        """Run the post-processing workflow (CLI action ``post-process``)."""
        config = self._config("post-process", **overrides)
        tools.log("Starting PhenoDigm2 [post-process]")
        prep.runDirPrep(config)
        post_process.run_post_process(config)
        tools.log("Done")
        return self

    def parquet(self, **overrides):
        """Export the Parquet bundle from the built db."""
        config = self._config("parquet", **overrides)
        tools.log("Starting PhenoDigm2 [parquet]")
        parquet_export.runParquetBundleBuild(config)
        tools.log("Done")
        return self

    # -- inspection / one-off actions ----------------------------------------
    # These mirror the CLI, which prints results to stdout. They return
    # whatever the underlying stage returns (mostly None today); returning
    # structured data is future work.

    def status(self, **overrides):
        """Report build status."""
        return status.statusPhenodigm(self._config("status", **overrides))

    def query(self, **overrides):
        """Run a manual query against the built db (see --disease/--model/...)."""
        return query.queryPhenodigm(self._config("query", **overrides))

    def compute(self, **overrides):
        """Compute scores for a single disease/model on demand."""
        return query.computeScores(self._config("compute", **overrides))

    def export(self, **overrides):
        """Export db tables (see --table/--where)."""
        return export.exportTables(self._config("export", **overrides))

    # -- generic dispatch ----------------------------------------------------

    def run(self, action, **overrides):
        """Run any stage by its CLI action name (hyphen or underscore)."""
        return getattr(self, action.replace("-", "_"))(**overrides)
