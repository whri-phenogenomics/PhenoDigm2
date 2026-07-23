"""Use contents of PhenoDigm2 db to create a Solr core.

@author: Tomasz Konopka
"""

import os.path
from importlib.resources import as_file
from pathlib import Path
from shutil import copy2, copytree, rmtree

import requests

from . import tools as pd2tools
from . import solrdata
from . import solrlinks
from . import solrsearch


SOLR_OUTPUT_DIRECTORY = Path("output") / "solr"
SOLR_COMPOSE_FILENAME = "dc-solr-7.5.yml"
SOLR_CORES_DIRECTORY = "solrcores7.5"


# ############################################################################
#


def prepareSolrBundle(config):
    """Prepare the Compose file and core storage inside a build bundle.

    Existing Compose files are preserved so release-specific edits are not
    overwritten. A custom ``--solr_cores_dir`` remains supported; otherwise
    the path matches the relative volume in the bundled Compose file.
    """

    bundle_dir = Path(config.db).expanduser().resolve()
    output_dir = bundle_dir / SOLR_OUTPUT_DIRECTORY
    compose_path = output_dir / SOLR_COMPOSE_FILENAME
    default_cores_dir = output_dir / SOLR_CORES_DIRECTORY

    output_dir.mkdir(parents=True, exist_ok=True)
    if not compose_path.exists():
        resources_dir = Path(pd2tools.getPD2dirs(config)[2])
        release_source = resources_dir / SOLR_COMPOSE_FILENAME
        if release_source.is_file():
            copy2(release_source, compose_path)
        else:
            bundled_source = pd2tools.getBundledResourcesDir() / SOLR_COMPOSE_FILENAME
            with as_file(bundled_source) as source_path:
                copy2(source_path, compose_path)

    configured_cores_dir = getattr(config, "solr_cores_dir", None)
    if configured_cores_dir:
        cores_dir = Path(configured_cores_dir).expanduser().resolve()
    else:
        cores_dir = default_cores_dir
        config.solr_cores_dir = str(cores_dir)
    cores_dir.mkdir(parents=True, exist_ok=True)

    return output_dir, compose_path, cores_dir


def initSolrCore(config):
    """Initialize a solr core (solr version 7)"""

    # identify file system directories
    t1, downdir, resdir, dbdir = pd2tools.getPD2dirs(config)
    corename = config.solr_corename

    if config.db.endswith("/"):
        config.db = config.db[:-1]
    instance = "phenodigm2_" + os.path.basename(config.db)

    # get paths to directories for this core
    coredir, confdir, datadir = pd2tools.getPD2coredir(config)

    # unload and remove the core if necessary
    if os.path.exists(coredir):
        pd2tools.log("Removing existing core", 2)
        unload = config.solr_url + "admin/cores?action=UNLOAD"
        unload += "&core=" + corename + "&deleteIndex=true"
        requests.get(url=unload)
        rmtree(coredir)

    pd2tools.log("Creating new solr core", 2)

    # create the core from scratch and let anyone read/write execute
    for onedir in [coredir, datadir]:
        if not os.path.exists(onedir):
            os.makedirs(onedir)
            os.chmod(onedir, 0o777)

    # copy the configuration files into the conf directory
    resdirsolr = os.path.join(resdir, "solr7")
    copytree(resdirsolr, confdir)
    os.chmod(confdir, 0o777)

    create = config.solr_url + "admin/cores?action=CREATE"
    create += "&name=" + corename + "&instanceDir=mycores/" + instance

    r = requests.get(url=create)
    return r.ok, r.text


# ############################################################################
# run this from outside of module


def runSolrCoreBuild(config):
    """Create and fill a Solr core from a Phenodigm2 db."""

    output_dir, compose_path, cores_dir = prepareSolrBundle(config)
    pd2tools.log(f"Solr bundle directory: {output_dir}", 2)
    pd2tools.log(f"Docker Compose file: {compose_path}", 2)
    pd2tools.log(f"Solr cores directory: {cores_dir}", 2)

    pd2tools.log("Initializing solr core")
    init_ok, init_result = initSolrCore(config)
    if not init_ok:
        raise RuntimeError("Failed to create Solr core: " + init_result)

    pd2tools.log("Transferring data to solr core")
    solrlinks.runSolrGeneGene(config)
    solrdata.runSolrOntologies(config)
    solrdata.runSolrGenes(config)
    solrdata.runSolrDiseases(config)
    solrdata.runSolrMouseModels(config)
    solrlinks.runSolrOntoOnto(config)
    solrsearch.runSolrDiseaseSearch(config)
    solrlinks.runSolrDiseaseGenes(config)
    solrlinks.runSolrDiseaseModels(config)
