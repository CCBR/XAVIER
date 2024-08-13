import argparse
import glob
import os
import tempfile
from ccbr_tools.pipeline.util import get_tmp_dir, get_hpcname
from ccbr_tools.pipeline.cache import get_sif_cache_dir
from ccbr_tools.shell import exec_in_context

from xavier.src.xavier.util import xavier_base
from xavier.src.xavier.run import run


def test_dryrun():
    if get_hpcname() == "biowulf":
        with tempfile.TemporaryDirectory() as tmp_dir:
            run_args = argparse.Namespace(
                runmode="init",
                input=list(glob.glob(xavier_base(".tests/*.fastq.gz"))),
                output=tmp_dir,
                genome="hg38",
                targets=xavier_base("resources/Agilent_SSv7_allExons_hg38.bed"),
                mode="local",
                job_name="pl:xavier",
                callers=["mutect2", "mutect", "strelka", "vardict", "varscan"],
                pairs=xavier_base(".tests/pairs.tsv"),
                ffpe=False,
                cnv=False,
                wait=False,
                create_nidap_folder=False,
                silent=False,
                singularity_cache=os.environ.get("SINGULARITY_CACHEDIR", None),
                sif_cache=get_sif_cache_dir(),
                tmp_dir=get_tmp_dir(None, tmp_dir),
                threads=2,
            )
            # init
            allout_1 = exec_in_context(run, run_args)
            run_args.runmode = "dryrun"
            # dryrun
            allout_2 = exec_in_context(run, run_args)
        assert all(
            [
                "--Initializing" in allout_1,
                "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
                in allout_2,
            ]
        )
