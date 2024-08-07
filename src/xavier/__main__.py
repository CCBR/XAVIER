#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""XAVIER: eXome Analysis and Variant explorER:
A highly reproducible and portable Whole Exome-seq data analysis pipeline

ABOUT: This is the main entry for the XAVIER pipeline.
REQUIRES:
  - python>=3.5
  - snakemake   (recommended>=6.0.0 **DO NOT USE snakemake v7.7.0 (bug reported!)**)
  - singularity (recommended==latest)
DISCLAIMER:
                    PUBLIC DOMAIN NOTICE
        CCR Collaborative Bioinformatics Resource (CCBR)
                National Cancer Institute (NCI)

This software/database is a "United  States Government Work" under
the terms of the United  States Copyright Act.  It was written as
part of the author's official duties as a United States Government
employee and thus cannot be copyrighted. This software is freely
available to the public for use.

Although all  reasonable  efforts have been taken  to ensure  the
accuracy and reliability of the software and data, CCBR do not and
cannot warrant the performance or results that may  be obtained by
using this software or data. CCBR and NCI disclaim all warranties,
express  or  implied,  including   warranties   of   performance,
merchantability or fitness for any particular purpose.

Please cite the author and the "NIH Biowulf Cluster" in any work or
product based on this material.
USAGE:
  $ xavier <run> [OPTIONS]
EXAMPLE:
  $ xavier run --input *.R?.fastq.gz --output xavier_hg38/ --genome hg38
"""

# Python standard library
from __future__ import print_function
import sys, os, subprocess, re, json, textwrap


# 3rd party imports from pypi
import argparse  # potential python3 3rd party package, added in python/3.5

# Local imports
from .run import init, setup, bind, dryrun, runner
from .shells import bash
from .options import genome_options
from .util import err, exists, fatal, permissions, check_cache, require, get_version

__version__ = get_version()
__email__ = "ccbr@mail.nih.gov"
__home__ = os.path.dirname(os.path.abspath(__file__))


def run(sub_args):
    """Initialize, setup, and run the XAVIER pipeline.
    Calls initialize() to create output directory and copy over pipeline resources,
    setup() to create the pipeline config file, dryrun() to ensure their are no issues
    before running the pipeline, and finally run() to execute the Snakemake workflow.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    """
    # Step 0. Check for required dependencies
    # The pipelines has only two requirements:
    # snakemake and singularity
    require(["snakemake", "singularity"], ["snakemake", "singularity"])

    # Optional Step. Initialize working directory,
    # copy over required resources to run
    # the pipeline
    git_repo = __home__
    if sub_args.runmode == "init":
        print("--Initializing")
        input_files = init(
            repo_path=git_repo, output_path=sub_args.output, links=sub_args.input
        )

    # Required Step. Setup pipeline for execution,
    # dynamically create config.json config
    # file from user inputs and base config
    # determine "nidap folder"
    create_nidap_folder_YN = "no"
    if sub_args.create_nidap_folder:
        create_nidap_folder_YN = "yes"

    # templates
    config = setup(
        sub_args,
        repo_path=git_repo,
        output_path=sub_args.output,
        create_nidap_folder_YN=create_nidap_folder_YN,
        links=sub_args.input,
    )

    # Required Step. Resolve docker/singularity bind
    # paths from the config file.
    bindpaths = bind(sub_args, config=config)

    # Optional Step: Dry-run pipeline
    # if sub_args.dry_run:
    if sub_args.runmode == "dryrun" or sub_args.runmode == "run":
        print("--Dry-Run")
        # Dryrun pipeline
        dryrun_output = dryrun(
            outdir=sub_args.output
        )  # python3 returns byte-string representation
        print(
            "\nDry-running XAVIER pipeline:\n{}".format(dryrun_output.decode("utf-8"))
        )

    # Optional Step. Orchestrate pipeline execution,
    # run pipeline in locally on a compute node
    # for debugging purposes or submit the master
    # job to the job scheduler, SLURM, and create
    # logging file
    if sub_args.runmode == "run":
        print("--Run full pipeline")
        if not exists(os.path.join(sub_args.output, "logfiles")):
            # Create directory for logfiles
            os.makedirs(os.path.join(sub_args.output, "logfiles"))
        if sub_args.mode == "local":
            log = os.path.join(sub_args.output, "logfiles", "snakemake.log")
        else:
            log = os.path.join(sub_args.output, "logfiles", "master.log")
        logfh = open(log, "w")
        wait = ""
        if sub_args.wait:
            wait = "--wait"
        mjob = runner(
            mode=sub_args.mode,
            outdir=sub_args.output,
            # additional_bind_paths = all_bind_paths,
            alt_cache=sub_args.singularity_cache,
            threads=int(sub_args.threads),
            jobname=sub_args.job_name,
            submission_script="runner",
            logger=logfh,
            additional_bind_paths=",".join(bindpaths),
            tmp_dir=sub_args.tmp_dir,
            wait=wait,
        )

        # Step 5. Wait for subprocess to complete,
        # this is blocking and not asynchronous
        if not sub_args.silent:
            print("\nRunning XAVIER pipeline in '{}' mode...".format(sub_args.mode))
        mjob.wait()
        logfh.close()

        # Step 6. Relay information about submission
        # of the master job or the exit code of the
        # pipeline that ran in local mode
        if sub_args.mode == "local":
            if int(mjob.returncode) == 0:
                print("XAVIER has successfully completed")
            else:
                fatal(
                    "XAVIER failed. Please see {} for more information.".format(
                        os.path.join(sub_args.output, "logfiles", "snakemake.log")
                    )
                )
        elif sub_args.mode == "slurm":
            jobid = (
                open(os.path.join(sub_args.output, "logfiles", "mjobid.log"))
                .read()
                .strip()
            )
            if not sub_args.silent:
                if int(mjob.returncode) == 0:
                    print("Successfully submitted master job: ", end="")
                else:
                    fatal("Error occurred when submitting the master job.")
            print(jobid)


def unlock(sub_args):
    """Unlocks a previous runs output directory. If snakemake fails ungracefully,
    it maybe required to unlock the working directory before proceeding again.
    This is rare but it does occasionally happen. Maybe worth add a --force
    option to delete the '.snakemake/' directory in the future.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for unlock sub-command
    """

    print("Unlocking the pipeline's output directory...")
    outdir = sub_args.output

    try:
        unlock_output = subprocess.check_output(
            ["snakemake", "--unlock", "--cores", "1", "--configfile=config.json"],
            cwd=outdir,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        # Unlocking process returned a non-zero exit code
        sys.exit("{}\n{}".format(e, e.output))

    print("Successfully unlocked the pipeline's working directory!")


def cache(sub_args):
    """Caches remote resources or reference files stored on DockerHub and S3.
    Local SIFs will be created from images defined in 'config/containers/images.json'.
    @TODO: add option to cache other shared S3 resources (i.e. kraken db and fqscreen indices)
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for unlock sub-command
    """
    print(sub_args)

    sif_cache = sub_args.sif_cache
    # Get absolute PATH to templates in XAVIER git repo
    repo_path = os.path.dirname(os.path.abspath(__file__))
    images = os.path.join(repo_path, "config", "containers", "images.json")

    # Create image cache
    if not exists(sif_cache):
        # Pipeline output directory does not exist on filesystem
        os.makedirs(sif_cache)
    elif exists(sif_cache) and os.path.isfile(sif_cache):
        # Provided Path for pipeline output directory exists as file
        raise OSError(
            """\n\tFatal: Failed to create provided sif cache directory!
        User provided --sif-cache PATH already exists on the filesystem as a file.
        Please {} cache again with a different --sif-cache PATH.
        """.format(
                sys.argv[0]
            )
        )

    # Check if local SIFs already exist on the filesystem
    with open(images, "r") as fh:
        data = json.load(fh)

    pull = []
    for image, uri in data["images"].items():
        sif = os.path.join(
            sif_cache, "{}.sif".format(os.path.basename(uri).replace(":", "_"))
        )
        if not exists(sif):
            # If local sif does not exist on in cache, print warning
            # and default to pulling from URI in config/containers/images.json
            print('Image will be pulled from "{}".'.format(uri), file=sys.stderr)
            pull.append(uri)

    if not pull:
        # Nothing to do!
        print("Singularity image cache is already up to update!")
    else:
        # There are image(s) that need to be pulled
        if not sub_args.dry_run:
            # submission_script for XAVIER cache is /path/to/output/resources/cacher
            # Quote user provided values to avoid shell injections
            masterjob = subprocess.Popen(
                "sbatch --parsable -J pl:cache --gres=lscratch:200  --time=10:00:00 --mail-type=BEGIN,END,FAIL "
                + str(os.path.join(repo_path, "resources", "cacher"))
                + " slurm "
                + " -s '{}' ".format(sif_cache)
                + " -i '{}' ".format(",".join(pull))
                + " -t '/lscratch/${SLURM_JOB_ID}/.singularity/' ",
                cwd=sif_cache,
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
            )

            masterjob.communicate()
            print(
                "XAVIER reference cacher submitted master job with exit-code: {}".format(
                    masterjob.returncode
                )
            )


def parsed_arguments():
    """Parses user-provided command-line arguments. Requires argparse and textwrap
    package. argparse was added to standard lib in python 3.5 and textwrap was added
    in python 3.5. To create custom help formatting for subparsers a docstring is
    used create the help message for required options. argparse does not support named
    subparser groups, which is normally what would be used to accomplish this reformatting.
    As so, the help message for require options must be suppressed. If a new required arg
    is added to a subparser, it must be added to the docstring and the usage statement
    also must be updated.
    """

    # Create a top-level parser
    parser = argparse.ArgumentParser(
        prog="xavier", description="XAVIER: eXome Analysis and Variant explorER:"
    )

    # Adding Version information
    parser.add_argument(
        "--version", action="version", version="%(prog)s {}".format(__version__)
    )

    # Create sub-command parser
    subparsers = parser.add_subparsers(help="List of available sub-commands")

    # Sub-parser for the "run" sub-command
    # Grouped sub-parser arguments are currently not supported by argparse.
    # https://bugs.python.org/issue9341
    # Here is a work around to create more useful help message for named
    # options that are required! Please note: if a required arg is added the
    # description below should be updated (i.e. update usage and add new option)
    required_run_options = textwrap.dedent(
        """\
        usage: xavier run [--help] \\
                              [--mode {local, slurm}] \\
                              [--job-name JOB_NAME] \\
                              [--callers {mutect2,mutect,strelka, ...}] \\
                              [--pairs PAIRS] \\
                              [--ffpe] \\
                              [--cnv] \\
                              [--silent] \\
                              [--singularity-cache SINGULARITY_CACHE] \\
                              [--sif-cache SIF_CACHE] \\
                              [--tmpdir TMP_DIR] \\
                              [--threads THREADS] \\
                              [--wait] \\
                              [--create-nidap-folder] \\
                              --runmode RUNMODE [init, dryrun, run]
                              --input INPUT [INPUT ...] \\
                              --output OUTPUT \\
                              --genome {hg38, mm10, ...} \\
                              --targets TARGETS

        required arguments:
          --runmode RUNMODE [init, dryrun, run ...]
                                Runmode. Determine the runmode to deploy for the pipeline
                                1) init. Initialize and prepare the output directory
                                2) dry-run. Dry-run the pipeline
                                3) run. Run the pipeline
                                Example: --runmode init
          --input INPUT [INPUT ...]
                                Input FastQ or BAM file(s) to process. One or more input
                                files can be provided.  The pipeline does NOT support
                                single-end WES data. Please provide either a set of
                                FastQ files or a set of BAM files. The pipeline does
                                NOT support processing a mixture of FastQ files and
                                BAM files.
                                Example: --input .tests/*.R?.fastq.gz
          --output OUTPUT
                                Path to an output directory. This location is where
                                the pipeline will create all of its output files, also
                                known as the pipeline's working directory. If the user
                                provided working directory has not been initialized,
                                it will be created automatically.
                                Example: --output /data/$USER/xavier_hg38
          --genome {hg38, mm10, ...}
                                Reference genome. This option defines the reference
                                genome of the samples. Currently, hg38 and mm10 are supported.
                                Support for additional and custom genomes will be added soon.
                                Example: --genome hg38
          --targets TARGETS
                                Path to exome targets BED file. This file can be
                                obtained from the manufacturer of the target capture
                                kit that was used.

        """
    )

    # Display example usage in epilog
    run_epilog = textwrap.dedent(
        """\
        example:
          # Step 1.) Grab an interactive node (do not run on head node)
          sinteractive --mem=8g --cpus-per-task=4
          module purge
          module load ccbrpipeliner

          # Step 2A.) Initialize the pipeline
          xavier run \\
                        --runmode init \\
                        --input .tests/*.R?.fastq.gz \\
                        --output /data/$USER/xavier_hg38 \\
                        --genome hg38 \\
                        --targets .tests/Agilent_SSv7_allExons_hg38.bed

          # Step 2B.) Dry-run the pipeline
          xavier run \\
                        --runmode dryrun \\
                        --input .tests/*.R?.fastq.gz \\
                        --output /data/$USER/xavier_hg38 \\
                        --genome hg38 \\
                        --targets Agilent_SSv7_allExons_hg38.bed \\
                        --mode slurm \\

          # Step 2C.) Run the XAVIER pipeline
          # The slurm mode will submit jobs to the cluster.
          # It is recommended running xavier in this mode.
          xavier run \\
                        --runmode run \\
                        --input .tests/*.R?.fastq.gz \\
                        --output /data/$USER/xavier_hg38 \\
                        --genome hg38 \\
                        --targets .tests/Agilent_SSv7_allExons_hg38.bed \\
                        --mode slurm

        version:
          {}
        """.format(
            __version__
        )
    )

    # Suppressing help message of required args to overcome no sub-parser named groups
    subparser_run = subparsers.add_parser(
        "run",
        help="Run the XAVIER  pipeline with input files.",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=required_run_options,
        epilog=run_epilog,
    )

    # Required Arguments
    # Input FastQ files
    subparser_run.add_argument(
        "--runmode",
        # Determines how to run the pipeline: init, run, or dry-run
        required=True,
        choices=["init", "run", "dryrun"],
        type=str,
        help=argparse.SUPPRESS,
    )

    # Input FastQ files
    subparser_run.add_argument(
        "--input",
        # Check if the file exists and if it is readable
        type=lambda file: permissions(parser, file, os.R_OK),
        required=True,
        nargs="+",
        help=argparse.SUPPRESS,
    )

    # Output Directory (analysis working directory)
    subparser_run.add_argument(
        "--output",
        type=lambda option: os.path.abspath(os.path.expanduser(option)),
        required=True,
        help=argparse.SUPPRESS,
    )

    # Reference Genome (to dynamically select reference files)
    subparser_run.add_argument(
        "--genome",
        required=True,
        # choices = ['hg38', 'mm10'],
        type=lambda option: str(
            genome_options(subparser_run, option, ["hg38", "mm10"])
        ),
        help=argparse.SUPPRESS,
    )

    # Exome TARGET BED file
    subparser_run.add_argument(
        "--targets",
        # Check if the file exists and if it is readable
        type=lambda file: permissions(parser, file, os.R_OK),
        required=True,
        help=argparse.SUPPRESS,
    )

    # Optional Arguments
    # Execution Method (run locally on a compute node, submit to SLURM job scheduler, etc.)
    subparser_run.add_argument(
        "--mode",
        type=str,
        required=False,
        default="slurm",
        choices=["slurm", "local"],
        help="Execution Method [Default: slurm]. Defines the mode or method of execution. \
        Valid mode options include: local or slurm. \
        local: uses local method of execution. local executions will run serially on \
        compute instance. This is useful for testing, debugging, or when a users does \
        not have access to a high performance computing environment. If this option is \
        not provided, it will default to a local execution mode. \
        slurm: uses slurm and singularity backend. The slurm execution method will submit \
        jobs to a cluster. It is recommended running xavier in this mode as execution \
        will be significantly faster in a distributed environment. \
        Example: --mode slurm",
    )

    # Name of master job
    subparser_run.add_argument(
        "--job-name",
        type=str,
        required=False,
        default="pl:xavier",
        help="Set the name of the pipeline's master job. \
        When submitting the pipeline to a job scheduler, like SLURM, \
        this option always you to set the name of the pipeline's master \
        job. By default, the name of the pipeline's master job \
        is set to \"pl:xavier\". \
        Example: --job-name xavier_hg38_main",
    )

    # Variant Callers
    subparser_run.add_argument(
        "--callers",
        type=str,
        required=False,
        nargs="+",
        metavar="CALLERS",
        default=["mutect2", "mutect", "strelka", "vardict", "varscan"],
        choices=["mutect2", "mutect", "strelka", "vardict", "varscan"],
        help="Variant Callers. List of variant callers to call mutations. Please select from one or \
        more of the following options: [mutect2, mutect, strelka, vardict, varscan]. Defaults to using all \
        variant callers. Example: --callers mutect2 strelka varscan",
    )

    # Tumor normal pairs file
    subparser_run.add_argument(
        "--pairs",
        # Check if the file exists and if it is readable
        type=lambda file: permissions(parser, file, os.R_OK),
        required=False,
        help='Tumor normal pairs file. This tab delimited file contains two columns with the names \
        of tumor and normal pairs, one per line. The header of the file needs to be "Tumor" for the \
        tumor column and "Normal" for the normal column.',
    )

    # Correction for FFPE samples
    subparser_run.add_argument(
        "--ffpe",
        action="store_true",
        required=False,
        default=False,
        help="FFPE correction. Runs an additional filtering step for Formalin-Fixed Paraffin-Embedded \
        (FFPE) samples. Do NOT use this option with non-FFPE samples.",
    )

    # Call CNVs
    subparser_run.add_argument(
        "--cnv",
        action="store_true",
        required=False,
        default=False,
        help="Call copy number variations or CNVs. CNVs will only be called from tumor-normal pairs. \
        If this option is provided without providing a --pairs file, CNVs will NOT be called.",
    )

    # wait until master job finishes ... required for HPC API execution
    subparser_run.add_argument(
        "--wait",
        action="store_true",
        required=False,
        default=False,
        help="Wait until master job completes. This is required if \
        the job is submitted using HPC API. If not provided \
        the API may interpret submission of master job as \
        completion of the pipeline!",
    )

    # create-nidap-folder create a folder called "NIDAP" to be moved back to NIDAP
    subparser_run.add_argument(
        "--create-nidap-folder",
        action="store_true",
        required=False,
        default=False,
        help='Create folder called "NIDAP" with file to-be-moved back to NIDAP \
              This makes it convenient to move only this folder (called NIDAP) and its content back \
              to NIDAP, rather than the entire pipeline output folder',
    )

    # Silent output mode
    subparser_run.add_argument(
        "--silent",
        action="store_true",
        required=False,
        default=False,
        help="Silence standard output. Reduces the amount of information directed \
        to standard output when submitting master job to the job scheduler. Only the \
        job id of the master job is returned.",
    )

    # Singularity cache directory (default uses output directory)
    subparser_run.add_argument(
        "--singularity-cache",
        type=lambda option: check_cache(
            parser, os.path.abspath(os.path.expanduser(option))
        ),
        required=False,
        help="Overrides the $SINGULARITY_CACHEDIR environment variable. Singularity will cache \
        image layers pulled from remote registries. By default, the cache is set to \
        '/path/to/output/directory/.singularity/'. \
        Please note that this cache cannot be shared across users.",
    )

    # Local SIF cache directory (default pull from Dockerhub)
    subparser_run.add_argument(
        "--sif-cache",
        type=lambda option: os.path.abspath(os.path.expanduser(option)),
        required=False,
        help="Path where a local cache of SIFs are stored. \
        This cache can be shared across users if permissions are \
        set correctly. If a SIF does not exist in the SIF cache, \
        the image will be pulled from Dockerhub. The xavier cache \
        subcommand can be used to create a local SIF cache. Please see \
        xavier cache for more information.",
    )

    # Base directory to write temporary files
    subparser_run.add_argument(
        "--tmp-dir",
        type=str,
        required=False,
        default="/lscratch/$SLURM_JOBID/",
        help="Path on the filesystem for writing intermediate, temporary output \
        files. By default, this variable is set to '/lscratch/$SLURM_JOBID' \
        for backwards compatibility with the NIH's Biowulf cluster; however, \
        if you are running the pipeline on another cluster, this option will \
        need to be specified. Ideally, this path should point to a dedicated \
        location on the filesystem for writing tmp files. On many systems, this \
        location is set to somewhere in /scratch. If you need to inject a variable \
        into this string that should NOT be expanded, please quote this options \
        value in single quotes. As an example, on the NCI/NIH FRCE cluster the \
        value of this option would be set to \
        --tmp-dir '/scratch/cluster_scratch/$USER/', \
        default:  '/lscratch/$SLURM_JOBID/'",
    )

    # Number of threads for the xavier pipeline's main proceess
    subparser_run.add_argument(
        "--threads",
        type=int,
        required=False,
        default=2,
        help="Max number of threads for local processes. It is recommended \
        setting this value to the maximum number of CPUs available on the host \
        machine, default: 2.",
    )

    # Sub-parser for the "unlock" sub-command
    # Grouped sub-parser arguments are currently not supported.
    # https://bugs.python.org/issue9341
    # Here is a work around to create more useful help message for named
    # options that are required! Please note: if a required arg is added the
    # description below should be updated (i.e. update usage and add new option)
    required_unlock_options = textwrap.dedent(
        """\
        usage: xavier unlock [-h] --output OUTPUT

        If the pipeline fails ungracefully, it maybe required to unlock the working
        directory before proceeding again. Please verify that the pipeline is not
        running before running this command. If the pipeline is still running, the
        workflow manager will report the working directory is locked. This is normal
        behavior. Do NOT run this command if the pipeline is still running.

        required arguments:
          --output OUTPUT
                                Path to a previous run's output directory to
                                unlock. This will remove a lock on the working
                                directory. Please verify that the pipeline is
                                not running before running this command.
                                Example: --output /data/$USER/xavier_hg38

        """
    )

    # Display example usage in epilog
    unlock_epilog = textwrap.dedent(
        """\
        example:
          # Unlock xavier output directory
          xavier unlock --output /scratch/$USER/xavier_hg38

        version:
          {}
        """.format(
            __version__
        )
    )

    # Suppressing help message of required args to overcome no sub-parser named groups
    subparser_unlock = subparsers.add_parser(
        "unlock",
        help="Unlocks a previous runs output directory.",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=required_unlock_options,
        epilog=unlock_epilog,
    )

    # Required Arguments
    # Output Directory (analysis working directory)
    subparser_unlock.add_argument(
        "--output", type=str, required=True, help=argparse.SUPPRESS
    )

    # Sub-parser for the "cache" sub-command
    # Grouped sub-parser arguments are currently not supported.
    # https://bugs.python.org/issue9341
    # Here is a work around to create more useful help message for named
    # options that are required! Please note: if a required arg is added the
    # description below should be updated (i.e. update usage and add new option)
    required_cache_options = textwrap.dedent(
        """\
        usage: xavier cache [-h] [-n] --sif-cache SIF_CACHE

        Creates a local cache resources hosted on DockerHub or AWS S3.
        These resources are normally pulled onto the filesystem when the
        pipeline runs; however, due to network issues or DockerHub pull
        rate limits, it may make sense to pull the resources once so a
        shared cache can be created. It is worth noting that a singularity
        cache cannot normally be shared across users. Singularity strictly
        enforces that a cache is owned by the user. To get around this
        issue, the cache subcommand can be used to create local SIFs on
        the filesystem from images on DockerHub.

        required arguments:
          --sif-cache SIF_CACHE
                      Path where a local cache of SIFs will be stored.
                      Images defined in config/containers/images.json
                      will be pulled into the local filesystem. The
                      path provided to this option can be passed to
                      the --sif-cache option of the run sub command.
                      Please see xavier run sub command for more
                      information.

                      Example: --sif-cache /scratch/$USER/cache

        """
    )

    # Display example usage in epilog
    cache_epilog = textwrap.dedent(
        """\
        example:
          # Cache xavier resources
          xavier cache --sif-cache /scratch/$USER/cache

        version:
          {}
        """.format(
            __version__
        )
    )

    # Suppressing help message of required args to overcome no sub-parser named groups
    subparser_cache = subparsers.add_parser(
        "cache",
        help="Cache remote resources locally.",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=required_cache_options,
        epilog=cache_epilog,
    )

    # Required Arguments
    # Output Directory (analysis working directory)
    subparser_cache.add_argument(
        "--sif-cache",
        type=lambda option: os.path.abspath(os.path.expanduser(option)),
        required=True,
        help=argparse.SUPPRESS,
    )

    # Optional Arguments
    # Dry-run xavier cache (do not pull any remote resources)
    subparser_cache.add_argument(
        "--dry-run",
        action="store_true",
        required=False,
        default=False,
        help="Only display what remote resources would be pulled.",
    )

    # Define handlers for each sub-parser
    subparser_run.set_defaults(func=run)
    subparser_unlock.set_defaults(func=unlock)
    subparser_cache.set_defaults(func=cache)

    # Parse command-line args
    args = parser.parse_args()
    return args


def main():
    # show helpful error message when no arguments given
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal("Invalid usage: xavier [-h] [--version] ...")

    # Collect args for sub-command
    args = parsed_arguments()

    # Display version information
    err("xavier ({})".format(__version__))

    # Mediator method to call sub-command's set handler function
    args.func(args)


if __name__ == "__main__":
    main()
