#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
import contextlib
import io
import os
import re
import json
import shutil
import sys
import subprocess
import datetime
from ccbr_tools.pipeline.util import (
    git_commit_hash,
    join_jsons,
    fatal,
    which,
    exists,
    err,
    require,
    get_hpcname,
)
from ccbr_tools.pipeline.cache import image_cache

# Local imports
from .util import get_version, xavier_base


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
    require(["snakemake", "singularity"], ["snakemake/7", "singularity"])

    # Optional Step. Initialize working directory,
    # copy over required resources to run
    # the pipeline
    git_repo = xavier_base()
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
    wait = ""
    if sub_args.wait:
        wait = "--wait"

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


def init(
    repo_path, output_path, links=[], required=["workflow", "resources", "config"]
):
    """Initialize the output directory. If user provides a output
    directory path that already exists on the filesystem as a file
    (small chance of happening but possible), a OSError is raised. If the
    output directory PATH already EXISTS, it will not try to create the directory.
    @param repo_path <str>:
        Path to RNA-seek source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @param links list[<str>]:
        List of files to symlink into output_path
    @param required list[<str>]:
        List of folder to copy over into output_path
    """
    if not exists(output_path):
        # Pipeline output directory
        # does not exist on filesystem
        os.makedirs(output_path)

    elif exists(output_path) and os.path.isfile(output_path):
        # Provided Path for pipeline
        # output directory exists as file
        raise OSError(
            """\n\tFatal: Failed to create provided pipeline output directory!
        User provided --output PATH already exists on the filesystem as a file.
        Please run {} again with a different --output PATH.
        """.format(
                sys.argv[0]
            )
        )

    # Copy over templates are other required resources
    copy_safe(source=repo_path, target=output_path, resources=required)

    # Create renamed symlinks for each rawdata
    # file provided as input to the pipeline
    inputs = sym_safe(input_data=links, target=output_path)

    return inputs

def _now():
    ct = datetime.datetime.now()
    now = ct.strftime("%y%m%d%H%M%S")
    return now

def copy_safe(source, target, resources=[]):
    """Private function: Given a list paths it will recursively copy each to the
    target location. If a target path already exists, it will NOT over-write the
    existing paths data.
    @param resources <list[str]>:
        List of paths to copy over to target location
    @params source <str>:
        Add a prefix PATH to each resource
    @param target <str>:
        Target path to copy templates and required resources
    """

    for resource in resources:
        destination = os.path.join(target, resource)
        if not exists(destination):
            # Required resources do not exist
            shutil.copytree(os.path.join(source, resource), destination)


def sym_safe(input_data, target):
    """Creates re-named symlinks for each FastQ file provided
    as input. If a symlink already exists, it will not try to create a new symlink.
    If relative source PATH is provided, it will be converted to an absolute PATH.
    @param input_data <list[<str>]>:
        List of input files to symlink to target location
    @param target <str>:
        Target path to copy templates and required resources
    @return input_fastqs list[<str>]:
        List of renamed input FastQs
    """
    input_fastqs = []  # store renamed fastq file names
    for file in input_data:
        filename = os.path.basename(file)
        renamed = os.path.join(target, rename(filename))
        input_fastqs.append(renamed)

        if not exists(renamed):
            # Create a symlink if it does not already exist
            # Follow source symlinks to resolve any binding issues
            os.symlink(os.path.abspath(os.path.realpath(file)), renamed)

    return input_fastqs


def rename(filename):
    """Dynamically renames FastQ file to have one of the following extensions: *.R1.fastq.gz, *.R2.fastq.gz
    To automatically rename the fastq files, a few assumptions are made. If the extension of the
    FastQ file cannot be inferred, an exception is raised telling the user to fix the filename
    of the fastq files.
    @param filename <str>:
        Original name of file to be renamed
    @return filename <str>:
        A renamed FastQ filename
    """
    # Covers common extensions from SF, SRA, EBI, TCGA, and external sequencing providers
    # key = regex to match string and value = how it will be renamed
    extensions = {
        # Matches: _R[12]_fastq.gz, _R[12].fastq.gz, _R[12]_fq.gz, etc.
        ".R1.f(ast)?q.gz$": ".R1.fastq.gz",
        ".R2.f(ast)?q.gz$": ".R2.fastq.gz",
        # Matches: _R[12]_001_fastq_gz, _R[12].001.fastq.gz, _R[12]_001.fq.gz, etc.
        # Capture lane information as named group
        ".R1.(?P<lane>...).f(ast)?q.gz$": ".R1.fastq.gz",
        ".R2.(?P<lane>...).f(ast)?q.gz$": ".R2.fastq.gz",
        # Matches: _[12].fastq.gz, _[12].fq.gz, _[12]_fastq_gz, etc.
        "_1.f(ast)?q.gz$": ".R1.fastq.gz",
        "_2.f(ast)?q.gz$": ".R2.fastq.gz",
    }

    if (
        filename.endswith(".R1.fastq.gz")
        or filename.endswith(".R2.fastq.gz")
        or filename.endswith(".bam")
    ):
        # Filename is already in the correct format
        return filename

    converted = False
    for regex, new_ext in extensions.items():
        matched = re.search(regex, filename)
        if matched:
            # regex matches with a pattern in extensions
            converted = True
            filename = re.sub(regex, new_ext, filename)
            break  # only rename once

    if not converted:
        raise NameError(
            """\n\tFatal: Failed to rename provided input '{}'!
        Cannot determine the extension of the user provided input file.
        Please rename the file list above before trying again.
        Here is example of acceptable input file extensions:
          sampleName.R1.fastq.gz      sampleName.R2.fastq.gz
          sampleName_R1_001.fastq.gz  sampleName_R2_001.fastq.gz
          sampleName_1.fastq.gz       sampleName_2.fastq.gz
        Please also check that your input files are gzipped?
        If they are not, please gzip them before proceeding again.
        """.format(
                filename, sys.argv[0]
            )
        )

    return filename


def setup(sub_args, repo_path, output_path, create_nidap_folder_YN="no", links=[]):
    """Setup the pipeline for execution and creates config file from templates
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param repo_path <str>:
        Path to RNA-seek source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @param create_nidap_folder_YN <str>:
        yes or no
    @return config <dict>:
         Config dictionary containing metadata to run the pipeline
    """
    # Check for mixed inputs,
    # inputs which are a mixture
    # of FastQ and BAM files
    ifiles = sym_safe(input_data=links, target=output_path)
    mixed_inputs(ifiles)

    shorthostname = get_hpcname()
    if not shorthostname:
        shorthostname = "biowulf"
        print(
            f"{shorthostname} unknown host. Configuration files for references may not be correct. Defaulting to Biowulf config"
        )
    else:
        print(f"Thank you for running XAVIER on {shorthostname.upper()}")

    genome_config = os.path.join(
        repo_path, "config", "genomes", get_hpcname(), sub_args.genome + ".json"
    )

    if sub_args.genome.endswith(".json"):
        # Provided a custom reference genome generated by rna-seek build
        genome_config = os.path.abspath(sub_args.genome)

    if not os.path.exists(genome_config):
        raise FileNotFoundError(f"Genome config file does not exist: {genome_config}")

    required = {
        # Base configuration file
        "base": os.path.join(repo_path, "config", "config.json"),
        # Template for project-level information
        "project": os.path.join(repo_path, "config", "templates", "project.json"),
        # Template for genomic reference files
        # User provided argument --genome is used to select the template
        "genome": genome_config,
        # Template for tool information
        "tools": os.path.join(repo_path, "config", "templates", "tools.json"),
    }

    cluster_config = os.path.join(
        repo_path, "config", "cluster" + "." + shorthostname + ".json"
    )
    cluster_output = os.path.join(output_path, "cluster.json")
    shutil.copyfile(cluster_config, cluster_output)

    # Global config file for pipeline, config.json
    config = join_jsons(required.values())  # uses templates in the rna-seek repo
    config = add_user_information(config)
    config = add_rawdata_information(sub_args, config, ifiles)

    # Resolves if an image needs to be pulled from an OCI registry or
    # a local SIF generated from the rna-seek cache subcommand exists
    config = image_cache(sub_args, config)

    # Add other cli collected info
    config["project"]["annotation"] = sub_args.genome
    config["project"]["version"] = get_version()
    config["project"]["workpath"] = os.path.abspath(sub_args.output)

    # Add optional cli workflow steps
    config["input_params"]["CNV_CALLING"] = str(sub_args.cnv).lower()
    config["input_params"]["FFPE_FILTER"] = str(sub_args.ffpe).lower()
    config["input_params"]["EXOME_TARGETS"] = (
        str(sub_args.targets)
        if sub_args.targets
        else os.path.join(
            config["project"]["workpath"], config["references"]["exome_targets"]
        )
    )
    if not os.path.exists(config["input_params"]["EXOME_TARGETS"]):
        raise FileNotFoundError(
            f"Exome targets file does not exist: {config['input_params']['EXOME_TARGETS']}"
        )
    config["input_params"]["VARIANT_CALLERS"] = sub_args.callers
    config["input_params"]["PAIRS_FILE"] = str(sub_args.pairs)
    config["input_params"]["BASE_OUTDIR"] = str(sub_args.output)
    config["input_params"]["tmpdisk"] = str(sub_args.tmp_dir)
    config["input_params"]["create_nidap_folder"] = str(create_nidap_folder_YN)

    # Get latest git commit hash
    git_hash = git_commit_hash(repo_path)
    config["project"]["git_commit_hash"] = git_hash
    config["project"]["pipehome"] = repo_path
    config["project"]["pairs"] = str(sub_args.pairs)

    if sub_args.runmode == "init":
        # Save config to output directory
        with open(os.path.join(output_path, "config.json"), "w") as fh:
            json.dump(config, fh, indent=4, sort_keys=True)

    return config


def unpacked(nested_dict):
    """Generator to recursively retrieves all values in a nested dictionary.
    @param nested_dict dict[<any>]:
        Nested dictionary to unpack
    @yields value in dictionary
    """
    # Iterate over all values of given dictionary
    for value in nested_dict.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over
            # all its values recursively
            for v in unpacked(value):
                yield v
        else:
            # If value is not dict type then
            # yield the value
            yield value


def get_fastq_screen_paths(fastq_screen_confs, match="DATABASE", file_index=-1):
    """Parses fastq_screen.conf files to get the paths of each fastq_screen database.
    This path contains bowtie2 indices for reference genome to screen against.
    The paths are added as singularity bind points.
    @param fastq_screen_confs list[<str>]:
        Name of fastq_screen config files to parse
    @param match <string>:
        Keyword to indicate a line match [default: 'DATABASE']
    @param file_index <int>:
        Index of line line containing the fastq_screen database path
    @return list[<str>]:
        Returns a list of fastq_screen database paths
    """
    databases = []
    for file in fastq_screen_confs:
        with open(file, "r") as fh:
            for line in fh:
                if line.startswith(match):
                    db_path = line.strip().split()[file_index]
                    databases.append(db_path)
    return databases


def resolve_additional_bind_paths(search_paths):
    """Finds additional singularity bind paths from a list of random paths. Paths are
    indexed with a compostite key containing the first two directories of an absolute
    file path to avoid issues related to shared names across the /gpfs shared network
    filesystem. For each indexed list of file paths, a common path is found. Assumes
    that the paths provided are absolute paths, the build sub command creates reference
    files with absolute filenames.
    @param search_paths list[<str>]:
        List of absolute file paths to find common bind paths from
    @return common_paths list[<str>]:
        Returns a list of common shared file paths to create additional singularity bind paths
    """
    common_paths = []
    indexed_paths = {}

    for ref in search_paths:
        # Skip over resources with remote URI and
        # skip over strings that are not file PATHS as
        # build command creates absolute resource PATHS
        if (
            ref.lower().startswith("sftp://")
            or ref.lower().startswith("s3://")
            or ref.lower().startswith("gs://")
            or not ref.lower().startswith(os.sep)
        ):
            continue

        # Break up path into directory tokens
        path_list = os.path.abspath(ref).split(os.sep)
        try:  # Create composite index from first two directories
            # Avoids issues created by shared /gpfs/ PATHS
            index = path_list[1:3]
            index = tuple(index)
        except IndexError:
            index = path_list[1]  # ref startswith /
        if index not in indexed_paths:
            indexed_paths[index] = []
        # Create an INDEX to find common PATHS for each root
        # child directory like /scratch or /data. This prevents
        # issues when trying to find the common path between
        # these two different directories (resolves to /)
        indexed_paths[index].append(str(os.sep).join(path_list))

    for index, paths in indexed_paths.items():
        # Find common paths for each path index
        p = os.path.dirname(os.path.commonprefix(paths))
        if p == os.sep:
            # Avoids adding / to bind list when
            # given /tmp or /scratch as input
            p = os.path.commonprefix(paths)
        common_paths.append(p)

    return list(set(common_paths))


def bind(sub_args, config):
    """Resolves bindpaths for singularity/docker images.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param configfile dict[<any>]:
        Config dictionary generated by setup command.
    @return bindpaths list[<str>]:
        List of singularity/docker bind paths
    """
    bindpaths = []
    for value in unpacked(config):
        if not isinstance(value, str):
            continue
        if exists(value):
            if os.path.isfile(value):
                value = os.path.dirname(value)
            if value not in bindpaths:
                bindpaths.append(value)

    # Get FastQ Screen Database paths
    # and other reference genome file paths
    rawdata_bind_paths = [
        os.path.realpath(p) for p in config["project"]["datapath"].split(",")
    ]
    working_directory = os.path.realpath(config["project"]["workpath"])
    fqscreen_cfg = config["references"]["FASTQ_SCREEN_CONFIG"]
    fq_screen_paths = get_fastq_screen_paths(
        [os.path.join(sub_args.output, fqscreen_cfg)]
    )
    kraken_db_path = config["references"]["KRAKENBACDB"]
    # Add Bindpath for VCF2maf
    vep_db_path = config["references"]["VCF2MAF"]["VEPRESOURCEBUNDLEPATH"]
    genome_bind_paths = resolve_additional_bind_paths(
        bindpaths + fq_screen_paths + [kraken_db_path] + [vep_db_path]
    )
    bindpaths = [working_directory] + rawdata_bind_paths + genome_bind_paths
    bindpaths = list(set([p for p in bindpaths if p != os.sep]))

    return bindpaths


def mixed_inputs(ifiles):
    """Check if a user has provided a set of input files which contain a
    mixture of FastQ and BAM files. The pipeline does not support processing
    a mix of FastQ and BAM files.
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    """
    bam_files, fq_files = [], []
    fastqs = False
    bams = False
    for file in ifiles:
        if file.endswith(".R1.fastq.gz") or file.endswith(".R2.fastq.gz"):
            fastqs = True
            fq_files.append(file)
        elif file.endswith(".bam"):
            bams = True
            bam_files.append(file)

    if fastqs and bams:
        # User provided a mix of FastQs and BAMs
        raise TypeError(
            """\n\tFatal: Detected a mixture of --input data types.
            A mixture of BAM and FastQ files were provided; however, the pipeline
            does NOT support processing a mixture of input FastQ and BAM files.
            Input FastQ Files:
                {}
            Input BAM Files:
                {}
            Please do not run the pipeline with a mixture of FastQ and BAM files.
            This feature is currently not supported within '{}', and it is not
            recommended to process samples in this way either. If this is a priority
            for your project, please run the set of FastQ and BAM files separately
            (in two separate output directories). If you feel like this functionality
            should exist, feel free to open an issue on Github.
            """.format(
                " ".join(fq_files), " ".join(bam_files), sys.argv[0]
            )
        )


def add_user_information(config):
    """Adds username and user's home directory to config.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """
    # Get PATH to user's home directory
    # Method is portable across unix-like OS and Windows
    home = os.path.expanduser("~")

    # Get username from home directory PATH
    username = os.path.split(home)[-1]

    # Update config with home directory and username
    config["project"]["userhome"] = home
    config["project"]["username"] = username

    return config


def add_rawdata_information(sub_args, config, ifiles):
    """Adds information about rawdata provided to pipeline.
    Determines whether the dataset is paired-end or single-end and finds the set of all
    rawdata directories (needed for -B option when running singularity). If a user provides
    paired-end data, checks to see if both mates (R1 and R2) are present for each sample.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """

    # Determine whether dataset is paired-end or single-ends
    # Updates config['project']['nends']: 1 = single-end, 2 = paired-end, -1 = bams
    convert = {1: "single-end", 2: "paired-end", -1: "bam"}
    nends = get_nends(ifiles)  # Checks PE data for both mates (R1 and R2)
    config["project"]["nends"] = nends
    config["project"]["filetype"] = convert[nends]

    # Finds the set of rawdata directories to bind
    rawdata_paths = get_rawdata_bind_paths(input_files=sub_args.input)
    config["project"]["datapath"] = ",".join(rawdata_paths)

    return config


def get_nends(ifiles):
    """Determines whether the dataset is paired-end or single-end.
    If paired-end data, checks to see if both mates (R1 and R2) are present for each sample.
    If single-end, nends is set to 1. Else if paired-end, nends is set to 2.
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    @return nends_status <int>:
         Integer reflecting nends status: 1 = se, 2 = pe, -1 = bams
    """
    # Determine if dataset contains paired-end data
    paired_end = False
    bam_files = False
    nends_status = 1
    for file in ifiles:
        if file.endswith(".bam"):
            bam_files = True
            nends_status = -1
            break
        elif file.endswith(".R2.fastq.gz"):
            paired_end = True
            nends_status = 2
            break  # dataset is paired-end

    # Check to see if both mates (R1 and R2)
    # are present paired-end data
    if paired_end:
        nends = {}  # keep count of R1 and R2 for each sample
        for file in ifiles:
            # Split sample name on file extension
            sample = re.split("\.R[12]\.fastq\.gz", os.path.basename(file))[0]
            if sample not in nends:
                nends[sample] = 0

            nends[sample] += 1

        # Check if samples contain both read mates
        missing_mates = [sample for sample, count in nends.items() if count == 1]
        if missing_mates:
            # Missing an R1 or R2 for a provided input sample
            raise NameError(
                """\n\tFatal: Detected pair-end data but user failed to provide
                both mates (R1 and R2) for the following samples:\n\t\t{}\n
                Please check that the basename for each sample is consistent across mates.
                Here is an example of a consistent basename across mates:
                  consistent_basename.R1.fastq.gz
                  consistent_basename.R2.fastq.gz

                Please do not run the pipeline with a mixture of single-end and paired-end
                samples. This feature is currently not supported within {}, and it is
                not recommended either. If this is a priority for your project, please run
                paired-end samples and single-end samples separately (in two separate output
                directories). If you feel like this functionality should exist, feel free to
                open an issue on Github.
                """.format(
                    missing_mates, sys.argv[0]
                )
            )
    elif not bam_files:
        # Provided only single-end data
        # not supported or recommended
        raise TypeError(
            """\n\tFatal: Single-end data detected.
            {} does not support single-end data. Calling variants from single-end
            data is not recommended either. If you feel like this functionality should
            exist, feel free to open an issue on Github.
            """.format(
                sys.argv[0]
            )
        )

    return nends_status


def get_rawdata_bind_paths(input_files):
    """
    Gets rawdata bind paths of user provided fastq files.
    @params input_files list[<str>]:
        List containing user-provided input fastq files
    @return bindpaths <set>:
        Set of rawdata bind paths
    """
    bindpaths = []
    for file in input_files:
        # Get directory of input file
        rawdata_src_path = os.path.dirname(os.path.abspath(os.path.realpath(file)))
        if rawdata_src_path not in bindpaths:
            bindpaths.append(rawdata_src_path)

    return bindpaths


def dryrun(
    outdir, config="config.json", snakefile=os.path.join("workflow", "Snakefile"),
     write_to_file=True,
):
    """Dryruns the pipeline to ensure there are no errors prior to running.
    @param outdir <str>:
        Pipeline output PATH
    @return dryrun_output <str>:
        Byte string representation of dryrun command
    """
    try:
        dryrun_output = subprocess.check_output(
            [
                "snakemake",
                "-npr",
                "--rerun-incomplete",
                "-s",
                str(snakefile),
                "--use-singularity",
                "--cores",
                str(1),
                "--configfile={}".format(config),
            ],
            cwd=outdir,
            stderr=subprocess.STDOUT,
        )
    except OSError as e:
        # Catch: OSError: [Errno 2] No such file or directory
        #  Occurs when command returns a non-zero exit-code
        if e.errno == 2 and not which("snakemake"):
            # Failure caused because snakemake is NOT in $PATH
            err(
                "\n\x1b[6;37;41mError: Are snakemake AND singularity in your $PATH?\x1b[0m"
            )
            fatal("\x1b[6;37;41mPlease check before proceeding again!\x1b[0m")
        else:
            # Failure caused by unknown cause, raise error
            raise e
    except subprocess.CalledProcessError as e:
        print(e, e.output)
        raise (e)
    
    if write_to_file:
        now = _now()
        with open(os.path.join(outdir, "dryrun." + str(now) + ".log"), "w") as outfile:
            outfile.write("{}".format(dryrun_output.decode("utf-8")))

    return dryrun_output


def runner(
    mode,
    outdir,
    alt_cache,
    logger,
    additional_bind_paths=None,
    threads=2,
    jobname="pl:xavier",
    submission_script="runner",
    tmp_dir="/lscratch/$SLURM_JOBID/",
    wait="",
):
    """Runs the pipeline via selected executor: local or slurm.
    If 'local' is selected, the pipeline is executed locally on a compute node/instance.
    If 'slurm' is selected, jobs will be submitted to the cluster using SLURM job scheduler.
    Support for additional job schedulers (i.e. PBS, SGE, LSF) may be added in the future.
    @param outdir <str>:
        Pipeline output PATH
    @param mode <str>:
        Execution method or mode:
            local runs serially a compute instance without submitting to the cluster.
            slurm will submit jobs to the cluster using the SLURM job scheduler.
    @param additional_bind_paths <str>:
        Additional paths to bind to container filesystem (i.e. input file paths)
    @param alt_cache <str>:
        Alternative singularity cache location
    @param logger <file-handle>:
        An open file handle for writing
    @param threads <str>:
        Number of threads to use for local execution method
    @param masterjob <str>:
        Name of the master job
    @param wait <str>:
        "--wait" or "" ... used only while submitting job via HPC API
    @return masterjob <subprocess.Popen() object>:
    """
    # Add additional singularity bind PATHs
    # to mount the local filesystem to the
    # containers filesystem, NOTE: these
    # PATHs must be an absolute PATHs
    outdir = os.path.abspath(outdir)
    # Add any default PATHs to bind to
    # the container's filesystem, like
    # tmp directories, /lscratch
    bindpaths = "{},{}".format(outdir, os.path.dirname(tmp_dir.rstrip("/")))
    # Set ENV variable 'SINGULARITY_CACHEDIR'
    # to output directory
    my_env = {}
    my_env.update(os.environ)
    cache = os.path.join(outdir, ".singularity")
    my_env["SINGULARITY_CACHEDIR"] = cache
    if alt_cache:
        # Override the pipeline's default
        # cache location
        my_env["SINGULARITY_CACHEDIR"] = alt_cache
        cache = alt_cache

    if additional_bind_paths:
        # Add Bind PATHs for rawdata directories
        bindpaths = "{},{}".format(additional_bind_paths, bindpaths)

    if not exists(os.path.join(outdir, "logfiles")):
        # Create directory for logfiles
        os.makedirs(os.path.join(outdir, "logfiles"))

    # Create .singularity directory for
    # installations of snakemake without
    # setuid which creates a sandbox in
    # the SINGULARITY_CACHEDIR
    if not exists(cache):
        # Create directory for sandbox
        # and image layers
        os.makedirs(cache)

    # Run on compute node or instance
    # without submitting jobs to a scheduler
    if mode == "local":
        # Run pipeline's main process
        # Look into later: it maybe worth
        # replacing Popen subprocess with a direct
        # snakemake API call: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html
        masterjob = subprocess.Popen(
            [
                "snakemake",
                "-pr",
                "--rerun-incomplete",
                "--use-singularity",
                "--singularity-args",
                "'-B {}'".format(bindpaths),
                "--cores",
                str(threads),
                "--configfile=config.json",
            ],
            cwd=outdir,
            stderr=subprocess.STDOUT,
            stdout=logger,
            env=my_env,
        )

    # Submitting jobs to cluster via SLURM's job scheduler
    elif mode == "slurm":
        # Run pipeline's main process
        # Look into later: it maybe worth
        # replacing Popen subprocess with a direct
        # snakemake API call: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html
        # CLUSTER_OPTS="'sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} \
        #   -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e $SLURMDIR/slurm-%j_{params.rname}.out \
        #   -o $SLURMDIR/slurm-%j_{params.rname}.out'"
        # sbatch --parsable -J "$2" --gres=lscratch:500 --time=5-00:00:00 --mail-type=BEGIN,END,FAIL \
        #   --cpus-per-task=32 --mem=96g --output "$3"/logfiles/snakemake.log --error "$3"/logfiles/snakemake.log \
        # snakemake --latency-wait 120 -s "$3"/workflow/Snakefile -d "$3" \
        #   --use-singularity --singularity-args "'-B $4'" --configfile="$3"/config.json \
        #   --printshellcmds --cluster-config "$3"/resources/cluster.json \
        #   --cluster "${CLUSTER_OPTS}" --keep-going --restart-times 3 -j 500 \
        #   --rerun-incomplete --stats "$3"/logfiles/runtime_statistics.json \
        #   --keep-remote --local-cores 30 2>&1 | tee -a "$3"/logfiles/master.log
        if wait == "":
            masterjob = subprocess.Popen(
                [
                    str(os.path.join(outdir, "resources", str(submission_script))),
                    mode,
                    "-j",
                    jobname,
                    "-b",
                    str(bindpaths),
                    "-o",
                    str(outdir),
                    "-c",
                    str(cache),
                    "-t",
                    "'{}'".format(tmp_dir),
                ],
                cwd=outdir,
                stderr=subprocess.STDOUT,
                stdout=logger,
                env=my_env,
            )
        else:
            masterjob = subprocess.Popen(
                [
                    str(os.path.join(outdir, "resources", str(submission_script))),
                    mode,
                    "-j",
                    jobname,
                    "-b",
                    str(bindpaths),
                    "-o",
                    str(outdir),
                    "-c",
                    str(cache),
                    str(wait),
                    "-t",
                    "'{}'".format(tmp_dir),
                ],
                cwd=outdir,
                stderr=subprocess.STDOUT,
                stdout=logger,
                env=my_env,
            )

    return masterjob
