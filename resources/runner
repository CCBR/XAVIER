#!/usr/bin/env bash
set -eu

function usage() { cat << EOF
runner: submits the master job of pipeline to the job scheduler.

USAGE:
  runner <MODE> [-h] [-c CACHE] \\
          -o OUTDIR \\
          -j MASTER_JOB_NAME \\
          -b SINGULARITY_BIND_PATHS \\
          -t TMP_DIR

SYNOPSIS:
  This script creates/submits the pipeline's master job  to  the
cluster. The master job acts as the pipeline's main controller or
its main process. This main job dictates how subsequent jobs are
submitted to the cluster via the job scheduler, SLURM. Support for
additional job schedulers (i.e. PBS, SGE, LSF, Tibanna) may be added
in the future.
  The main entry point of the pipeline 'xavier' calls this job
submission wrapper script. As so, this script can be used to by-pass
'xavier' for a previously failed run; meaning, it can be used to
re-run the pipeline to pick back off where the last failure occurred
or re-start the pipeline.
  Please Note: it is highly recommended to use 'xavier'; it is the
main entry point and preferred entry point of the xavier pipeline.
If you are experiencing an error, it maybe due to improperly mounting
singularity bind paths which 'xavier' will internally handle.

Required Positional Argument:
  [1] MODE  [Type: Str] Defines the snakemake execution mode.
                        Valid mode options include: <slurm, ...>
                          slurm: uses slurm and singularity/snakemake
                             backend. This EXECUTOR will submit child
                             jobs to the cluster. It is recommended
                             running xavier in this mode as most
                             of the steps are computationally intensive.

Required Arguments:
  -o,  --outdir  [Type: Path]   Path to output directory of the pipeline.
                                 This is the pipeline's working directory
                                 where all output files will be generated.
  -j, --job-name [Type: Str]    Name of pipeline's master job.
  -b, --bind-paths [Type:Path]  Singularity bind paths. The pipeline uses
                                 singularity images for execution. Bind
                                 paths are used to mount the host file
                                 system to the container's file system.
                                 Multiple bind paths can be provided as
                                 a comma-separated list. The main entry
                                 point of the pipeline internally collects
                                 and aggregates bind paths to mount to the
                                 container's filesystem.
                                 If you are manually running this script
                                 or by-passing xavier, you will need
                                 to provide the bind paths of the rawdata
                                 directory(s) along with the pipeline's
                                 output directory and any directories for
                                 reference files. Please see example usage
                                 below.
 -t, --tmp-dir [Type:Path]      Temporary directory. The pipeline generates
                                 intermediate, temporary output files. Any
                                 temporary output files will be written to
                                 this location. On Biowulf, it should be
                                 set to '/lscratch/\$SLURM_JOBID/'. On FRCE,
                                 this value should be set to the following:
                                 '/scratch/cluster_scratch/\$USER/'.

OPTIONS:
  -c,  --cache  [Type: Path]   Path to singularity cache. If not provided,
                                the path will default to the current working
                                directory of this script.
                                [Default: $(dirname  "$0")/.singularity/]
  -h, --help    [Type: Bool]  Displays usage and help information.
  -w, --wait    [Type: Bool]  Wait until master job completes. This is required if
                                the job is submitted using HPC API. If not provided
                                the API may interpret submission of master job as
                                completion of the pipeline!

Example:
  $ runner slurm -h
  $ runner slurm -j xavier_hg38 -b "/data/$USER/,/lscratch"

Version:
  1.0.0
EOF
}


# Functions
function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; usage; exit 1; }
function abspath() { readlink -e "$1"; }
function parser() {
  # Adds parsed command-line args to GLOBAL $Arguments associative array
  # + KEYS = short_cli_flag ("j", "o", ...)
  # + VALUES = parsed_user_value ("MasterJobName" "/scratch/hg38", ...)
  # @INPUT "$@" = user command-line arguments
  # @CALLS check() to see if the user provided all the required arguments

  while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
      -h  | --help) usage && exit 0;;
      -j  | --job-name)   provided "$key" "${2:-}"; Arguments["j"]="$2"; shift; shift;;
      -b  | --bind-paths) provided "$key" "${2:-}"; Arguments["b"]="$2"; shift; shift;;
      -t  | --tmp-dir) provided "$key" "${2:-}"; Arguments["t"]="$2"; shift; shift;;
      -o  | --outdir)  provided "$key" "${2:-}"; Arguments["o"]="$2"; shift; shift;;
      -c  | --cache)  provided "$key" "${2:-}"; Arguments["c"]="$2"; shift; shift;;
      -w  | --wait)  Arguments["w"]="--wait"; shift;;
      -*  | --*) err "Error: Failed to parse unsupported argument: '${key}'."; usage && exit 1;;
      *) err "Error: Failed to parse unrecognized argument: '${key}'. Do any of your inputs have spaces?"; usage && exit 1;;
    esac
  done

  # Check for required args
  check
}


function provided() {
  # Checks to see if the argument's value exists
  # @INPUT $1 = name of user provided argument
  # @INPUT $2 = value of user provided argument
  # @CALLS fatal() if value is empty string or NULL

  if [[ -z "${2:-}" ]]; then
     fatal "Fatal: Failed to provide value to '${1}'!";
  fi
}


function check(){
  # Checks to see if user provided required arguments
  # @INPUTS $Arguments = Global Associative Array
  # @CALLS fatal() if user did NOT provide all the $required args

  # List of required arguments
  local required=("j" "b" "o")
  #echo -e "Provided Required Inputs"
  for arg in "${required[@]}"; do
    value=${Arguments[${arg}]:-}
    if [[ -z "${value}" ]]; then
      fatal "Failed to provide all required args.. missing ${arg}"
    fi
  done
}


function require(){
  # Requires an executable is in $PATH, as a last resort it will attempt to load
  # the executable or dependency as a module
  # @INPUT $@ = List of dependencies or executables to check

  for exe in "${@}"; do
    # Check if executable is in $PATH
    command -V ${exe} &> /dev/null && continue;
    # Try to load exe as lua module
    module load ${exe} &> /dev/null || \
      fatal "Failed to find or load '${exe}', not installed on target system."
  done
}


function submit(){
  # Submit jobs to the defined job scheduler or executor (i.e. slurm)
  # INPUT $1 = Snakemake Mode of execution
  # INPUT $2 = Name of master/main job or process (pipeline controller)
  # INPUT $3 = Pipeline output directory
  # INPUT $4 = Singularity Bind paths
  # INPUT $5 = Singularity cache directory
  # INPUT $6 = Temporary directory for output files
  # INPUT $7 = Wait ("--wait") or no wait ("--nowait".. default)

  # Check if singularity and snakemake are in $PATH
  # If not, try to module load singularity as a last resort
  require singularity snakemake

  # Snakemake executor, or target job scheduler
  # more maybe added in the future, TBA
  executor=${1}

  # Goto Pipeline Output directory
  # Create a local singularity cache in output directory
  # cache can be re-used instead of re-pulling from DockerHub every time
  cd "$3" && export SINGULARITY_CACHEDIR="${5}"

  # unsetting XDG_RUNTIME_DIR to avoid some unsighly but harmless warnings
  unset XDG_RUNTIME_DIR

  # Run the workflow with specified executor
  case "$executor" in
    slurm)
          # Create directory for logfiles
          mkdir -p "$3"/logfiles/slurmfiles/
          CLUSTER_OPTS="sbatch --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"
          snakemakeVer=$(snakemake --version 2>/dev/null)
          # Submit the master job to the cluster
          # sbatch --parsable -J {jobname} --time=5-00:00:00 --mail-type=BEGIN,END,FAIL
          # --cpus-per-task=24 --mem=96g --gres=lscratch:500
          # --output os.path.join({outdir}, 'logfiles', 'snakemake.log') --error os.path.join({outdir}, 'logfiles', 'snakemake.log')
          # snakemake -pr --latency-wait 120 -d {outdir} --configfile=config.json
          # --cluster-config os.path.join({outdir}, 'config', 'cluster.json')
          # --cluster {CLUSTER_OPTS} --stats os.path.join({outdir}, 'logfiles', 'runtime_statistics.json')
          # --printshellcmds --keep-going --rerun-incomplete
          # --keep-remote --restart-times 3 -j 500 --use-singularity
          # --singularity-args -B {}.format({bindpaths}) --local-cores 24
          # Check if NOT running on Biowulf
          # Assumes other clusters do NOT
          # have GRES for local node disk,
          # long term it might be worth
          # adding a new option to allow
          # a user to decide whether to
          # use GRES at job submission,
          # trying to infer this because
          # most users will not even know
          # what GRES is and how or why
          # it should be used and by default
          # SLURM is not configured to use
          # GRES, remove prefix single quote
          #if [[ ${6#\'} != /lscratch* ]]; then
          #  CLUSTER_OPTS="sbatch --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e $SLURM_DIR/slurm-%j_{params.rname}.out -o $SLURM_DIR/slurm-%j_{params.rname}.out"
          #fi
          verlte() {
        [  "$1" = "`echo -e "$1\n$2" | sort -V | head -n1`" ]
          }
          verlt() {
        [ "$1" = "$2" ] && return 1 || verlte $1 $2
          }
          snakemakeOld=$(verlt $snakemakeVer 7.8 && echo "yes" || echo "no") # check if snakemake is older than 7.8
          if [ "$snakemakeOld" == "no" ];then
            triggeroptions="--rerun-triggers mtime"
          else
            triggeroptions=""
          fi
          if [[ "$HOSTNAME" == "biowulf.nih.gov" ]];then
            CLUSTER_OPTS="${CLUSTER_OPTS} --gres {cluster.gres}"
          elif [[ "$HOSTNAME" == cn[0-9][0-9][0-9][0-9] ]];then
            CLUSTER_OPTS="${CLUSTER_OPTS} --gres {cluster.gres}"
          fi
          # Create sbatch script to build index
    cat << EOF > kickoff.sh
#!/usr/bin/env bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=96g
#SBATCH --time=5-00:00:00
#SBATCH --parsable
#SBATCH -J "$2"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output "$3/logfiles/snakemake.log"
#SBATCH --error "$3/logfiles/snakemake.log"

set -euo pipefail

snakemake --latency-wait 120 -s "$3/workflow/Snakefile" -d "$3" \\
  --use-singularity --singularity-args "'-B $4'" --configfile="$3/config.json" \\
  --printshellcmds --cluster-config "$3/cluster.json" \\
  --cluster "${CLUSTER_OPTS}" --keep-going --restart-times 3 -j 500 \\
  $triggeroptions \\
  --rerun-incomplete --stats "$3/logfiles/runtime_statistics.json" \\
  --keep-remote --local-cores 30 2>&1
EOF
    chmod +x kickoff.sh
    wait="$7"
    if [ "$wait" != "--wait" ];then
      job_id=$(sbatch kickoff.sh | tee -a "$3"/logfiles/master.log)
    else
      job_id=$(sbatch --wait kickoff.sh | tee -a "$3"/logfiles/master.log)
    fi
        ;;
      *)  echo "${executor} is not available." && \
          fatal "Failed to provide valid execution backend: ${executor}. Please use slurm."
        ;;
    esac

  # Return exit-code of pipeline sumbission
  echo "$job_id"
}


function main(){
  # Parses args and submits master job of pipeline to the target job scheduler
  # @INPUT "$@" = command-line arguments
  # @CALLS parser(), initialize(), setup(), cromwell()

  if [ $# -eq 0 ]; then usage; exit 1; fi

  # Associative array to store parsed args
  declare -Ag Arguments

  # Positional Argument for Snakemake Executor
  case $1 in
    slurm) Arguments["e"]="$1";;
    -h    | --help | help) usage && exit 0;;
    -*    | --*) err "Error: Failed to provide required positional argument: <slurm>."; usage && exit 1;;
    *) err "Error: Failed to provide valid positional argument. '${1}' is not supported. Valid option(s) are slurm"; usage && exit 1;;
  esac

  # Parses remaining user provided command-line arguments
  parser "${@:2}" # Remove first item of list
  outdir="$(abspath "$(dirname  "${Arguments[o]}")")"
  Arguments[o]="${Arguments[o]%/}" # clean outdir path (remove trailing '/')

  # Setting defaults for non-required arguments
  # If singularity cache not provided, default to ${outdir}/.singularity
  cache="${Arguments[o]}/.singularity"
  Arguments[c]="${Arguments[c]:-$cache}"
  Arguments[c]="${Arguments[c]%/}" # clean outdir path (remove trailing '/')

  default_wait="--nowait"
  Arguments[w]="${Arguments[w]:-$default_wait}" # wait is a blank string if not provided.

  # Print pipeline metadata prior to running
  echo -e "xavier\t$(date)"
  echo -e "Running pipeline with the following parameters:"
  for key in "${!Arguments[@]}"; do echo -e "\t${key}\t${Arguments["$key"]}"; done

  # Run pipeline and submit jobs to cluster using the defined executor
  mkdir -p "${Arguments[o]}/logfiles/"
  job_id=$(submit "${Arguments[e]}" "${Arguments[j]}" "${Arguments[o]}" "${Arguments[b]}" "${Arguments[c]}" "${Arguments[t]}" "${Arguments[w]}")
  echo -e "XAVIER pipeline submitted to cluster.\nMaster Job ID: $job_id"
  echo "${job_id}" > "${Arguments[o]}/logfiles/mjobid.log"

}


# Main: check usage, parse args, and run pipeline
main "$@"
