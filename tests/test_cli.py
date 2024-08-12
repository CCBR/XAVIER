import json
import os
import subprocess
import tempfile
from xavier.src.xavier.__main__ import main

xavier_run = (
    "xavier run "
    "--input .tests/*.fastq.gz "
    "--pairs .tests/pairs.tsv "
    "--mode local "
)


def run_in_temp(command_str):
    with tempfile.TemporaryDirectory() as tmp_dir:
        outdir = os.path.join(tmp_dir, "testout")
        run_command = f"{command_str} --output {outdir}"
        output = subprocess.run(
            f"{run_command} --runmode init && {run_command} --runmode dryrun",
            capture_output=True,
            shell=True,
            text=True,
        )
        if os.path.exists(os.path.join(outdir, "config.json")):
            with open(os.path.join(outdir, "config.json"), "r") as infile:
                config = json.load(infile)
        else:
            config = None
    return output, config


def test_help():
    assert (
        "XAVIER"
        in subprocess.run(
            "./bin/xavier --help", capture_output=True, shell=True, text=True
        ).stdout
    )


def test_dryrun_targets():
    output_human, config_human = run_in_temp(f"{xavier_run} --genome hg38")
    output_mouse, config_mouse = run_in_temp(f"{xavier_run} --genome mm10")
    output_custom, config_custom = run_in_temp(
        f"{xavier_run} --genome mm10 --targets resources/Agilent_SSv7_allExons_hg38.bed"
    )
    output_invalid, config_invalid = run_in_temp(
        f"{xavier_run} --genome hg38 --target not/a/file.txt"
    )
    assert all(
        [
            "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
            in output_human.stdout,
            "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
            in output_mouse.stdout,
            "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
            in output_custom.stdout,
            "error: Path 'not/a/file.txt' does not exists! Failed to provide valid input."
            in output_invalid.stderr,
            config_human["input_params"]["EXOME_TARGETS"].endswith(
                "resources/Agilent_SSv7_allExons_hg38.bed"
            ),
            config_mouse["input_params"]["EXOME_TARGETS"].endswith(
                "resources/SureSelect_mm10_sorted.bed"
            ),
            config_custom["input_params"]["EXOME_TARGETS"].endswith(
                "resources/Agilent_SSv7_allExons_hg38.bed"
            ),
            not config_invalid,
        ]
    )
