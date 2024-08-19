#!/usr/bin/env python3
import argparse
import os
import sys
import glob
import PySimpleGUI as sg

from ccbr_tools.pipeline.util import (
    get_genomes_dict,
    get_tmp_dir,
    get_hpcname,
    check_python_version,
)
from ccbr_tools.pipeline.cache import get_sif_cache_dir
from ccbr_tools.shell import exec_in_context

from .util import xavier_base, get_version
from .run import run


def launch_gui(DEBUG=True):
    check_python_version()
    # get drop down genome options
    jsons = get_genomes_dict(repo_base=xavier_base)
    genome_annotation_combinations = list(jsons.keys())
    genome_annotation_combinations.sort()
    if DEBUG:
        print(jsons)
    if DEBUG:
        print(genome_annotation_combinations)

    # Create different layouts
    tumorPair_layout = [
        [
            sg.Text("Pairs file", size=(20, 1)),
            sg.InputText(key="-PAIRS-"),
            sg.FileBrowse(target="-PAIRS-"),
        ],
        [
            sg.Text("Copy Number Variants (CNV):"),
            sg.Radio("No", "CNVRADIO", enable_events=True, default=True, key="-NOCNV-"),
            sg.Radio("Yes", "CNVRADIO", enable_events=True, key="-CNV-"),
        ],
    ]

    tumorOnly_layout = [
        [sg.T("Copy Number Variants (CNVs) can only be analyzed in Tumor-Normal mode.")]
    ]

    analysis_layout = [
        [
            sg.Radio(
                "Tumor-normal pair", "TUMORADIO", enable_events=True, key="-TUMNORM-"
            ),
            sg.Radio("Tumor-only", "TUMORADIO", enable_events=True, key="-TUMONLY-"),
        ],
        [
            sg.Frame(
                "Tumor-Normal",
                tumorPair_layout,
                font=("Helvetica", 12, "bold"),
                key="-PAIROPTS-",
                visible=False,
            )
        ],
        [
            sg.Frame(
                "Tumor-Only",
                tumorOnly_layout,
                font=("Helvetica", 12, "bold"),
                key="-ONLYOPTS-",
                visible=False,
            )
        ],
    ]

    targets_layout = [
        [
            sg.Text("Targets .BED file", size=(20, 1)),
            sg.InputText(key="-TARGETS-"),
            sg.FileBrowse(target="-TARGETS-"),
        ]
    ]

    settings_layout = [
        [sg.T("Please read the Documentation before changing")],
        [
            sg.T("Apply FFPE correction?"),
            sg.Radio(
                "No", "FFPERADIO", enable_events=True, default=True, key="-NOFFPE-"
            ),
            sg.Radio("Yes", "FFPERADIO", enable_events=True, key="-FFPE-"),
        ],
        [sg.T("Targets (.BED file):")],
        [
            sg.Radio(
                "Default", "BEDRADIO", enable_events=True, default=True, key="-DEFTARG-"
            ),
            sg.Radio("Custom", "BEDRADIO", enable_events=True, key="-CUSTARG-"),
        ],
        [sg.Frame("Custom Targets", targets_layout, visible=False, key="-BED-")],
        [sg.Button(button_text="Discard", key="-DISCSET-", button_color="#3864AB")],
    ]

    default_values = {
        "-NOCNV-": True,
        "-CNV-": False,
        "-NOFFPE-": True,
        "-FFPE-": False,
        "-DEFTARG-": True,
        "-CUSTARG-": False,
    }
    textKeys = [
        "-INDIR-",
        "-OUTDIR-",
        "-PAIRS-",
        "-TARGETS-",
        "-JOBNAME-",
        "-ANNOTATION-",
    ]
    # create main layout
    logo = sg.Image(xavier_base(os.path.join("resources", "CCBRlogo.png")))
    layout = [
        [sg.Column([[logo]], justification="center")],
        [
            sg.Text(
                "XAVIER - eXome Analysis and Variant explorER",
                font=("Helvetica", 12, "bold"),
            )
        ],
        [
            sg.Text("Input Fastqs folder", size=(20, 1)),
            sg.InputText(key="-INDIR-"),
            sg.FolderBrowse(target="-INDIR-"),
        ],
        [
            sg.Text("Output folder", size=(20, 1)),
            sg.InputText(key="-OUTDIR-"),
            sg.FolderBrowse(target="-OUTDIR-"),
        ],
        [
            sg.Text("Genome", size=(20, 1)),
            sg.Combo(
                values=genome_annotation_combinations,
                key="-ANNOTATION-",
                tooltip="hg38: Homo sapiens GRCh38.p14; mm10: Mus musculus GRCm38.p6",
            ),
        ],
        [
            sg.Text("Job name", size=(20, 1)),
            sg.InputText(
                key="-JOBNAME-",
                tooltip="Name of the job for this run. All output files will be stored under this folder name in the output folder.",
            ),
        ],
        [sg.Frame("Analysis Mode", analysis_layout, visible=True)],
        [sg.Button(button_text="Additional Settings", key="-SETTINGS-")],
        [sg.Frame("", settings_layout, key="-SET-", visible=False)],
        [
            sg.Submit(key="-SUBMIT-"),
            sg.Button(button_text="Documentation", key="--DOC--"),
            sg.Button(button_text="Help", key="--HELP--"),
            sg.Cancel(key="--CANCEL--", button_color="tomato"),
        ],
    ]

    if DEBUG:
        print("layout is ready!")

    window = sg.Window(
        f"XAVIER {get_version()}", layout, location=(0, 500), finalize=True
    )
    if DEBUG:
        print("window created!")

    # Event loop:
    while True:
        event, values = window.read()
        # if DEBUG: print(event,values) ## Turn on for debugging

        # if any((event != 'Submit')):
        if event in ("--CANCEL--", sg.WINDOW_CLOSED):
            sg.popup_auto_close(
                "Thank you for running XAVIER. GoodBye!", location=(0, 500), title=""
            )
            sys.exit(69)
        if event == "-TUMNORM-":
            window["-PAIROPTS-"].update(visible=True)
            window["-ONLYOPTS-"].update(visible=False)
        elif event == "-TUMONLY-":
            window["-PAIROPTS-"].update(visible=False)
            window["-ONLYOPTS-"].update(visible=True)
            values["-CNV-"] = False
        if event == "--DOC--":
            copy_to_clipboard("https://ccbr.github.io/XAVIER/")
            sg.Popup(
                "Visit https://ccbr.github.io/XAVIER/ for links to complete documentation. The link has been copied to your clipboard. Please paste it in your favorite web browser.",
                font=("Arial", 12, "bold"),
                title="",
                location=(0, 500),
            )
            continue
        if event == "--HELP--":
            copy_to_clipboard("ccbr_pipeliner@mail.nih.gov")
            sg.Popup(
                "Email ccbr_pipeliner@mail.nih.gov for help. The email id has been copied to your clipboard. Please paste it in your emailing software.",
                font=("Arial", 12, "bold"),
                title="",
                location=(0, 500),
            )
            continue
        if event == "-SETTINGS-":
            window["-SET-"].update(visible=True)
        if event == "-DEFTARG-":
            window["-BED-"].update(visible=False)
        if event == "-CUSTARG-":
            window["-BED-"].update(visible=True)
        if event == "-DISCSET-":
            window["-SET-"].update(visible=False)
            window["-BED-"].update(visible=False)
            for key, value in default_values.items():
                window[key].Update(value)
        if event == "-SUBMIT-":
            # check for correct inputs
            if values["-INDIR-"] == "":
                sg.PopupError("Input folder must be provided!!", location=(0, 500))
                continue
            elif not os.path.exists(values["-INDIR-"]):
                sg.PopupError("Input folder doesn't exist!!", location=(0, 500))
                continue
            elif len(get_fastqs(values["-INDIR-"])) == 0:
                sg.PopupError("Input folder has no fastqs!!", location=(0, 500))
                continue
            else:
                inputfastqs = get_fastqs(values["-INDIR-"])
                if DEBUG:
                    print(inputfastqs)
                if len(inputfastqs) == 0:
                    sg.PopupError(
                        "Input folder has no fastqs!!",
                        location=(0, 500),
                        title="ERROR!",
                        font=("Arial", 12, "bold"),
                    )
                    window.Element("-INDIR-").update("")
                    continue
            if values["-OUTDIR-"] == "":
                sg.PopupError("Output folder must be provided!!", location=(0, 500))
                continue
            outputfolder = values["-OUTDIR-"] + "/" + values["-JOBNAME-"]
            if os.path.exists(outputfolder):
                ch = sg.popup_yes_no(
                    "Output folder name exists... this is probably a re-run ... proceed?",
                    title="Rerun??",
                    location=(0, 500),
                )
                if ch == "No":
                    window.Element("-OUTDIR-").update("")
                    continue
            if values["-CUSTARG-"] == True:
                if not values["-TARGETS-"]:
                    sg.PopupError(
                        "Custom Targets BED file selected but not provided!!",
                        location=(0, 500),
                    )
                    continue
            if values["-TUMNORM-"] == "" and values["-TUMONLY-"] == "":
                sg.PopupError("Select an analysis mode", location=(0, 500))
                continue
            if values["-TUMNORM-"] == True:
                if values["-PAIRS-"] == "":
                    sg.PopupError(
                        "Tumor-normal mode selected. Need Pairs file to continue",
                        location=(0, 500),
                    )
                    continue
            genome = values["-ANNOTATION-"]
            output_dir = os.path.join(values["-OUTDIR-"], values["-JOBNAME-"])
            run_args = argparse.Namespace(
                runmode="init",
                input=list(glob.glob(os.path.join(values["-INDIR-"], "*.fastq.gz"))),
                output=output_dir,
                genome=genome,
                targets=values["-TARGETS-"],
                mode="slurm",
                job_name="pl:xavier",
                callers=["mutect2", "mutect", "strelka", "vardict", "varscan"],
                pairs=values.get("-PAIRS-", None),
                ffpe=values["-FFPE-"],
                cnv=values["-CNV-"],
                wait=False,
                create_nidap_folder=False,
                silent=False,
                singularity_cache=os.environ.get("SINGULARITY_CACHEDIR", None),
                sif_cache=get_sif_cache_dir(),
                tmp_dir=get_tmp_dir(None, output_dir),
                threads=2,
            )
            allout_init = exec_in_context(run, run_args)
            run_args.runmode = "dryrun"
            allout_dryrun = exec_in_context(run, run_args)
            allout = "\n".join([allout_init, allout_dryrun])
            if DEBUG:
                print(allout)
            sg.popup_scrolled(
                allout, title="Dryrun:STDOUT/STDERR", location=(0, 500), size=(80, 30)
            )
            if "error" in allout or "Error" in allout or "ERROR" in allout:
                continue
            ch = sg.popup_yes_no(
                "Submit run to slurm?", title="Submit??", location=(0, 500)
            )
            if ch == "Yes":
                run_args.runmode = "run"
                allout = exec_in_context(run, run_args)
                sg.popup_scrolled(
                    allout,
                    title="Slurmrun:STDOUT/STDERR",
                    location=(0, 500),
                    size=(80, 30),
                )
                rerun = sg.popup_yes_no(
                    "Submit another XAVIER job?", title="", location=(0, 500)
                )
                if rerun == "Yes":
                    for key in window.read():
                        window[key].Update(value="")
                        window["-PAIROPTS-"].update(visible=False)
                        window["-ONLYOPTS-"].update(visible=False)
                        window["-TUMNORM-"].update(value=False)
                        window["-TUMONLY-"].update(value=False)
                if rerun == "No":
                    sg.popup_auto_close(
                        "Thank you for running XAVIER. GoodBye!",
                        location=(0, 500),
                        title="",
                    )
                    break
            elif ch == "No":
                for key in textKeys:
                    window[key].Update(value="")
                for key, value in default_values.items():
                    window[key].Update(value)
                window["-PAIROPTS-"].update(visible=False)
                window["-ONLYOPTS-"].update(visible=False)
                window["-TUMNORM-"].update(value=False)
                window["-TUMONLY-"].update(value=False)
                continue
    window.close()


def copy_to_clipboard(string):
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(string)
    r.update()
    r.destroy()


def fixpath(p):
    return os.path.abspath(os.path.expanduser(p))


def get_fastqs(inputdir):
    inputdir = fixpath(inputdir)
    inputfastqs = glob.glob(inputdir + os.sep + "*.fastq.gz")
    inputfqs = glob.glob(inputdir + os.sep + "*.fq.gz")
    inputfastqs.extend(inputfqs)
    return inputfastqs


def delete_files(files):
    for f in files:
        if os.path.exists(f):
            os.remove(f)


if __name__ == "__main__":
    main()
