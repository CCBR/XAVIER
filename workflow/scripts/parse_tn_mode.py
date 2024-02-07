#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function, division
import sys, gzip


# USAGE
# sys.argv[1] = sample_name.R1.fastq.gz
# sys.argv[2] = sample_name (name without PATH and .R?.fastq.gz extension)
# Example
# $ python get_flowcell_lanes.py input.R1.fastq.gz input > flowcell_lanes.txt


def usage(message="", exitcode=0):
    """Displays help and usage information. If provided invalid usage
    returns non-zero exit-code. Additional message can be displayed with
    the 'message' parameter.
    """
    print(
        "Usage: python {} paired  [/path/to/pairs] [ sample_names ] > sampleName.flowcell_lanes.txt".format(
            sys.argv[0]
        )
    )
    if message:
        print(message)
    sys.exit(exitcode)


def read_pairsfile(tn_mode="auto", pairs_filepath="", sample_names=[]):
    ## Make sure tn_mode is valid
    if not tn_mode in ["auto", "paired", "tumor_only"]:
        raise NameError(
            """\n\tFatal: tn_mode must be one of 'auto', 'paired', or 'tumor_only'
        Argument received: {}
        """.format(
                tn_mode, sys.argv[0]
            )
        )

    ## Initialize some empty variables
    tumor_ids = []
    normal_ids = []
    paired_ids = {}

    ## If pairs file exists, try to use it
    if os.path.isfile(pairs_filepath):
        ## Read pairs file as data frame
        df = pd.read_csv(pairs_filepath, header=0, sep="\t")
        df.columns = df.columns.str.lower()  ## Make column names case-insensitive

        ## Make sure it contains a "tumor" column
        if not "tumor" in df:
            raise NameError(
                """\n\tFatal: Pairs file must contain at least a 'tumor' column
            Columns found: {}
            """.format(
                    df.columns.tolist(), sys.argv[0]
                )
            )

        df = df[pd.notna(df["tumor"])]  ## Remove rows where tumor id is empty/na
        tumor_ids = df["tumor"]

        if "normal" in df:
            normal_ids = df["normal"]

        ## Make sure normal ids are not empty/na
        if any(pd.notna(normal_ids)):
            t_pair = tumor_ids[pd.notna(normal_ids)]
            n_pair = normal_ids[pd.notna(normal_ids)]
            paired_ids = dict(zip(t_pair.tolist(), n_pair.tolist()))

    ## If pairs file not found, try to use provided sample names as tumor-only IDs
    else:
        if tn_mode == "paired":
            print("WARNING: Paired mode selected without a valid pairs file!!!")

        if not sample_names:
            raise NameError(
                """\n\tFatal: Either a valid pairs file or sample names must be provided.
            Pairs file path provided: {}
            Sample names provided: {}
            """.format(
                    pairs_filepath, sample_names, sys.argv[0]
                )
            )
        else:
            tumor_ids = sample_names

    ## Overlap with given sample names
    if sample_names:
        overlapped_pairs = {k: paired_ids[k] for k in sample_names if k in paired_ids}
        overlapped_tumors = list(set(tumor_ids) & set(sample_names))

        print(
            str(len(overlapped_pairs))
            + " of "
            + str(len(paired_ids))
            + " pairs in pairs file matched given sample names"
        )
        print(
            str(len(overlapped_tumors))
            + " of "
            + str(len(tumor_ids))
            + " tumors in pairs file matched given sample names"
        )

        paired_ids = overlapped_pairs
        tumor_ids = overlapped_tumors

    out_dict = {"pairs": paired_ids, "tumors": set(tumor_ids)}

    if tn_mode == "paired":
        out_dict["tumors"] = []
    elif tn_mode == "tumor_only":
        out_dict["pairs"] = []

    return out_dict


def md5sum(filename, blocksize=65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks defined by the blocksize parameter. This is a safer
    option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return hasher.hexdigest() <str>:
        MD5 checksum of the file's contents
    """
    import hashlib

    hasher = hashlib.md5()
    with open(filename, "rb") as fh:
        buf = fh.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = fh.read(blocksize)

    return hasher.hexdigest()


if __name__ == "__main__":
    # Check Usage
    if "-h" in sys.argv or "--help" in sys.argv or "-help" in sys.argv:
        usage(exitcode=0)
    elif len(sys.argv) != 3:
        usage(
            message="Error: failed to provide all required positional arguments!",
            exitcode=1,
        )

    # Get file name and sample name prefix
    filename = sys.argv[1]
    sample = sys.argv[2]
    # Get md5 checksum
    md5 = md5sum(filename)

    # Get Flowcell and Lane information
    handle = reader(filename)
    meta = {"flowcell": [], "lane": [], "flowcell_lane": []}
    i = 0  # keeps track of line number
    with handle(filename, "r") as file:
        print(
            "sample_name\ttotal_read_pairs\tflowcell_ids\tlanes\tflowcell_lanes\tmd5_checksum"
        )
        for line in file:
            line = line.strip()
            if i % 4 == 0:  # read id or sequence identifier
                fc, lane = get_flowcell_lane(line)
                fc = fc.lstrip("@")
                fc_lane = "{}_{}".format(fc, lane)
                if fc not in meta["flowcell"]:
                    meta["flowcell"].append(fc)
                if lane not in meta["lane"]:
                    meta["lane"].append(lane)
                if fc_lane not in meta["flowcell_lane"]:
                    meta["flowcell_lane"].append(fc_lane)
            i += 1

    print(
        "{}\t{}\t{}\t{}\t{}\t{}".format(
            sample,
            int(i / 4),
            ",".join(sorted(meta["flowcell"])),
            ",".join(sorted(meta["lane"])),
            ",".join(sorted(meta["flowcell_lane"])),
            md5,
        )
    )
