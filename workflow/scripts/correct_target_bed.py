#!/usr/bin/env python3

__desc__ = """
# FREEC errors out with
##################################################################################
# Error: your BED file with coordinates of targeted regions may contain duplicates
##################################################################################
# This script tries to fix the issue by attempting to fix the target BED file
# parsed to FREEC
# reference: https://github.com/BoevaLab/FREEC/issues/43
# @parameters:
# @param1: input BED file ... file which causes FREEC to error out
# @param2: output BED file .. file which may fix FREEC
"""

import sys
import subprocess
import os

if len(sys.argv) != 3:  # 2 arguments are required
    print(__desc__)
    print("Usage: python3 " + sys.argv[0] + " <param1> <param2>")
    exit()

inputBed = sys.argv[1]
outputBed = sys.argv[2]

iBed = open(inputBed, "r")
oBed = open(outputBed + ".tmp", "w")

# collapse redundant regions
# regiongs with identical chrom/start/end but with possibly different annotations
# i.e. differences in column 4 or 5 are collapsed

annotations = dict()
for line in iBed.readlines():
    line = line.strip().split("\t")
    region = line[0] + "##" + line[1] + "##" + line[2]
    if not region in annotations:
        annotations[region] = dict()
        annotations[region]["1"] = list()
        annotations[region]["2"] = list()
    a1 = line[3].split(",")
    a2 = line[4].split(",")
    annotations[region]["1"].extend(a1)
    annotations[region]["2"].extend(a2)
iBed.close()

for k, v in annotations.items():
    region = k.split("##")
    oBed.write(
        "%s\t%s\t%s\t%s\t%s\t.\n"
        % (region[0], region[1], region[2], ",".join(v["1"]), ",".join(v["2"]))
    )
oBed.close()

# sort the collapsed file

sortcmd = "sort -k1,1 -k2,2n -k3,3n " + outputBed + ".tmp" + " > " + outputBed + ".tmp2"
subprocess.run(sortcmd, shell=True, check=True)

# Based off of the suggestions here --> https://github.com/BoevaLab/FREEC/issues/43
# FREEC treats repeated chrom+start as duplicates even though they may have
# different ends and annotations
# Solution: increment start by 1 to create "non-duplicate" entry for FREEC

incrementcmd = (
    'awk -F"\\t" -v OFS="\\t" \'{seen[$1"##"$2]+=1;if (seen[$1"##"$2]!=1){$2=$2+seen[$1"##"$2]-1};print}\' '
    + outputBed
    + ".tmp2"
    + ">"
    + outputBed
    + ".tmp3"
)
subprocess.run(incrementcmd, shell=True, check=True)

# resort the corrected file

sortcmd = "sort -k1,1 -k2,2n -k3,3n " + outputBed + ".tmp3" + " > " + outputBed
subprocess.run(sortcmd, shell=True, check=True)

# delete intermediate files

os.remove(outputBed + ".tmp")
os.remove(outputBed + ".tmp2")
os.remove(outputBed + ".tmp3")
