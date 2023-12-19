########################################################
## This script is designed to take a "covered" targets bed file for exome seq
## and make it compatible for downstream analysis
## In particular, it does two things to work with FREEC
## (1) Add a 'chr' prefix to chromosome names if genome is hg19
## (2) Remove repeated regions with the same start site (the first one is kept)

from __future__ import print_function
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_bed", help="Input BED file to be reformatted")
# parser.add_argument("-g","--genome", help="Only used if equals 'hg19'; otherwise ignored", default="hg38")
parser.add_argument("-o","--output_bed", help="Reformatted output BED file", default="exome_targets.bed")
parser.add_argument("-f","--output_fields", help="Number of fields in output BED ('3' or '6')", default="6")
args = parser.parse_args()

infile=args.input_bed
outfile=args.output_bed
nfield=str(args.output_fields)

last_start="-1" ## Position of the last start site

### Open input bed file for reading
with open(infile, 'r') as inputFile:
    ### Open output exome targets bed file for writing
    with open(outfile, 'w') as exome_bed:
        ### Step through each line of input
        for line in inputFile:
            ### Skip comments
            if not line.startswith("#") and not line.startswith("track") and not line.startswith("browser"):
                curr_cols=line.strip().split("\t")
                if (len(curr_cols) < 3):
                    sys.exit("Targets BED file must contain at least three columns: chr, start, end")
                
                
                if (len(curr_cols) < 4):
                    curr_cols[3]="."
                
                min_output=curr_cols[0] + "\t" + curr_cols[1]  + "\t" + curr_cols[2]
                
                extra_fields=""
                if (not(nfield == "3")):
                    extra_fields="\t" + curr_cols[3] + "\t0\t.\n"
                    
                bed_output=min_output + extra_fields
                
                ### Add 'chr' prefix if genome is hg19
                # if (args.genome=="hg19"):
                #     freec_bed_output="chr" + curr_cols[0].lstrip("chr") + "\t" + curr_cols[1]  + "\t" + curr_cols[2] + "\n"
                # else:
                #     freec_bed_output=curr_cols[0] + "\t" + curr_cols[1]  + "\t" + curr_cols[2] + "\n"
                
                ### If current start location is same as previous, output empty string
                if (curr_cols[1] == last_start):
                    print(curr_cols[1] + " IS equal to " + last_start + " so skipping it...")
                    bed_output=""

                ### Write to both files
                exome_bed.write(bed_output)
                # freec_bed.write(freec_bed_output)
                
                ### Update loop variables
                last_start=curr_cols[1]



