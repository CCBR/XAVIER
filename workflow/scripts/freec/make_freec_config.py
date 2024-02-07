#! /usr/local/bin/python3

import argparse
import os

parser = argparse.ArgumentParser(
    description="Make a config file for running CONTROL-Freec."
)
parser.add_argument("-t", "--tumor", help="Tumor BAM file")
parser.add_argument("-n", "--normal", help="Normal BAM file")
parser.add_argument("-g", "--genome-fasta", help="Full genome fasta")
parser.add_argument(
    "-l", "--chrom-lengths", help="File with chromosome lengths (an .fai file)"
)
parser.add_argument("-s", "--chrom-seqs", help="Folder with chromosome fastqs")
parser.add_argument(
    "-r", "--capture-regions", help="BED file with exome capture regions", default=None
)
parser.add_argument(
    "-o",
    "--output-config",
    help="Filename of the config file to be generated",
    default="freec_config.txt",
)
parser.add_argument(
    "--snps-file", help="Passed to 'minCNAlength' argument", default="3"
)
parser.add_argument(
    "--degree",
    help="Polynomial degree used for GC normalization used for 'degree' argument",
    default="3",
)
parser.add_argument(
    "--min-length", help="Passed to 'minCNAlength' argument", default="3"
)
parser.add_argument(
    "--min-read-count", help="Passed to 'readCountThreshold' argument", default="50"
)
parser.add_argument("--pileup-vcf", help="Passed to 'makePileup' argument", default="")
parser.add_argument("--ploidy", help="Estimated ploidy", default="")
parser.add_argument("--contamination", help="Estimated contamination", default="")
parser.add_argument("--threads", help="Number of threads to use", default="2")
args = parser.parse_args()

if not os.path.exists(args.chrom_lengths):
    print("Can't find chr lengths file: " + args.chrom_lengths)
    exit

if not os.path.exists(os.path.dirname(args.output_config)):
    print("Can't output directory: " + os.path.dirname(args.output_config))
    exit

config_contents = list()
config_contents.append("[general]")
config_contents.append("")
config_contents.append("degree = " + args.degree)
config_contents.append("minCNAlength = " + args.min_length)
config_contents.append("readCountThreshold = " + args.min_read_count)
config_contents.append("chrLenFile = " + args.chrom_lengths)
config_contents.append("ploidy = " + args.ploidy)
config_contents.append("contamination = " + args.contamination)
config_contents.append("BedGraphOutput = TRUE")
config_contents.append("forceGCcontentNormalization = 1")
config_contents.append("noisyData = TRUE")
config_contents.append("breakPointThreshold = 0.8")
config_contents.append("window = 0")
config_contents.append("chrFiles = " + args.chrom_seqs)
config_contents.append("minimalSubclonePresence = 30")
config_contents.append("printNA = FALSE")
config_contents.append("contaminationAdjustment = TRUE")
config_contents.append("maxThreads = " + args.threads)
config_contents.append("numberOfProcesses = " + args.threads)
config_contents.append("outputDir = " + os.path.dirname(args.output_config))
config_contents.append("")

config_contents.append("[sample]")
config_contents.append("")
config_contents.append("mateFile = " + args.tumor)
config_contents.append("inputFormat = BAM")
config_contents.append("mateOrientation = FR")
config_contents.append("")

config_contents.append("[control]")
config_contents.append("")
config_contents.append("mateFile = " + args.normal)
config_contents.append("inputFormat = BAM")
config_contents.append("mateOrientation = FR")
config_contents.append("")

config_contents.append("[BAF]")
config_contents.append("")
config_contents.append("makePileup = " + args.pileup_vcf)
config_contents.append("fastaFile = " + args.genome_fasta)
config_contents.append("minimalCoveragePerPosition = 20")
config_contents.append("nminimalQualityPerPosition = 20")
config_contents.append("SNPfile =" + args.snps_file)
config_contents.append("")

if args.capture_regions is not None:
    if os.path.exists(args.capture_regions):
        config_contents.append("[target]")
        config_contents.append("")
        config_contents.append("captureRegions = " + args.capture_regions)
        config_contents.append("")


with open(args.output_config, "w") as config_out:
    config_out.write("\n".join(config_contents))


########################################################################################################
#### Justin's original perl script: /data/CCBR_Pipeliner/4.0.2/Pipeliner/Results-template/Scripts/make_freec_pass2_exome_tn_config.pl
# #!/usr/bin/perl -w
# use strict;
# use List::Util 'shuffle';
#
# #INPUT
#
# #my $mergedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_merged.maf'; #to fix...
# #open C, ">$mergedmaf";
#
# my $outfile = $ARGV[0] . '/freec_exome_config.txt';
# my $chrLenFile = $ARGV[1];
# my $chrFiles = $ARGV[2];
# my $tumormateFile = $ARGV[3];
# my $controlmateFile = $ARGV[4];
# my $makePileup = $ARGV[5];
# my $fastaFile = $ARGV[6];
# my $SNPfile = $ARGV[7];
# my $targets = $ARGV[8];
# my @line=();
# my $contamination='';
# my $ploidy='';
# my $rep=0;
#
# my $infile=$ARGV[9];
# open G, "<$infile";
# while (<G>){
# 	chomp;
#   	last if m/^$/;
#   	@line = split;
#   	next if ($line[0] =~ m'cellularity');
#   	if ($rep==0) {
# 		$contamination=(1-$line[0]);
# 		$ploidy=$line[1];
# 		$rep++;
# 	}
# }
#
# open C, ">$outfile";
#
# print C '[general]' . "\n\n";
#
# print C "BedGraphOutput = TRUE\ndegree = 1\nforceGCcontentNormalization = 1\nminCNAlength = 3\nnoisyData = TRUE\nreadCountThreshold = 50\n";
# print C "chrLenFile = $chrLenFile\n";
# print C "ploidy = $ploidy\ncontamination=$contamination\nbreakPointThreshold = 0.8\nwindow = 0\n";
# print C "chrFiles = $chrFiles\n";
# print C "minimalSubclonePresence = 30\nprintNA = FALSE\ncontaminationAdjustment = TRUE\nmaxThreads = 24\nnumberOfProcesses = 24\n";
# print C "outputDir = $ARGV[0]\n\n";
#
# print C '[sample]' . "\n\n";
#
# print C "mateFile = $tumormateFile\n";
# print C "inputFormat = BAM\nmateOrientation = FR\n\n";
#
# print C '[control]' . "\n\n";
#
# print C "mateFile = $controlmateFile\n";
# print C "inputFormat = BAM\nmateOrientation = FR\n\n";
#
# print C '[target]' . "\n\n";
#
# print C "captureRegions = $targets\n\n";
#
# print C '[BAF]' . "\n\n";
#
# print C "makePileup = $makePileup\n";
# print C "fastaFile = $fastaFile\n";
# print C "minimalCoveragePerPosition = 20\nminimalQualityPerPosition = 20\n";
# print C "SNPfile = $SNPfile";
