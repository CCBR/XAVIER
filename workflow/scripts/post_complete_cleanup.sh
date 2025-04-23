#!/usr/bin/env bash

# Description
echo "Post success cleanup script"

# Parse arguments
for i in "$@"; do
  case $i in
    --outdir=*)
      OUTDIR="${i#*=}"
      shift
      ;;
    *)
      echo "Unknown argument: $i"
      exit 1
      ;;
  esac
done

# Check if OUTDIR is set
if [ -z "$OUTDIR" ]; then
  echo "Error: --outdir is required."
  exit 1
fi

# Echo the directory being cleaned
echo "Cleaning up XAVIER output directories in $OUTDIR..."

# Change to the output directory
cd "$OUTDIR" || { echo "Failed to cd into $OUTDIR"; exit 1; }

# Remove the directories
rm -r bams/chrom_split
rm -r somatic_paired/SNP_Indels/*/chrom_split
rm -r QC/kraken

echo "Cleanup complete!"