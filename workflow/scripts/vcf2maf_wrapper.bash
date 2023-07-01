#!/usr/bin/env bash
set -e -o pipefail

# Source argparse bash wrapper for more
# advanced argument parsing in bash
source /opt2/argparse.bash \
        || exit 1

ARGPARSE_DESCRIPTION="vcf2maf wrapper script"
argparse "$@" <<EOF || exit 1
parser.add_argument('--vcf', required=True, help='input vcf file')
parser.add_argument('--maf', required=True, help='output maf file')
parser.add_argument('--genomebuild', required=True, help='hg19/hg38/mm10')
parser.add_argument('--tid', required=False, help='TumorID in VCF ... sample id of tumor sample [first column name of VCF used by default]')
parser.add_argument('--nid', required=False, help='NormalID in VCF ... sample id of normal sample ... do not provide in case of tumor-only')
parser.add_argument('--maf_tid', required=False, default='', help='TumorID for MAF... output sample id of tumor sample [--tid used by default]')
parser.add_argument('--maf_nid', required=False, default='', help='NormalID for MAF ... output sample id of normal sample [--nid used by default]')
parser.add_argument('--info', required=False, default='', help='Comma-delimited names of INFO fields to retain as extra columns in MAF [default \'\']')
parser.add_argument('--addchr', required=False, default=True, help='Whether or not to add a chr prefix to chromosome names [default=True]')
parser.add_argument('--threads', required=False, default=0, help='Number of forks to use for VEP [default=1, or $SLURM_CPUS_PER_TASK-1 if on Biowulf]')
parser.add_argument('--vepresourcebundlepath', required=False, default='/data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources',help='location of the resources required by vep')
parser.add_argument('--genomefasta', required=False, help='Fasta location')
EOF

# Set Genome aliases, vaild choices = hg19/hg38/mm10
ncbi_build=""
species=""
if [ $GENOMEBUILD == "hg38" ]; then 
    ncbi_build="GRCh38" 
    species="homo_sapiens"
elif [ $GENOMEBUILD == "hg19" ]; then 
    ncbi_build="GRCh37"
    species="homo_sapiens"
elif [ $GENOMEBUILD == "mm10" ] || [ $GENOMEBUILD == "GRCm38" ]; then
    ncbi_build="GRCm38"
    species="mus_musculus"
else
    echo "Unsupport value to option: --genomebuild"
    echo "Please select from: hg19/hg38/mm10"
    exit 1
fi

# Set paths to VEP resources based on paths 
dotvep="${VEPRESOURCEBUNDLEPATH}"

# to species and ncbi build name
#if [ $GENOME == "hg38" ]; then 
#fa="${VEPRESOURCEBUNDLEPATH}/${ncbi_build}.fa"
#elif [ $GENOME == "mm10" ]; then
#    fa="${VEPRESOURCEBUNDLEPATH}/mouse.fa"
#fi


# Detect file extension and unzip if necessary
OUTPUT_DIR=$(dirname $MAF)
VCF_BASENAME=$(basename $VCF)
vcf_file_extension=${VCF_BASENAME##*.}
if [ $vcf_file_extension == "gz" ]; then
    ungz_vcf="$OUTPUT_DIR/$(basename -s .gz $VCF)"
    echo "Decompressing VCF to $ungz_vcf..."
    gzip -dc $VCF > $ungz_vcf
    VCF=$ungz_vcf
fi
INPUT_DIR=$(dirname $VCF)

# Add chr prefix if requested and missing
chr_text=""
if [ $GENOMEBUILD == "hg38" ] || [ $GENOMEBUILD == "mm10" ] ; then 
    chr_text="chr"
fi
echo "Adding text '$chr_text' to chrom names..."
awk -v pre=$chr_text \
    '{out_text=$0; gsub("^chr","",out_text); if( /^[^#]/ ) { out_text=pre out_text}; print out_text }' \
    $VCF > $VCF.fixed

# Detect sample names in VCF
VCF_SAMPLE_IDS=($(bcftools query -l $VCF.fixed))

# By default, the first column will be
# the tumor (in case of single column)
TID_IDX=0
NID_IDX=""
VCF_NID=""
NORM_VCF_ID_ARG=""
NSAMPLES=${#VCF_SAMPLE_IDS[@]}
if [ $NSAMPLES -gt 1 ]; then
    # Assign tumor, normal IDs 
    # Look through column names and 
    # see if they match provided IDs
    for (( i = 0; i < $NSAMPLES; i++ )); do
        echo "${VCF_SAMPLE_IDS[$i]}"
        if [ "${VCF_SAMPLE_IDS[$i]}" == "$TID" ]; then
            TID_IDX=$i
        fi
        
        if [ "${VCF_SAMPLE_IDS[$i]}" == "$NID" ]; then
            NID_IDX=$i
        fi
    done

    if [ ! -z $NID_IDX ]; then
        VCF_NID=${VCF_SAMPLE_IDS[$NID_IDX]}
        NORM_VCF_ID_ARG="--vcf-normal-id $VCF_NID"
    fi
fi

VCF_TID=${VCF_SAMPLE_IDS[$TID_IDX]}
echo ""
echo "___________________"
echo $VCF_TID
echo $MAF_TID
echo $VCF_NID
echo $MAF_NID
echo "___________________"
echo ""

if [ "$MAF_TID" == "" ]; then
    MAF_TID=$VCF_TID
fi
if [ "$MAF_NID" == "" ]; then
    MAF_NID=$VCF_NID
fi

NORM_MAF_ID_ARG=""
if [ ! "$MAF_NID" == "" ]; then
    NORM_MAF_ID_ARG="--normal-id $MAF_NID"
fi

# Set option for multiple threads 
if [ $THREADS==0 ]; then
    THREADS=1
    if [ ! -z $SLURM_CPUS_PER_TASK ]; then
        # Use threads if on Biowulf
        THREADS=$((SLURM_CPUS_PER_TASK-1))
    fi
fi

echo "Running vcf2maf.pl... "
vcf2maf.pl \
    --vep-forks "$THREADS" \
    --input-vcf "${VCF}.fixed" \
    --vcf-tumor-id "$VCF_TID" $NORM_VCF_ID_ARG \
    --output-maf "${MAF}" \
    --tumor-id "$MAF_TID" $NORM_MAF_ID_ARG \
    --vep-path "/opt/vep/src/ensembl-vep" \
    --vep-data "$dotvep" \
    --ncbi-build "$ncbi_build" \
    --species "$species" \
    --retain-info "$INFO" \
    --ref-fasta "$GENOMEFASTA" \
    --vep-overwrite
    #--filter-vcf "$filtervcf" \


