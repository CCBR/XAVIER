# Output files

## XAVIER

The output files and their locations are broken down here for the XAVIER pipeline. Pre-processing and germline variant calling steps are common but somatic variant calling is dependent on whether the pipeline was run in either (A) tumor-normal pair or (B) tumor-only analysis mode. All file locations are relative to the output directory specified during the job submission.

The output directory after a complete XAVIER run should look like:

```bash
xavier_output/
├── bams
├── cluster.json # cluster info for the run
├── config
├── config.json  # config file for the run
├── fastqs
├── germline
├── indels.vcf.gz[.tbi] # raw germline INDELs
├── input_files
├── intervals.list
├── {sample1}-normal.R1.fastq.gz -> /path/to/{sample1}-normal.R1.fastq.gz
├── {sample1}-normal.R2.fastq.gz -> /path/to/{sample1}-normal.R2.fastq.gz
├── {sample1}-tumor.R1.fastq.gz -> /path/to/{sample1}-tumor.R1.fastq.gz
├── {sample1}-tumor.R2.fastq.gz -> /path/to/{sample1}-tumor.R2.fastq.gz
.
.
.
├── kickoff.sh
├── logfiles
├── QC
├── resources
├── snps.vcf.gz[.tbi] # raw germline SNPs
├── somatic_paired # in case of tumor-normal paired run
├── somatic_tumor_only # in case of tumor-only run
└── workflow

```

Below we describe the different folders that contain specific outputs obtained for all samples from the XAVIER pipeline

### 1. `QC`

The `QC` folder contains all the Quality-Control analyses performed at different steps of the pipeline for each sample to assess sequencing quality before and after adapter trimming, microbial taxonomic composition, contamination, variant calling, etc. The final summary report and data is available `finalQC` folder. \
The MultiQC report also contains results from other analysis like mapping statistics, ancestry and relatedness, etc. It is recommended to study the MultiQC report first to get a birds eye view of the sequence data quality.

```bash
QC/
├── exome_targets.bed
├── finalQC/
│   ├── MultiQC_Report_data
│   └── MultiQC_Report.html
├── FQscreen
├── sample1-normal
├── sample1-normal_fastqc.html
.
.
.
├── sample1-tumor
├── sample1-tumor_fastqc.html
.
.
.
├── kraken
├── raw_variants.het
├── raw_variants.variant_calling_detail_metrics
└── raw_variants.variant_calling_summary_metrics
```

### 2. `bams`

The `bams` folder contain two subfolders `chrom_split` and `final_bams`. `final_bams` contains the final processed BAM files for each sample in the run and the `chrom_split` folder contains all the sample BAM files split by each chromosome.

```bash
bams/
├── chrom_split
│   ├── {sample1}-normal.chr1.split.bam[.bai]
│   ├── {sample1}-normal.chr2.split.bam[.bai]
.   .
.   .
.   .
│   └── {sampleN}-tumor.chrN.split.bam[.bai]
└── final_bams
    ├── {sample1}-normal.bam[.bai]
    ├── {sample1}-tumor.bam[.bai]
    .
    .
    .
    └── {sampleN}-tumor.bam.bai
```

### 3. `germline`

This folder contains the output from the GATK Best Practices pipeline to obtain germline variants with a few alterations detailed below. Briefly, joint SNP and INDEL variant detection is conducted across all samples included in a pipeline run using the GATK Haplotypcaller under default settings. Raw variants are then subsequently filtered based on several GATK annotations: \
A strict set of criteria (QD < 2.0, FS > 60.0, MQ < 40.0, MQRankSum < -12.5, ReadPosRankSum < -8.0 for SNPs; QD < 2.0, FS > 200.0, ReadPosRankSum < -20.0 for INDELs) generates the 'combined.strictFilter.vcf'.

This call set is highly stringent, maximizing the true positive rate at the expense of an elevated false negative rate. This call set is really only intended for more general population genetic scale analyses (e.g., burden tests, admixture, linkage/pedigree based analysis, etc.) where false positives can be significantly confounding.

In case of human sequence data, a basic analyses of sample relatedness and ancestry (e.g., % European, African, etc.) is also performed using somalier.

The output folder looks like:

```bash
germline/
├── gVCFs
.
.
.
├── somalier # only for hg38 genome
└── VCF
```

The `VCF` folder contains the final filtered germline variants (SNPs and INDELs) for all samples combined. The folder also contains raw variants for each sample, all samples combined, and also combined raw variants split by chromosome.

```bash
VCF/
├── by_chrom
│   ├── raw_variants_byChrom.list
│   ├── raw_variants.chr1.vcf.gz[.tbi]
.   .
.   .
.   .
│   └── raw_variants.chrN.vcf.gz[.tbi]
├── indel.filterd.vcf.gz[.tbi]
├── {sample1}-normal.germline.vcf.gz[.tbi]
├── {sample1}-tumor.germline.vcf.gz[.tbi]
.
.
.
├── {sampleN}-normal.germline.vcf.gz[.tbi]
├── {sampleN}-tumor.germline.vcf.gz[.tbi]
├── raw_variants.vcf.gz[.tbi]
├── snp.filtered.vcf.gz[.tbi]
└── snp_indel.filtered.vcf.gz[.tbi]
```

### 4. `logfiles`

This folder contains the snakemake log files and computational statistics for the XAVIER run. All the log files (i.e., standard output and error) for each individual step are in the `slurmfiles` folder. These logfiles are important to diagnose errors in case the pipeline fails.

```bash
logfiles/
├── master.log
├── mjobid.log
├── runtime_statistics.json
├── slurmfiles
├── snakemake.log
├── snakemake.log.jobby
└── snakemake.log.jobby.short
```

## Tumor-normal pair

### `somatic_paired`

This workflow calls somatic SNPs and INDELs using multiple variant detection algorithms. For each of these tools, variants are called in a paired tumor-normal fashion, with default settings. See **Pipeline Details** for more information about the tools used and their parameter settings.

For each sample, the resulting VCF is fully annotated using VEP and converted to a MAF file using the vcf2maf tool. Resulting MAF files are found in `maf` folder within each caller's results directory (i.e., `mutect2_out`, `strelka_out`, etc.). Individual sample MAF files are then merged and saved in `merged_somatic_variants` directory.

For Mutect2, we use a panel of normals (PON) developed from the ExAC (excluding TCGA) dataset, filtered for variants <0.001 in the general population, and also including and in-house set of blacklisted recurrent germline variants that are not found in any population databases.

For Copy Number Variants (CNVs), two tools are employed in tandem. First, Control-FREEC is run with default parameters. This generates pileup files that can be used by Sequenza, primarily for jointly estimating contamination and ploidy. These value are used to run Freec a second time for improved performance.

The output directory should look like:

```bash
somatic_paired/
├── CNV # only if CNVs analyzed
│   ├── freec_out
│   └── sequenza_out
├── ffpe_filter # only if FFPE filter applied
├── qc
└── SNP_Indels
    ├── merged_somatic_variants
    │   ├── cohort_summary
    │   ├── maf # Final merged MAFs for each sample
    │   └── vcf
    ├── mutect2_out
    │   ├── chrom_split
    │   ├── cohort_summary
    │   ├── maf
    │   ├── pileup_summaries
    │   ├── read_orientation_data
    │   └── vcf
    ├── mutect_out
    │   ├── maf
    │   ├── .
    │   ├── .
    │   └── vcf
    ├── strelka_out
    ├── vardict_out
    └── varscan_out
```

## Tumor-only

### `somatic_tumor_only`

In general, the tumor-only pipeline is a stripped down version of the tumor-normal pipeline. We only run MuTect2, Mutect, and VarDict for somatic variant detection, with the same PON and filtering as described above for the tumor-normal pipeline.

```bash
somatic_tumor_only/
├── ffpe_filter # only if FFPE filter applied
├── qc
└── SNP_Indels
    ├── merged_somatic_variants
    │   ├── cohort_summary
    │   ├── maf # Final merged MAFs for each sample
    │   └── vcf
    ├── mutect2_out
    │   ├── chrom_split
    │   ├── cohort_summary
    │   ├── maf
    │   ├── read_orientation_data
    │   └── vcf
    ├── mutect_out
    │   ├── maf
    │   ├── .
    │   ├── .
    │   └── vcf
    ├── vardict_out
    └── varscan_out
```
