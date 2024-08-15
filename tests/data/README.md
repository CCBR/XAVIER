# About

These input files are used for continuous integration purposes, specifically to dry run the pipeline whenever commits have been made to the main, master, or unified branches.

Human whole exome sequence reads from the Sequencing Quality Control Phase 2 (SEQC2) Consortium has been subsampled and added. 

The tumor-normal paired reads were downloaded from the [seqc2](https://sites.google.com/view/seqc2/home/sequencing) server that were sequenced by the NCI (WES_NC_T_1 vs. WES_NC_N_1) which corresponds to NCBI SRA accession no. [SRX4728524](https://www.ncbi.nlm.nih.gov/sra/SRX4728524) and [SRX4728523](https://www.ncbi.nlm.nih.gov/sra/SRX4728523) respectively. 

Next, the reads were subsampled to 0.1% using `seqtk` and gzipped as follows:

```bash
seqtk sample -s100 {input}.R[1/2].fastq.gz 0.001 > {input}.R[1/2]_sub.R2.fastq
gzip *.fastq
```

Similarly, the BAM files were created by first mapping to the hg38 genome and then subsampled using `samtools`:

```bash
samtools view -s 0.00125 -b WES_NC_[T/N]_1.bam -o WES_NC_[T/N]_1_sub.bam
```