# Quality-control Non-Human related rules
localrules: somalier_extract
rule somalier_extract:
    """
    To estimate ancestry, Somalier first extracts known sites from mapped reads
    @Input:
        Mapped and pre-processed BAM file
    @Output:
        Exracted sites in (binary) somalier format
    """
    input:
        bam = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
        bai = os.path.join(output_bamdir,"final_bams","{samples}.bai"),
    output: 
        somalierOut = os.path.join(output_germline_base,"somalier","{samples}.somalier")
    params:
        sites_vcf = config['references']['SOMALIER']['SITES_VCF'],
        genomeFasta = config['references']['GENOME'],
        rname = 'somalier_extract'
    container: config['images']['wes_base']
    shell: """ 
    echo "Extracting sites to estimate ancestry"
    somalier extract \\
        -d "$(dirname {output.somalierOut})" \\
        --sites {params.sites_vcf} \\
        -f {params.genomeFasta} \\
        {input.bam}
    """

localrules: somalier_analysis
rule somalier_analysis:
    """
    To estimate relatedness, Somalier uses extracted site information to
    compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Exracted sites in (binary) somalier format for ALL samples in the cohort
    @Output:
        Separate tab-separated value (TSV) files with relatedness and ancestry outputs
    """
    input:
        somalier=expand(os.path.join(output_germline_base,"somalier","{samples}.somalier"), samples=samples),
    output:
        relatedness = os.path.join(output_germline_base,"somalier","relatedness.pairs.tsv"),
        relatednessSamples = os.path.join(output_germline_base,"somalier","relatedness.samples.tsv"),
        ancestry = os.path.join(output_germline_base,"somalier","ancestry.somalier-ancestry.tsv"),
        finalFileGender = os.path.join(output_germline_base,"predicted.genders.tsv"),
        finalFilePairs = os.path.join(output_germline_base,"predicted.pairs.tsv"),
        ancestoryPlot = os.path.join(output_germline_base,"sampleAncestryPCAPlot.html"),
        pairAncestoryHist = os.path.join(output_germline_base,"predictedPairsAncestry.pdf"),
    params:
        sites_vcf = config['references']['SOMALIER']['SITES_VCF'],
        genomeFasta = config['references']['GENOME'],
        script_path_gender = config['scripts']['genderPrediction'],
        script_path_samples = config['scripts']['combineSamples'],
        script_path_pca = config['scripts']['ancestry'],
        rname = 'somalier_analysis'
    container: config['images']['wes_base']
    shell: """ 
    echo "Estimating relatedness"
    somalier relate \\
        -o "$(dirname {output.relatedness})/relatedness" \\
        {input.somalier}

    Rscript {params.script_path_gender} \\
        {output.relatednessSamples} \\
        {output.finalFileGender}    
    
    Rscript {params.script_path_samples} \\
        {output.relatedness} \\
        {output.finalFilePairs}
    
    """



rule multiqc:
    """
    Reporting step to aggregate sample summary statistics and quality-control
    information across all samples. This will be one of the last steps of the 
    pipeline. The inputs listed here are to ensure that this step runs last. 
    During runtime, MultiQC will recurively crawl through the working directory
    and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    input:  
        expand(os.path.join(output_fqdir,"{samples}.fastq.info.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"FQscreen","{samples}.R2.trimmed_screen.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.krona.html"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}_fastqc.zip"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}","genome_results.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.samtools_flagstat.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.bcftools_stats.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.eval.grp"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.html"), samples=samples),
        os.path.join(output_qcdir,"raw_variants.het"), 
        os.path.join(output_qcdir,"raw_variants.variant_calling_detail_metrics"),
    output: 
        report  = os.path.join(output_qcdir,"finalQC","MultiQC_Report.html"),
    params: 
        rname  = "multiqc",
        workdir = os.path.join(BASEDIR)
    envmodules: 'multiqc/1.11'
    container: config['images']['multiqc']
    shell: """
    multiqc --ignore '*/.singularity/*' \\
        --ignore '*/*/*/*/*/*/*/*/pyflow.data/*' \\
        --ignore 'slurmfiles/' \\
        -f --interactive \\
        -n {output.report} \\
        {params.workdir}
    """
