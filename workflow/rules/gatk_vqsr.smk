rule gatk_vqsr:
    """
    ##Hold off on implementing due to mouse genomes not having all the same resources for VQSR and the sample sizes required.
    ##Implement Deepvariant instead
    Run GATK VQSR on the SNP and INDEls
    @Input:
        Multi-sample gVCF with all chromosomes combined
    @Output:
        Variants scored by VQSLOD
    """
    input:
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output:
       indelvcf = os.path.join(output_germline_base,"VCF","indel.recalibrated.vcf.gz"),
       snpindelvcf = os.path.join(output_germline_base,"VCF","snp_indel.recalibrated.vcf.gz")
    params:
        genome=config['references']['GENOME'],
        mills=config['references']['MILLS'],
        axiom=config['references']['AXIOM'],
        dbsnp=config['references']['DBSNP'],
        hapmap=config['references']['HAPMAP'],
        omni=config['references']['OMNI'],
        onekgp=config['references']['1000GSNP'],
        rname="vqsr",
        ver_gatk=config['tools']['gatk4']['version']
    message: "Running GATK4 VQSR on Cohort VCF input file"
    envmodules: config['tools']['gatk4']['modname']
    container: config['images']['wes_base']
    shell:
        """
        gatk --java-options VariantRecalibrator \\
        -R {params.genome}
        -V {input.vcf} \\
        --trust-all-polymorphic \\
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \\
        -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \\
        -mode INDEL \\
        --max-gaussians 4 \\
        -resource:mills,known=false,training=true,truth=true,prior=12 {params.mills}\\
        -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {params.axiom} \\
        -resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp} \\
        --tranches-file cohort_indels.tranches \\
        -O cohort_indels.recal

        gatk --java-options VariantRecalibrator \\
        -R {params.genome} \\
        -V {input.vcf} \\
        --trust-all-polymorphic \\
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \\
        --resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.onekgp} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \\
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode SNP \\
        --max-gaussians 6 \\
        -O cohort_snps.recal \\
        --tranches-file output_snp.tranches \\
        --rscript-file output.plots.SNP.R

        gatk --java-options ApplyVQSR \\
        -V {input.vcf} \\
        --recal-file  cohort_indels.recal \\
        --tranches-file cohort_indels.tranches \\
        --truth-sensitivity-filter-level 99.7 \\
        --create-output-variant-index true \\
        -mode INDEL \\
        -O {output.indelvcf}

        gatk --java-options ApplyVQSR \\
        -V indel.recalibrated.vcf.gz \\
        --recal-file cohort_snps.recal \\
        --tranches-file output_snp.tranches \\
        --truth-sensitivity-filter-level 99.7 \\
        --create-output-variant-index true \\
        -mode SNP \\
        -O {output.snpindelvcf}

        """
