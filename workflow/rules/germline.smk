# Rules for germline variant calling
rule haplotypecaller:
    """
    Germline variant calling. This can be done independently across the
    genome, so we're splitting it up by chromosome.
    @Input:
        Aligned reads in BAM format
    @Output:
        Single-sample gVCF
    """
    input: 
        bam = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
        bai = os.path.join(output_bamdir,"final_bams","{samples}.bai"),
    output:
        gzvcf = temp(os.path.join(output_germline_base,"gVCFs","{samples}.{chroms}.g.vcf.gz")),
        index = temp(os.path.join(output_germline_base,"gVCFs","{samples}.{chroms}.g.vcf.gz.tbi")),
    params: 
        sample = "{samples}",
        genome = config['references']['GENOME'],
        snpsites=config['references']['DBSNP'],
        chrom="{chroms}",
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "hapcaller"
    message: "Running GATK4 HaplotypeCaller on '{input.bam}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        myoutdir="$(dirname {output.gzvcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
         
        gatk --java-options '-Xmx24g' HaplotypeCaller \\
            --reference {params.genome} \\
            --input {input.bam} \\
            --use-jdk-inflater \\
            --use-jdk-deflater \\
            --emit-ref-confidence GVCF \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --dbsnp {params.snpsites} \\
            --output {output.gzvcf} \\
            --intervals {params.chrom} \\
            --max-alternate-alleles 3
        """


rule mergegvcfs:
    """
    Merge germline variants across samples
    @Input:
        Single-sample gVCFs, scattered across chromosomes
    @Output:
        Multi-sample gVCF, scattered across chromosomes
    """
    input: gzvcf = expand(os.path.join(output_germline_base,"gVCFs","{samples}.{{chroms}}.g.vcf.gz"),samples=samples),
           index = expand(os.path.join(output_germline_base,"gVCFs","{samples}.{{chroms}}.g.vcf.gz.tbi"),samples=samples),
           # list = "gVCFs/gVCFs.{chroms}.list",
    output:
        gzvcf = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz"),
        index = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz.tbi"),
    params: 
        genome = config['references']['GENOME'],
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "mergegvcfs"
    message: "Running GATK4 CombineGVCFs on '{input.gzvcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        input_str="--variant $(echo "{input.gzvcf}" | sed -e 's/ / --variant /g')"
        
        gatk --java-options '-Xmx24g' CombineGVCFs \\
            --reference {params.genome} \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            $input_str \\
            --output {output.gzvcf} \\
            --intervals {wildcards.chroms} \\
            --use-jdk-inflater \\
            --use-jdk-deflater
        """


rule genotype:
    """
    Joint genotyping of germline variants
    @Input:
        Multi-sample gVCF, scattered across chromosomes
    @Output:
        Multi-sample gVCF, scattered across chromosomes (with joint genotyping updates)
    """
    input: 
        gzvcf = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz"),
        index = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz.tbi"),
    output:
        vcf = os.path.join(output_germline_base,"VCF","by_chrom","raw_variants.{chroms}.vcf.gz"),
    params:
        genome = config['references']['GENOME'],
        snpsites=config['references']['DBSNP'],
        chr="{chroms}",
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "genotype"
    message: "Running GATK4 GenotypeGVCFs on '{input.gzvcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        myoutdir="$(dirname {output.vcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
        
        gatk --java-options '-Xmx96g' GenotypeGVCFs \\
            --reference {params.genome} \\
            --use-jdk-inflater \\
            --use-jdk-deflater \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --dbsnp {params.snpsites} \\
            --output {output.vcf} \\
            --variant {input.gzvcf} \\
            --intervals {params.chr}
        """


rule germline_merge_chrom:
    """
    Combine joint genotyping from all chromosomes
    @Input:
        Multi-sample gVCF for all chromosomes
    @Output:
        Multi-sample gVCF with all chromosomes combined
    """
    input:
        expand(os.path.join(output_germline_base,"VCF","by_chrom","raw_variants.{chroms}.vcf.gz"), chroms=chroms),
    output:
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
        clist = os.path.join(output_germline_base,"VCF","by_chrom","raw_variants_byChrom.list"),
    params:
        rname = "merge_chrom", genome = config['references']['GENOME']
    message: "Running GATK4 MergeVcfs on all chrom split VCF files"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        # Avoids ARG_MAX issue which limits max length of a command
        ls --color=never -d $(dirname "{output.clist}")/raw_variants.*.vcf.gz > "{output.clist}"

        gatk MergeVcfs \\
            -R {params.genome} \\
            --INPUT {output.clist} \\
            --OUTPUT {output.vcf}
        """


rule gatk_vqsr: 
    """
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
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        #gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \\
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

#        gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \\        
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

rule Gatk_SelectVariants:
    """
    Make individual VCFs with variant sites in that sample
    @Input:
        Multi-sample gVCF with all chromosomes combined
    @Output:
        Single-sample VCF with unfiltered germline variants
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","snp_indel.recalibrated.vcf.gz"),
    output: 
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz")
    params: 
        genome=config['references']['GENOME'], 
        Sname = "{samples}", 
        rname="varselect",
        ver_gatk=config['tools']['gatk4']['version'],
        targets=exome_targets_bed
    message: "Running GATK4 SelectVariants on '{input.vcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        gatk SelectVariants \\
            -R {params.genome} \\
            --intervals {params.targets} \\
            --variant {input.vcf} \\
            --sample-name {params.Sname} \\
            --exclude-filtered \\
            --exclude-non-variants \\
            --output {output.vcf}
        """
