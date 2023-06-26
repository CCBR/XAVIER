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
    envmodules: config['tools']['gatk4']['modname']
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
    envmodules: config['tools']['gatk4']['modname']
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
    envmodules: config['tools']['gatk4']['modname']
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
    envmodules: config['tools']['gatk4']['modname']
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


rule Gatk_Variantfilter: 
    """
    Hard filters on
    @Input:
        Multi-sample gVCF with all chromosomes combined
    @Output:
       Variants filtered by QD, QUAL, SOR, FS, MQ, MQRankSum, ReadPosRankSum, Indels: QD, QUAL, FS, ReadPosRankSum
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output: 
       indelvcf = os.path.join(output_germline_base,"VCF","indel.filterd.vcf.gz"),
       snpvcf = os.path.join(output_germline_base,"VCF","snp.filtered.vcf.gz"),
       vcf = os.path.join(output_germline_base,"VCF","snp_indel.filtered.vcf.gz")

    params: 
        genome=config['references']['GENOME'], 
        rname="gatk_hardfilters",
        ver_gatk=config['tools']['gatk4']['version']
    message: "Running GATK4 hard filters on Cohort VCF input file"
    envmodules: config['tools']['gatk4']['modname']
    container: config['images']['wes_base']
    shell:
        """
        gatk SelectVariants \\
        -V {input.vcf} \\
        -select-type SNP \\
        -O snps.vcf.gz
        
        gatk SelectVariants \\
        -V {input.vcf} \\
        -select-type INDEL \\
        -O indels.vcf.gz

        gatk VariantFiltration \\
        -V snps.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "SOR > 3.0" --filter-name "SOR3" \\
        -filter "FS > 60.0" --filter-name "FS60" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -O {output.snpvcf}

        gatk VariantFiltration \\
        -V indels.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "FS > 200.0" --filter-name "FS200" \\
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        -O {output.indelvcf}
        
        gatk MergeVcfs \\
        -R {params.genome} \\
        --INPUT {output.indelvcf} \\
        --INPUT {output.snpvcf} \\
        --OUTPUT {output.vcf}
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
        vcf = os.path.join(output_germline_base,"VCF","snp_indel.filtered.vcf.gz"),
    output: 
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz")
    params: 
        genome=config['references']['GENOME'], 
        Sname = "{samples}", 
        rname="varselect",
        ver_gatk=config['tools']['gatk4']['version'],
        targets=exome_targets_bed
    message: "Running GATK4 SelectVariants on '{input.vcf}' input file"
    envmodules: config['tools']['gatk4']['modname']
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