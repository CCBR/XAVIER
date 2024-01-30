# Rules for primary processing of raw data: trim, align, and recal
rule bam2fastq:
    """
    Convert BAM files to paired FASTQ.
    A few of the QC tools run directly on fastq files (kraken and
    fastqscreen, maybe others?).  When starting the pipeline from BAMs,
    this rule ensures that a fastq is also available.  The 'ancient' tag
    on the input BAM ensures that this rule is not run unnecessarily when
    re-starting the pipeline after an error.
    @Input:
        BAM file (scatter)
    @Output:
        Paired FASTQ files
    """
    input:
        bam = ancient(os.path.join(input_bamdir, "{samples}.input.bam")),
    output:
        r1 = os.path.join(input_fqdir, "{samples}.R1.fastq.gz"),
        r2 = os.path.join(input_fqdir, "{samples}.R2.fastq.gz"),
        orphans = temp(os.path.join(input_fqdir, "{samples}.orphans.fastq.gz")),
    params:
        genome = config['references']['GENOME'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'bam2fastq',
        set_tmp = set_tmp()
    envmodules:
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    mkdir -p fastqs
    gatk SamToFastq \\
        --INPUT {input.bam} \\
        --FASTQ {output.r1} \\
        --SECOND_END_FASTQ {output.r2} \\
        --UNPAIRED_FASTQ {output.orphans} \\
        --TMP_DIR ${{tmp}} \\
        -R {params.genome}
    """


rule trimmomatic:
    """
    Data-processing step to remove adapter sequences and perform quality trimming
    prior to alignment the reference genome.  Adapters are composed of synthetic
    sequences and should be removed prior to alignment.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Trimmed FastQ file
    """
    input:
        r1 = os.path.join(input_fqdir, "{samples}.R1.fastq.gz"),
        r2 = os.path.join(input_fqdir, "{samples}.R2.fastq.gz")
    output:
        one = temp(os.path.join(output_fqdir, "{samples}.R1.trimmed.fastq.gz")),
        two = temp(os.path.join(output_fqdir, "{samples}.R1.trimmed.unpair.fastq.gz")),
        three = temp(os.path.join(output_fqdir, "{samples}.R2.trimmed.fastq.gz")),
        four = temp(os.path.join(output_fqdir, "{samples}.R2.trimmed.unpair.fastq.gz"))
    params:
        adapterfile = config['references']['trimmomatic.adapters'],
        ver = config['tools']['trimmomatic']['version'],
        rname = 'trimmomatic'
    envmodules:
        config['tools']['trimmomatic']['modname']
    container:
        config['images']['wes_base']
    threads: 24
    shell: """
    myoutdir="$(dirname {output.one})"
    if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
    trimmomatic PE \\
        -threads {threads} \\
        -phred33 \\
        {input.r1} {input.r2} \\
        {output.one} {output.two} \\
        {output.three} {output.four} \\
        ILLUMINACLIP:{params.adapterfile}:3:30:10 \\
        LEADING:10 \\
        TRAILING:10 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:20
    """


rule bwa_mem:
    """
    Map trimmed paired reads to the reference genome using the 'bwa' aligner.
    @Input:
        One pair of FASTQ files (scatter)
    @Output:
        Aligned reads in BAM format
    """
    input:
        os.path.join(output_fqdir, "{samples}.R1.trimmed.fastq.gz"),
        os.path.join(output_fqdir, "{samples}.R2.trimmed.fastq.gz")
    output:
        temp(os.path.join(output_bamdir, "preprocessing", "{samples}.raw_map.bam"))
    params:
        genome = config['references']['BWAGENOME'],
        sample = "{samples}",
        ver_samtools = config['tools']['samtools']['version'],
        ver_bwa = config['tools']['bwa']['version'],
        rname = 'bwamem'
    envmodules:
       config['tools']['samtools']['modname'],
       config['tools']['bwa']['modname'],
       config['tools']['samblaster']['modname']
    container:
        config['images']['wes_base']
    threads: 24
    shell: """
    myoutdir="$(dirname {output})"
    if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
    bwa mem -M \\
        -R \'@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:illumina\\tLB:{params.sample}\\tPU:{params.sample}\\tCN:hgsc\\tDS:wes\' \\
        -t {threads} \\
        {params.genome} \\
        {input} | \\
    samblaster -M | \\
    samtools sort -@12 -m 4G - -o {output}
    """

localrules: raw_index
rule raw_index:
    """
    Make index of mapped BAM file.
    Honestly, it makes more sense to bundle this with the map step in most cases.
    The reason to keep it separate is if you want to allow custom mapped files (before
    recalibration) to be input into the pipeline. But we're not doing that here.
    @Input:
        One pair of FASTQ files (scatter)
    @Output:
        Aligned reads in BAM format
    """
    input:
        bam = os.path.join(output_bamdir, "preprocessing", "{samples}.raw_map.bam")
    output:
        bai = temp(os.path.join(output_bamdir,"preprocessing","{samples}.raw_map.bai")),
    params:
        ver_samtools = config['tools']['samtools']['version'],
        rname = 'raw_index'
    envmodules:
        config['tools']['samtools']['modname']
    container:
        config['images']['wes_base']
    shell: """
    samtools index -@ 2 {input.bam} {output.bai}
    """


rule gatk_recal:
    """
    Base quality recalibration (BQSR), part of the GATK Best Practices.
    The idea is that each sequencer/run will have systematic biases.
    BQSR learns about these biases using known sites of variation (common SNPs),
    and uses it to adjust base quality on all sites, including novel sites of
    variation.  Since base quality is taken into account during variant calling,
    this will help pick up real variants in low depth or otherwise noisy loci.
    @Input:
        Aligned reads in BAM format (scatter)
    @Output:
        Aligned reads in BAM format, with altered quality scores
    """
    input:
        bam = os.path.join(output_bamdir, "preprocessing", "{samples}.raw_map.bam"),
        bai = os.path.join(output_bamdir, "preprocessing", "{samples}.raw_map.bai"),
    output:
        bam = os.path.join(input_bamdir, "{samples}.input.bam"),
        re = temp(os.path.join(output_bamdir, "preprocessing", "{samples}_recal_data.grp"))
    params:
        genome = config['references']['GENOME'],
        knowns = config['references']['KNOWNRECAL'],
        ver_gatk = config['tools']['gatk4']['version'],
        chrom = chroms,
        intervals = intervals_file,
        rname = 'recal'
    envmodules:
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    threads: 24
    shell: """
    gatk --java-options '-Xmx48g' BaseRecalibrator \\
        --input {input.bam} \\
        --reference {params.genome} \\
        {params.knowns} \\
        --output {output.re} \\
        --intervals {params.intervals}

    gatk --java-options '-Xmx48g' ApplyBQSR \\
        --reference {params.genome} \\
        --input {input.bam} \\
        --bqsr-recal-file {output.re} \\
        --output {output.bam} \\
        --use-jdk-inflater \\
        --use-jdk-deflater
    """


rule bam_check:
    """
    This is a checkpoint to make sure BAMs are ready for variant calling.
    The read group (RG) tags are checked to make sure they match the sample ID
    inferred from the file name, and the bam is indexed. This rule needs to be
    refactored at some point. It is not making great use of snake-thonic.
    @Input:
        Aligned reads in BAM format (scatter)
    @Output:
        Aligned reads in BAM format, with appropriate RG tags and index file
    """
    input:
        bam = os.path.join(input_bamdir, "{samples}.input.bam")
    output:
        bam = os.path.join(output_bamdir, "final_bams", "{samples}.bam"),
        bai = os.path.join(output_bamdir, "final_bams", "{samples}.bai"),
        bai2 = os.path.join(output_bamdir, "final_bams", "{samples}.bam.bai"),
    params:
        ver_samtools = config['tools']['samtools']['version'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'bam_check'
    envmodules:
        config['tools']['samtools']['modname'],
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    shell: """
    sample={wildcards.samples}
    ID=$sample
    PL="ILLUMINA"  # exposed as a config param
    LB="na"        # exposed as a config param

    # Check if there is no header or any of the info
    HEADER=`samtools view -H {input.bam} | grep ^@RG`
    if [[ "$HEADER" != "" ]]; then
        t=(${{HEADER//\t/ }})
        echo ${{t[1]}}
        ID=`printf '%s\n' "${{t[@]}}" | grep -P '^ID' | cut -d":" -f2` #(${{t[1]//:/ }})
        PL=`printf '%s\n' "${{t[@]}}" | grep -P '^PL' | cut -d":" -f2` #(${{t[3]//:/ }})
        LB=`printf '%s\n' "${{t[@]}}" | grep -P '^LB' | cut -d":" -f2` #(${{t[2]//:/ }})
        if [[ "$ID" == "$sample" ]]; then
            echo "The header of the BAM file is correct"
        else
            ID=$sample
        fi
    fi

    gatk AddOrReplaceReadGroups \\
        --INPUT {input.bam} \\
        --OUTPUT {output.bam} \\
        --RGID ${{ID}} \\
        --RGLB ${{LB}} \\
        --RGPL ${{PL}} \\
        --RGSM ${{ID}} \\
        --RGPU na

    samtools index -@ 2 {output.bam} {output.bai}
    cp {output.bai} {output.bai2}
    """
