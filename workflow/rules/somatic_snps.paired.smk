# Somatic SNP calling rules for tumor/normal pairs
rule gatk_mutect2:
    input:
        normal = lambda w: [os.path.join(output_bamdir, "chrom_split", pairs_dict[w.samples] + ".{chroms}.split.bam")],
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam")
    output:
        vcf = os.path.join(output_somatic_snpindels,"mutect2_out", "chrom_split", "{samples}.{chroms}.vcf"),
        read_orientation_file = os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{samples}.{chroms}.f1r2.tar.gz"),
        statsfiles = os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{samples}.{chroms}.vcf.stats")
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        chrom = '{chroms}',
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        germsource = config['references']['GERMLINERESOURCE'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'mutect2'
    threads: 2
    envmodules:
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    shell: """
    if [ ! -d "$(dirname {output.vcf})" ];
        then mkdir -p "$(dirname {output.vcf})";
    fi
    gatk Mutect2 \\
        -R {params.genome} \\
        -I {input.tumor} \\
        -I {input.normal} \\
        -normal {params.normalsample} \\
        --panel-of-normals {params.pon} \\
        {params.germsource} \\
        -L {params.chrom} \\
        -O {output.vcf} \\
        --f1r2-tar-gz {output.read_orientation_file} \\
        --independent-mates
    """


rule pileup_paired:
    input:
        tumor = os.path.join(output_bamdir, "final_bams", "{samples}.bam"),
        tumorbai = os.path.join(output_bamdir, "final_bams", "{samples}.bai"),
        normal = lambda w: [os.path.join(output_bamdir, "final_bams", pairs_dict[w.samples] + ".bam")],
        normalbai = lambda w: [os.path.join(output_bamdir, "final_bams", pairs_dict[w.samples] + ".bai")],
        intervals = intervals_file
    output:
        tumor_summary = os.path.join(output_somatic_snpindels, "mutect2_out", "pileup_summaries", "{samples}_tumor.pileup.table"),
        normal_summary = os.path.join(output_somatic_snpindels, "mutect2_out", "pileup_summaries", "{samples}_normal.pileup.table")
    params:
        genome = config['references']['GENOME'],
        germsource = config['references']['KNOWNSNPS'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'pileup'
    envmodules:
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Run GetPileupSummaries in bg concurrently for a tumor/normal pair
    gatk --java-options '-Xmx48g' GetPileupSummaries \\
        -I {input.tumor} \\
        -V {params.germsource} \\
        -L {input.intervals} \\
        -O {output.tumor_summary} & \\
    gatk --java-options '-Xmx48g' GetPileupSummaries \\
        -I {input.normal} \\
        -V {params.germsource} \\
        -L {input.intervals} \\
        -O {output.normal_summary} & \\
    wait
    """

rule contamination_paired:
    input:
        tumor = os.path.join(output_somatic_snpindels, "mutect2_out", "pileup_summaries", "{samples}_tumor.pileup.table"),
        normal = os.path.join(output_somatic_snpindels, "mutect2_out", "pileup_summaries", "{samples}_normal.pileup.table"),
    output:
        tumor_summary = os.path.join(output_somatic_base, "qc", "gatk_contamination", "{samples}.contamination.table"),
        normal_summary = os.path.join(output_somatic_base, "qc", "gatk_contamination", "{samples}_normal.contamination.table")
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'contamination'
    envmodules:
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    shell: """
    gatk CalculateContamination \\
        -I {input.tumor} \\
        --matched-normal {input.normal} \\
        -O {output.tumor_summary}
    gatk CalculateContamination \\
        -I {input.normal} \\
        -O {output.normal_summary}
    """



rule strelka:
    input:
        normal = lambda w: [os.path.join(output_bamdir, "chrom_split", pairs_dict[w.samples] + ".{chroms}.split.bam")],
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam")
    output:
        vcf = os.path.join(output_somatic_snpindels, "strelka_out", "chrom_split", "{samples}.{chroms}.vcf"),
    params:
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        basedir = BASEDIR,
        ver_strelka = config['tools']['strelka']['version'],
        rname = 'strelka',
        set_tmp = set_tmp(),
    envmodules:
        config['tools']['strelka']['modname'],
        config['tools']['gatk3']['modname'],
        config['tools']['java']['modname']
    container:
        config['images']['wes_base']
    threads: 16
    envmodules:
        config['tools']['strelka']['modname'],
        config['tools']['gatk3']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    workdir={params.basedir}
    myoutdir="$(dirname {output.vcf})/{wildcards.samples}/{wildcards.chroms}"
    if [ -d "$myoutdir" ]; then rm -r "$myoutdir"; fi
    mkdir -p "$myoutdir"

    configureStrelkaSomaticWorkflow.py \\
        --ref={params.genome} \\
        --tumor={input.tumor} \\
        --normal={input.normal} \\
        --runDir="$myoutdir" \\
        --exome
    cd "$myoutdir"
    ./runWorkflow.py -m local -j {threads}

    java8 -Xmx12g -Djava.io.tmpdir=${{tmp}} -XX:ParallelGCThreads={threads} \\
        -jar $GATK_JAR -T CombineVariants \\
        -R {params.genome} \\
        --variant results/variants/somatic.snvs.vcf.gz \\
        --variant results/variants/somatic.indels.vcf.gz \\
        --assumeIdenticalSamples \\
        --filteredrecordsmergetype KEEP_UNCONDITIONAL \\
        -o "$(basename {output.vcf})"

    cd $workdir
    mv "$myoutdir/$(basename {output.vcf})" "{output.vcf}"
    """


rule strelka_filter:
    input:
        vcf = os.path.join(output_somatic_snpindels, "strelka_out", "vcf", "{samples}.collected.vcf"),
    output:
        filtered = temp(os.path.join(output_somatic_snpindels, "strelka_out", "vcf", "{samples}.filtered.vcf")),
        samplesfile = temp(os.path.join(output_somatic_snpindels, "strelka_out", "vcf", "{samples}.FINAL.vcf.samples")),
        final = os.path.join(output_somatic_snpindels, "strelka_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "strelka_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        basedir = BASEDIR,
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'strelka_filter',
        set_tmp = set_tmp()
    threads: 4
    envmodules:
        config['tools']['gatk4']['modname'],
        config['tools']['bcftools']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {input.vcf} \\
        --discordance {params.pon} \\
        --exclude-filtered \\
        --output {output.filtered}

    echo -e "TUMOR\t{params.tumorsample}\nNORMAL\t{params.normalsample}" > "{output.samplesfile}"

    echo "Reheading VCFs with sample names..."
    bcftools reheader \\
        -o "{output.final}" \\
        -s "{output.samplesfile}" "{output.filtered}"

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools sort \\
        -T ${{tmp}} "{output.final}" \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """


rule mutect_paired:
    input:
        normal = lambda w: [os.path.join(output_bamdir, "chrom_split", pairs_dict[w.samples] + ".{chroms}.split.bam")],
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
    output:
        vcf = os.path.join(output_somatic_snpindels, "mutect_out", "chrom_split", "{samples}.{chroms}.vcf"),
        stats = os.path.join(output_somatic_snpindels, "mutect_out", "chrom_split", "{samples}.{chroms}.stats.out"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        pon = config['references']['PON'],
        genome = config['references']['GENOME'],
        dbsnp_cosmic = config['references']['DBSNP_COSMIC'],
        ver_mutect = config['tools']['mutect']['version'],
        rname = 'mutect',
        set_tmp = set_tmp(),
    envmodules:
        config['tools']['mutect']['modname']
    container:
        config['images']['mutect']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    if [ ! -d "$(dirname {output.vcf})" ]; then mkdir -p "$(dirname {output.vcf})"; fi

    java -Xmx8g -Djava.io.tmpdir=${{tmp}} -jar ${{MUTECT_JAR}} \\
        --analysis_type MuTect \\
        --reference_sequence {params.genome} \\
        --normal_panel {params.pon} \\
        --vcf {output.vcf} \\
        {params.dbsnp_cosmic} \\
        --disable_auto_index_creation_and_locking_when_reading_rods \\
        --input_file:normal {input.normal} \\
        --input_file:tumor {input.tumor} \\
        --out {output.stats} \\
        -rf BadCigar
    """


rule mutect_filter:
    input:
        vcf = os.path.join(output_somatic_snpindels, "mutect_out", "vcf", "{samples}.collected.vcf"),
    output:
        final = os.path.join(output_somatic_snpindels, "mutect_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "mutect_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'mutect_filter',
        set_tmp = set_tmp(),
    threads: 4
    envmodules:
        config['tools']['gatk4']['modname'],
        config['tools']['bcftools']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {input.vcf} \\
        --exclude-filtered \\
        --output {output.final}

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools view -s "{params.tumorsample}","{params.normalsample}" {output.final} \\
        | bcftools sort -T ${{tmp}} - \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """


rule vardict_paired:
    input:
        normal = lambda w: [os.path.join(output_bamdir, "chrom_split", pairs_dict[w.samples] + ".{chroms}.split.bam")],
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
    output:
        vcf = os.path.join(output_somatic_snpindels, "vardict_out", "chrom_split", "{samples}.{chroms}.vcf"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = "{samples}",
        genome = config['references']['GENOME'],
        targets = exome_targets_bed,
        pon = config['references']['PON'],
        rname = 'vardict'
    envmodules:
        config['tools']['R']['modname'],
        config['tools']['samtools']['modname']
    container:
        config['images']['wes_base']
    shell: """
    if [ ! -d "$(dirname {output.vcf})" ]; then mkdir -p "$(dirname {output.vcf})"; fi
    VarDict \\
        -G {params.genome} \\
        -f 0.05 \\
        -N \"{params.tumorsample}|{params.normalsample}\" \\
        --nosv \\
        -b \"{input.tumor}|{input.normal}\" \\
        -t \\
        -Q 20 \\
        -c 1 \\
        -S 2 \\
        -E 3 {params.targets} \\
        | testsomatic.R \\
        | var2vcf_paired.pl \\
            -S \\
            -Q 20 \\
            -d 10 \\
            -M \\
            -N \"{params.tumorsample}|{params.normalsample}\" \\
            -f 0.05 > {output.vcf}
    """


rule vardict_filter:
    input:
        vcf = os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.collected.vcf"),
    output:
        final = os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.FINAL.vcf"),
        filtered = temp(os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.filtered.vcf")),
        norm = os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        targets = exome_targets_bed,
        pon = config['references']['PON'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'vardict',
        set_tmp = set_tmp(),
    threads: 4
    envmodules:
        config['tools']['bcftools']['modname'],
        config['tools']['gatk4']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    bcftools filter \\
        --exclude \'STATUS=\"Germline\" | STATUS=\"LikelyLOH\" | STATUS=\"AFDiff\"\' \\
        {input.vcf} > {output.filtered}

    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {output.filtered} \\
        --discordance {params.pon} \\
        --exclude-filtered \\
        --output {output.final}

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools sort -T ${{tmp}} "{output.final}" \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """


rule varscan_paired:
    """Note: Refactor formatting of shell command for readability to
    be more snake-thonic."""
    input:
        normal = lambda w: [os.path.join(output_bamdir, "chrom_split", pairs_dict[w.samples] + ".{chroms}.split.bam")],
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
        tumor_summary = os.path.join(output_somatic_base, "qc", "gatk_contamination", "{samples}.contamination.table"),
        normal_summary = lambda w: [os.path.join(output_somatic_base, "qc", "gatk_contamination", "{samples}_normal.contamination.table")],
    output:
        vcf = os.path.join(output_somatic_snpindels, "varscan_out", "chrom_split", "{samples}.{chroms}.vcf"),
    params:
        genome = config['references']['GENOME'],
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        ver_varscan = config['tools']['varscan']['version'],
        rname = 'varscan',
        set_tmp = set_tmp(),
    threads: 4
    envmodules:
        config['tools']['varscan']['modname'],
        config['tools']['gatk3']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    if [ ! -d "$(dirname {output.vcf})" ]; then mkdir -p "$(dirname {output.vcf})"; fi

    tumor_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 {input.tumor_summary} | cut -f2 ))" | bc -l)
    normal_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 {input.normal_summary} | cut -f2 ))" | bc -l)
    varscan_opts="--strand-filter 1 --min-var-freq 0.01 --min-avg-qual 30 --somatic-p-value 0.05 --output-vcf 1 --normal-purity $normal_purity --tumor-purity $tumor_purity"
    dual_pileup="samtools mpileup -d 10000 -q 15 -Q 15 -f {params.genome} {input.normal} {input.tumor}"
    varscan_cmd="varscan somatic <($dual_pileup) {output.vcf} $varscan_opts --mpileup 1"
    eval "$varscan_cmd"

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    awk '{{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}}' {output.vcf}.snp \\
        | sed '/^$/d' > {output.vcf}.snp_temp
    awk '{{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}}' {output.vcf}.indel \\
        | sed '/^$/d' > {output.vcf}.indel_temp

    java8 -Xmx12g -Djava.io.tmpdir=${{tmp}} -XX:ParallelGCThreads={threads} \\
        -jar $GATK_JAR -T CombineVariants \\
        -R {params.genome} \\
        --variant {output.vcf}.snp_temp \\
        --variant {output.vcf}.indel_temp \\
        --assumeIdenticalSamples \\
        --filteredrecordsmergetype KEEP_UNCONDITIONAL \\
        -o {output.vcf}
    """


rule varscan_filter:
    input:
        vcf = os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.collected.vcf"),
    output:
        filtered = temp(os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.filtered.vcf")),
        filtered1 = temp(os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.filtered.1.vcf")),
        samplesfile = temp(os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.FINAL.vcf.samples")),
        final = os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        basedir = BASEDIR,
        somatic_filter_settings = config['tools']['varscan']['somatic_filter_settings'],
        ver_varscan = config['tools']['varscan']['version'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'varscan_filter',
        set_tmp = set_tmp(),
    threads: 4
    envmodules:
        config['tools']['varscan']['modname'],
        config['tools']['gatk4']['modname'],
        config['tools']['bcftools']['modname']
    container:
        config['images']['wes_base']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    {params.set_tmp}

    varscan somaticFilter \\
        {input.vcf} \\
        {params.somatic_filter_settings} --output-file {output.filtered1}

    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {output.filtered1} \\
        --discordance {params.pon} \\
        --exclude-filtered \\
        --output {output.filtered}

    samplesFile="{output.samplesfile}"
    echo -e "TUMOR\t{params.tumorsample}\nNORMAL\t{params.normalsample}" > "{output.samplesfile}"

    bcftools reheader \\
        -o "{output.final}" \\
        -s "{output.samplesfile}" \\
        "{output.filtered}"

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools sort -T ${{tmp}} "{output.final}" \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """
