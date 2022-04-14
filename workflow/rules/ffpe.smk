# Rules for correcting stand orientation bias in FFPE samples
rule sobdetect_get:
    input: 
    output: 
        SOBDetector_jar = SOBDetector_JARFILE
    params:
        rname = 'get_sobdetector'
    shell: """
    wget https://github.com/mikdio/SOBDetector/releases/download/v1.0.2/SOBDetector_v1.0.2.jar \\
        -O {output.SOBDetector_jar}
    """


rule sobdetect_pass1:
    input:
        vcf = os.path.join(output_somatic_snpindels, "{vc_outdir}", "vcf", "{samples}.FINAL.norm.vcf"),
        bam = os.path.join(output_bamdir, "final_bams", "{samples}.bam"),
        SOBDetector_jar = SOBDetector_JARFILE
    output:
        pass1_vcf = os.path.join(SOBDetector_out, "{vc_outdir}", "pass1", "{samples}.sobdetect.vcf"),
        pass1_info = os.path.join(SOBDetector_out, "{vc_outdir}", "pass1", "{samples}.info")
    params:
        chrom = chroms,
        rname = 'sobdetect1'
    envmodules:
        'samtools/1.12',
        'bcftools/1.9'
    container:
       config['images']['wes_base'] 
    shell: """
    if [ ! -d "$(dirname {output.pass1_vcf})" ]; then 
        mkdir -p "$(dirname {output.pass1_vcf})"
    fi

    echo "Running SOBDetector..."
    # Try/catch for running SOB Dectetor
    # with an empty input VCF file 
    java -jar {input.SOBDetector_jar} \\
        --input-type VCF \\
        --input-variants "{input.vcf}" \\
        --input-bam {input.bam} \\
        --output-variants {output.pass1_vcf} \\
        --only-passed false || {{
    # Compare length of VCF header to 
    # the total length of the file
    header_length=$(grep '^#' "{input.vcf}" | wc -l)
    file_length=$(cat "{input.vcf}" | wc -l)
    if [ $header_length -eq $file_length ]; then
        # VCF file only contains header
        # File contains no variants, catch
        # problem so pipeline can continue
        cat "{input.vcf}" > {output.pass1_vcf}
    else
        # SOB Dectector failed for another reason
        echo "SOB Detector Failed... exiting now!" 1>&2
        exit 1
    fi
    }}

    bcftools query \\
        -f '%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\n' \\
        {output.pass1_vcf} \\
        | awk '{{if ($1 != "."){{tum_alt=$1+$2; tum_depth=$1+$2+$3+$4+$5+$6; if (tum_depth==0){{tum_af=1}} else {{tum_af=tum_alt/tum_depth }}; print tum_alt,tum_depth,tum_af,$7}}}}' > {output.pass1_info} 
    """


rule sobdetect_cohort_params:
    input:
        info_files = expand(os.path.join(SOBDetector_out, "{{vc_outdir}}", "pass1", "{samples}.info"), samples=ffpe_sample_list)
    output:
        all_info_file = os.path.join(SOBDetector_out, "{vc_outdir}", "pass1", "all_samples.info"),
        params_file = os.path.join(SOBDetector_out, "{vc_outdir}", "cohort_params.txt")
    params:
        rname = 'sobdetect_params'
    container:
       config['images']['wes_base'] 
    shell: """
    echo -e "#TUMOR.alt\\tTUMOR.depth\\tTUMOR.AF\\tSOB\\tFS\\tSOR\\tTLOD\\tReadPosRankSum" > {output.all_info_file}
    cat {input.info_files} >> {output.all_info_file}
    
    # Try/catch for running calculating
    # mean and standard deviation with 
    # with a set of empty input VCF files
    all_length=$(tail -n+2 {output.all_info_file} | wc -l)
    if [ $all_length -eq 0 ]; then 
        echo 'WARNING: All SOB Dectect pass1 samples contained no variants.' \\
        | tee {output.params_file}
    else
        # Calculate mean and standard deviation
        grep -v '^#' {output.all_info_file} \\
        | awk '{{ total1 += $1; ss1 += $1^2; total2 += $2; ss2 += $2^2; total3 += $3; ss3 += $3^2; total4 += $4; ss4 += $4^2 }} END {{ print total1/NR,total2/NR,total3/NR,total4/NR; print sqrt(ss1/NR-(total1/NR)^2),sqrt(ss2/NR-(total2/NR)^2),sqrt(ss3/NR-(total3/NR)^3),sqrt(ss4/NR-(total4/NR)^2) }}' > {output.params_file}
    fi
    """

  
rule sobdetect_pass2:
    input:
        vcf = os.path.join(output_somatic_snpindels, "{vc_outdir}", "vcf", "{samples}.FINAL.norm.vcf"),
        bam = os.path.join(output_bamdir, "final_bams", "{samples}.bam"),
        SOBDetector_jar = SOBDetector_JARFILE,
        params_file = os.path.join(SOBDetector_out, "{vc_outdir}", "cohort_params.txt")
    output:
        pass2_vcf = os.path.join(SOBDetector_out, "{vc_outdir}", "pass2", "{samples}.sobdetect.vcf"),
        pass2_info = os.path.join(SOBDetector_out, "{vc_outdir}", "pass2", "{samples}.info"),
        filtered_vcf = os.path.join(SOBDetector_out, "{vc_outdir}", "pass2", "{samples}.artifact_filtered.vcf.gz")
    params:
        chrom=chroms,
        ver_bcftools=config['tools']['bcftools']['version'],
        rname="sobdetect2",
    threads: 4
    envmodules:
        'samtools/1.12',
        'bcftools/1.9'
    container:
       config['images']['wes_base'] 
    shell: """
    if [ ! -d "$(dirname {output.pass2_vcf})" ]; then 
        mkdir -p "$(dirname {output.pass2_vcf})"
    fi

    echo "Running SOBDetector..."
    # Try/catch for running SOB Dectetor
    # with an empty input VCF file
    bcf_annotate_option="-e 'INFO/pArtifact < 0.05' "
    java -jar {input.SOBDetector_jar} \\
        --input-type VCF \\
        --input-variants "{input.vcf}" \\
        --input-bam "{input.bam}" \\
        --output-variants "{output.pass2_vcf}" \\
        --only-passed true \\
        --standardization-parameters "{input.params_file}" || {{
    # Compare length of VCF header to 
    # the total length of the file
    header_length=$(grep '^#' "{input.vcf}" | wc -l)
    file_length=$(cat "{input.vcf}" | wc -l)
    if [ $header_length -eq $file_length ]; then
        # VCF file only contains header
        # File contains no variants, catch
        # problem so pipeline can continue
        cat "{input.vcf}" > {output.pass2_vcf}
    else
        # SOB Dectector failed for another reason
        echo "SOB Detector Failed... exiting now!" 1>&2
        exit 1
    fi
    }}
    
    echo "Making info table..."
    bcftools query \\
        -f '%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\n' \\
        "{output.pass2_vcf}" \\
        | awk '{{if ($1 != "."){{tum_alt=$1+$2; tum_depth=$1+$2+$3+$4+$5+$6; if (tum_depth==0){{tum_af=1}} else {{tum_af=tum_alt/tum_depth }}; print tum_alt,tum_depth,tum_af,$7}}}}' > "{output.pass2_info}"

    echo "Filtering out artifacts..."
    if [ "{wildcards.vc_outdir}" == "{config[output_params][MERGED_SOMATIC_OUTDIR]}" ]; then
        echo "Adding 'set' annotation back from merged variants..."
        bgzip --threads {threads} -c "{input.vcf}" > "{input.vcf}.gz"
        bcftools index -f -t "{input.vcf}.gz"
        bgzip --threads {threads} -c "{output.pass2_vcf}" > "{output.pass2_vcf}.gz"
        bcftools index -f -t "{output.pass2_vcf}.gz"
        bcftools annotate \\
            -a "{input.vcf}.gz" \\
            -c "INFO/set" \\
            "$bcf_annotate_option" \\
            -Oz \\
            -o {output.filtered_vcf} {output.pass2_vcf}.gz
    else
        bcftools filter \\
            "$bcf_annotate_option" \\
            -Oz \\
            -o {output.filtered_vcf} {output.pass2_vcf}
        bcftools index -f -t {output.filtered_vcf}
    fi
    """


rule sobdetect_metrics:
    input:
        pass1_vcf = expand(os.path.join(SOBDetector_out, "{{vc_outdir}}", "pass1", "{samples}.sobdetect.vcf"), samples=ffpe_sample_list),
        pass2_vcf = expand(os.path.join(SOBDetector_out, "{{vc_outdir}}", "pass2", "{samples}.sobdetect.vcf"), samples=ffpe_sample_list)
    output:
        count_table = os.path.join(SOBDetector_out, "{vc_outdir}", "metrics", "variant_count_table.txt"),
        full_metric_table = os.path.join(SOBDetector_out, "{vc_outdir}", "metrics", "all_metrics.txt")
    params:
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'sobdetect_metrics',
    envmodules:
        'bcftools/1.9'
    container:
        config['images']['wes_base'] 
    shell: """
    echo -e "#ID\\tDefaultParam\\tCohortParam\\tTotalVariants" > {output.count_table}
    echo -e "#SAMPLE_ID\\tParam\\tCHROM\\tPOS\\tnumF1R2Alt\\tnumF2R1Alt\\tnumF1R2Ref\\tnumF2R1Ref\\tnumF1R2Other\\tnumF2R1Other\\tSOB\\tpArtifact\\tFS\\tSOR\\tTLOD\\tReadPosRankSum" > {output.full_metric_table}
    
    P1FILES=({input.pass1_vcf})
    P2FILES=({input.pass2_vcf})
    for (( i=0; i<${{#P1FILES[@]}}; i++ )); do
        MYID=$(basename -s ".sobdetect.vcf" ${{P1FILES[$i]}})
        echo "Collecting metrics from $MYID..."

        # grep may fail if input files do not contain any variants 
        total_count=$(grep -v ^# ${{P1FILES[$i]}} | wc -l) || total_count=0
        count_1p=$(bcftools query -f '%INFO/pArtifact\n' ${{P1FILES[$i]}} | awk '{{if ($1 != "." && $1 < 0.05){{print}}}}' | wc -l)
        count_2p=$(bcftools query -f '%INFO/pArtifact\n' ${{P2FILES[$i]}} | awk '{{if ($1 != "." && $1 < 0.05){{print}}}}' | wc -l)
     
        echo -e "$MYID\\t$count_1p\\t$count_2p\\t$total_count" >> {output.count_table}
      
        bcftools query -f '%CHROM\\t%POS\\t%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\t%INFO/pArtifact\n' ${{P1FILES[$i]}} | awk -v id=$MYID 'BEGIN{{OFS="\t"}}{{print id,"PASS_1",$0}}' >> {output.full_metric_table}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\t%INFO/pArtifact\n' ${{P2FILES[$i]}} | awk -v id=$MYID 'BEGIN{{OFS="\t"}}{{print id,"PASS_2",$0}}' >> {output.full_metric_table}
    done
    """


rule ffpefilter_mafs:
    input:
        filtered_vcf = os.path.join(SOBDetector_out, "{vc_outdir}", "pass2", "{samples}.artifact_filtered.vcf.gz")
    output:
        maf = os.path.join(output_somatic_base, SOBDetector_out, "{vc_outdir}", "maf", "{samples}.maf")
    params:
        tumorsample = '{samples}',
        genome = config['references']['MAF_GENOME'],
        filtervcf = config['references']['MAF_FILTERVCF'],
        bundle = config['references']['VCF2MAF']['VEPRESOURCEBUNDLEPATH'],
        rname = 'vcf2maf',
        vcf2maf_script = VCF2MAF_WRAPPER
    threads: 4
    container:
        config['images']['vcf2maf'] 
    shell: """
    echo "Converting to MAF..."
    bash {params.vcf2maf_script} \\
        --vcf {input.filtered_vcf} \\
        --maf {output.maf} \\
        --tid {params.tumorsample} \\
        --genome {params.genome} \\
        --threads {threads} \\
        --vepresourcebundlepath {params.bundle} \\
        --info "set"
    echo "Done converting to MAF..."
    """


rule collect_ffpefilter_mafs:
    input:
        mafs = expand(os.path.join(SOBDetector_out, "{{vc_outdir}}", "maf", "{samples}"+".maf"), samples=ffpe_sample_list)
    output:
        maf = os.path.join(output_somatic_base, SOBDetector_out, "{vc_outdir}", "cohort_summary", "all_somatic_variants.maf")
    params:
        rname = "combine_maf"
    container:
        config['images']['wes_base'] 
    shell: """    
    echo "Combining MAFs..."
    head -2 {input.mafs[0]} > {output.maf}
    awk 'FNR>2 {{print}}' {input.mafs} >> {output.maf}
    """