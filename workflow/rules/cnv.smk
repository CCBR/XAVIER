# Rules to predict copy number variation
rule cnvkit_amplicon:
    input:
        normal = lambda w: [os.path.join(output_bamdir,"final_bams", pairs_dict[w.samples] + ".bam")],
        tumor = os.path.join(output_bamdir,"final_bams","{samples}.bam")
    output:
        call = os.path.join(output_somatic_cnv, "cnvkit", "{samples}.call.cns"),
        bintest = os.path.join(output_somatic_cnv, "cnvkit", "{samples}.bintest.cns"),
        cnr = os.path.join(output_somatic_cnv, "cnvkit", "{samples}.cnr"),
        cns = os.path.join(output_somatic_cnv, "cnvkit", "{samples}.cns")
    params:
        targets=exome_targets_bed,
        outdir = os.path.join(output_somatic_cnv,"cnvkit"),
        genome=config['references']['GENOME'], 
        refFlat=config['references']['REFFLAT'],
        rname = 'cnvkit_amplicon'
    envmodules:
        config['tools']['cnvkit']['modname']
    shell: """
    cnvkit.py batch {input.tumor} -n {input.normal} \\
        -m amplicon \\
        --segment-method hmm --diagram \\
        --targets {params.targets} --fasta {params.genome} \\
        --annotate {params.refFlat} \\
        --output-dir {params.outdir}
    """

