localrules: nidap
rule nidap:
    input:
        get_nidap_folder_input_files
    output: # these should exactly match the files in get_nidap_folder_output_files
        expand(os.path.join(NIDAP_OUTDIR,"{vc_outdir}_{samples}.artifact_filtered.vcf.gz"), samples=ffpe_sample_list, vc_outdir=ffpe_caller_list),
        expand(os.path.join(NIDAP_OUTDIR,"{vc_outdir}_all_somatic_variants.maf"), vc_outdir=ffpe_caller_list),
        expand(os.path.join(NIDAP_OUTDIR,"{vc_outdir}_all_metrics.txt"), vc_outdir=ffpe_caller_list),
        expand(os.path.join(NIDAP_OUTDIR,"{samples}.recal.bam_CNVs.p.value.txt"), samples=cnv_sample_list),
        expand(os.path.join(NIDAP_OUTDIR,"{samples}.contamination.table"), samples=samples_for_caller_merge),
        os.path.join(NIDAP_OUTDIR,"MultiQC_Report.html")
    params:
        outdir=os.path.join(NIDAP_OUTDIR)
    shell:"""
set -exo pipefail
if [ -d {params.outdir} ];then rm -rf {params.outdir};fi
mkdir -p {params.outdir}
cd {params.outdir}
# last file in inputs is NIDAP_files.tsv ... col1 is file ... col2 is the same file hardlinked in the NIDAP folder
# this file is created in get_nidap_folder_input_files function
linking_file=$(echo {input}|awk '{{print $NF}}')
while read a b;do
    ln $a $b
done < $linking_file
"""
