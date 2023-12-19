# Version changes

## development version

## v3.0.2
 - Documentation updates for GUI and other genome references 
 - Additional IUPAC changes to play nice with GATK tools (non ACGTN codes convert to N)
   
## v3.0.1

Hotfixes for all!
 - Increased memory for somatic merge rule
 - Fixed spooker/runner calls for run info
 - Fixed vcf2maf which missed the Intersection of calls

## v3.0.0

- added [FRCE](https://ncifrederick.cancer.gov/staff/frce/welcome) support.
- adding `$triggeroptions` to _snakemake_ cli.
- using `SLURM_SUBMIT_HOST` as a secondary environmental variable to extract cluster name information on compute nodes.
- `func set_tmp()` added to repeat `tmpdir` assignments accross multiple rules.
- `vcf2maf` fix applied. 
