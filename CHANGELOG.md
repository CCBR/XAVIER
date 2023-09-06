# Version changes:

## v3.0

- added [FRCE](https://ncifrederick.cancer.gov/staff/frce/welcome) support.
- adding `$triggeroptions` to _snakemake_ cli.
- using `SLURM_SUBMIT_HOST` as a secondary environmental variable to extract cluster name information on compute nodes.
- `func set_tmp()` added to repeat `tmpdir` assignments accross multiple rules.
- `vcf2maf` fix applied. 
