## development version

- You can now cite XAVIER with the DOI [10.5281/zenodo.12727315](https://doi.org/10.5281/zenodo.12727315). (#88, @kelly-sovacool)
- Minor documentation improvements. (#92, @kelly-sovacool)
- Minor documentation rendering improvements (#93, @samarth8392)
- The docs website now has a dropdown menu to select which version to view. The latest release is shown by default. (#150, @kelly-sovacool)

## XAVIER 3.0.3

- Bug fixes:
  - Use `--rerun--incomplete` Snakemake argument and fix bed. (#66, @dnousome)
  - Fix final somatic merge (#70, @dnousome)
  - Fix Varscan/Mutect2 output to keep Tumor Normal ordering. (#84, @dnousome)
  - Reduce memory allocation. (#86, @samarth8392)
  - Provide a more helpful error message when `xavier` is called with no arguments. (#75, @kelly-sovacool)
- Documentation improvements: (#78, @kelly-sovacool)
  - Document the release process for developers. (#63, @kelly-sovacool)
  - Create `CITATION.cff` to describe how to cite XAVIER. (#68, @kelly-sovacool)
- Sync GitHub Actions with other CCBR pipelines. (#65, #67, @kelly-sovacool; #73, @kopardev)

## XAVIER 3.0.2

- Documentation updates for GUI and other genome references
- Additional IUPAC changes to play nice with GATK tools (non ACGTN codes convert to N)

## XAVIER 3.0.1

Hotfixes for all!

- Increased memory for somatic merge rule
- Fixed spooker/runner calls for run info
- Fixed vcf2maf which missed the Intersection of calls

## XAVIER 3.0.0

- added [FRCE](https://ncifrederick.cancer.gov/staff/frce/welcome) support.
- adding `$triggeroptions` to _snakemake_ cli.
- using `SLURM_SUBMIT_HOST` as a secondary environmental variable to extract cluster name information on compute nodes.
- `func set_tmp()` added to repeat `tmpdir` assignments across multiple rules.
- `vcf2maf` fix applied.
