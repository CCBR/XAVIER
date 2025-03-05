## XAVIER development version

## XAVIER 3.1.3

- Fix breaking bug in varscan somaticFilter (#131, @dnousome)

## XAVIER 3.1.2

- Fix minor bug in `assess_significance.R` script associated with rule `freec_exome_somatic`. (#120, @samarth8392)
- Fix bug in multiqc docker container, which caused an error when running xavier from ccbrpipeliner/7. (#123, @kelly-sovacool)
- Fix `cache` subcommand to correctly read the container images config file. (#124, @kelly-sovacool)

## XAVIER 3.1.1

- New contributing guide available on GitHub and the documentation website. (#114, @kelly-sovacool)
- New `xavier debug` subcommand to determine the base directory for debugging purposes. (#114, @kelly-sovacool)
- Upgraded `ccbr_wes_base` docker to v1.1.0 with updated GATK version to v4.6.0.0 (#116, @samarth8392)
- Upgrade multiqc container to use v1.15. (#117, @kelly-sovacool)
- Upgrade memory for rule "bwa_mem" to 100G (#118, @samarth8392)

## XAVIER 3.1.0

### new features

- Add `xavier gui` subcommand to launch the graphical user interface. (#99, @kelly-sovacool)
  - Previously, `xavier_gui` (with an underscore) was a command in the `ccbrpipeliner` module.
- Provide default exome targets for hg38 and mm10, which can be overridden by the optional `--targets` argument. (#102, @kelly-sovacool)
  - Previously, the `--targets` argument was required with no defaults.
- Add new human test dataset (#27, @samarth8392)

### bug fixes

- Fix bug in the GUI that always used the Agilent targets file by default, instead of picking based on the genome. (#108, @samarth8392)
- Fix bug in the driver script that caused the snakemake module not to be loaded on biowulf in some cases. (#111, @kelly-sovacool)
- Increased memory for rules: BWA mem, qualimap, kraken. gatk_contamination is no longer a localrule. (#89, @samarth8392)

### documentation updates

- You can now cite XAVIER with the DOI [10.5281/zenodo.12727315](https://doi.org/10.5281/zenodo.12727315). (#88, @kelly-sovacool)
- The docs website now has a dropdown menu to select which version to view. The latest release is shown by default. (#150, @kelly-sovacool)
- Other minor documentation improvements. (#92, @kelly-sovacool; #93, @samarth8392)

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
