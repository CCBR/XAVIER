# CCBR XAVIER Pipeline
***XAVIER: eXome Analysis and Variant explorER***


**Before**:
```
CCBR_GATK4_Exome_Seq_Pipeline/
├── config
│   ├── cluster.json
│   ├── config.json
│   └── pairs.mixed.tsv
├── pipeline
│   └── run.sh
├── README.md
├── resources
│   └── fastq_screen.conf
└── workflow
    ├── rules
    │   ├── cnv.smk
    │   ├── ffpe.smk
    │   ├── germline.smk
    │   ├── qc.smk
    │   ├── somatic_snps.common.smk
    │   ├── somatic_snps.paired.smk
    │   ├── somatic_snps.tumor_only.smk
    │   └── trim_map_preprocess.smk
    ├── scripts
    │   ├── freec
    │   │   ├── make_freec_config.py
    │   │   └── py_config
    │   │       └── config.txt
    │   ├── get_flowcell_lanes.py
    │   ├── parse_tn_mode.py
    │   ├── reformat_bed.py
    │   └── RScripts
    │       ├── combineAllSampleCompareResults.R
    │       ├── combineVerifyBAMIDResults.R
    │       ├── predictGender.R
    │       └── sampleCompareAncestoryPlots.R
    └── Snakefile

9 directories, 24 files
```

**After**:
```

```

