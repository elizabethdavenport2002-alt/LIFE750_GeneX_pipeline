# LIFE750_GeneX_pipeline

# LIFE750 Cycle 2 — Command Record

# Contents
- `setup_conda_environment.txt` — conda installation and environment setup
- `track1_variant_calling.txt` — BWA alignment and freebayes variant calling pipeline
- `track2_differential.R` — DESeq2 analysis in RStudio
- `track3_binding_site_integration.txt` — bedtools intersect pipeline

# Expected File Structure 
```
assignment2/
├── data/
│   ├── GeneX_reference.fa
│   ├── GeneX.gff3
│   ├── genes.gff3
│   ├── genome.fa
│   ├── gene_counts.tsv
│   ├── sample_metadata.tsv
│   ├── gene_x_binding_sites.bed
│   ├── normal_GeneX_R1.fastq
│   ├── normal_GeneX_R2.fastq
│   ├── mutant_GeneX_R1.fastq
│   └── mutant_GeneX_R2.fastq
├── track1_results_variant/
├── track2_results/
└── track3_results/
```

## Software Used
- BWA v0.7.19
- SAMtools/BCFtools (Danecek et al., 2021)
- Freebayes (Garrison & Marth, 2012)
- BEDTools (Quinlan & Hall, 2010)
- R/DESeq2 (Love et al., 2014)
- IGV (Robinson et al., 2011)
