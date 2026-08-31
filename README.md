![version](https://img.shields.io/badge/version-2.0.0-blue)

# RNA-seq Preprocessing Pipeline - Conda (Mus musculus)

This repository contains a fully automated Snakemake workflow for preprocessing **single-end and paired-end RNA-seq data from *Mus musculus* (mouse)**.

The pipeline performs quality control, adapter trimming and quality filtering, genome alignment, gene-level quantification, BAM indexing, RNA-seq QC, and MultiQC reporting using reproducible conda-based environments and SLURM execution.

Version 2.0.0 adds support for single-end RNA-seq data while retaining paired-end functionality.

---

## Quick start

```bash
conda activate snakemake-pipeline

bash scripts/build_data_norm.sh
bash scripts/generate_samples_tsv.sh

bash scripts/run_snakemake.sh -n
bash scripts/run_snakemake.sh
````

---

## Overview

The pipeline performs the following steps:

1. **Input validation and sample detection**
2. **Raw FASTQ QC (FastQC)**
3. **Adapter trimming and quality filtering (fastp)**
4. **Post-trimming QC (FastQC)**
5. **Genome alignment (STAR)**
6. **Gene-level quantification (STAR GeneCounts)**
7. **BAM indexing (samtools)**
8. **Strandedness QC (RSeQC infer_experiment)**
9. **Read distribution QC (RSeQC read_distribution)**
10. **MultiQC aggregation (raw + trimmed)**

The workflow automatically determines whether each sample is single-end or paired-end based on the normalized FASTQ files.

Both sequencing layouts are supported:

* Single-end (SE)
* Paired-end (PE)

A dataset may contain either SE or PE samples. The workflow is also designed to classify samples individually, allowing SE and PE samples to coexist in the same project.

---

## Version 2.0.0

Version 2.0.0 adds single-end RNA-seq support while retaining the original paired-end workflow.

Major changes:

* Added single-end FASTQ recognition to `scripts/build_data_norm.sh`
* Updated `scripts/generate_samples_tsv.sh` to support single-end and paired-end samples
* Added automatic SE/PE sample classification to the Snakemake workflow
* Added single-end processing with FastQC, fastp, and STAR
* Retained paired-end processing with FastQC, fastp, and STAR
* Updated MultiQC targets to include both sequencing layouts
* Updated documentation for SE and PE input data
* Retained SLURM-based execution and conda environment isolation

The updated workflow has been tested with both single-end and paired-end RNA-seq data.

---

## Requirements

### 1. Conda environment system

This pipeline is designed to work with Miniforge3, but any Conda-compatible installation (Miniconda or Anaconda) is acceptable.

Conda must be available in your environment:

```bash
conda --version
which conda
```

If needed, add Miniforge to your PATH:

```bash
export PATH="$HOME/miniforge3/bin:$PATH"
```

To make this persistent, add it to your `.bashrc` or `.bash_profile`.

---

### 2. Snakemake execution environment

Create a dedicated environment for workflow execution:

```bash
conda create -n snakemake-pipeline \
    snakemake \
    snakemake-executor-plugin-slurm
```

Activate before running:

```bash
conda activate snakemake-pipeline
```

---

### 3. SLURM HPC configuration

This pipeline is designed for SLURM-based HPC systems and uses the Snakemake SLURM executor plugin.

Execution is handled via:

```text
scripts/run_snakemake.sh
```

The pipeline uses:

* Snakemake with the SLURM executor
* `profiles/slurm/` for cluster-specific configuration
* Per-rule resource definitions in the Snakefile
* Automatic SLURM job submission and orchestration

The Snakemake executable is configured in:

```bash
SNAKEMAKE_BIN="snakemake"
```

If the `snakemake` executable is not on your PATH, provide its full path:

```bash
SNAKEMAKE_BIN="$HOME/miniforge3/envs/snakemake-pipeline/bin/snakemake"
```

The appropriate setting is cluster-dependent.

---

## Input data requirements

### 1. Reference genome

The pipeline requires a local reference genome and annotation.

The default configuration expects:

```text
reference/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa
reference/annotation/Mus_musculus.GRCm39.115.gtf
```

These are not downloaded automatically and must be provided by the user.

Recommended sources include:

* Ensembl GRCm39 reference genome
* Corresponding Ensembl GRCm39 annotation

The FASTA and GTF paths are specified in:

```text
config.yaml
```

---

### 2. Raw FASTQ files

Place raw FASTQ files in the `data/` directory.

The pipeline supports both paired-end and single-end FASTQ data.

#### Paired-end input

Paired-end files may use common naming conventions such as:

```text
data/
├── SAMPLE1_R1.fastq.gz
├── SAMPLE1_R2.fastq.gz
├── SAMPLE2_R1.fastq.gz
└── SAMPLE2_R2.fastq.gz
```

or:

```text
data/
├── SAMPLE1_1.fastq.gz
├── SAMPLE1_2.fastq.gz
├── SAMPLE2_1.fastq.gz
└── SAMPLE2_2.fastq.gz
```

The normalization script also recognizes Illumina-style names such as:

```text
SAMPLE1_R1_001.fastq.gz
SAMPLE1_R2_001.fastq.gz
```

#### Single-end input

Single-end files can use simple sample-based names such as:

```text
data/
├── SAMPLE1.fastq.gz
├── SAMPLE2.fastq.gz
└── SAMPLE3.fastq.gz
```

The normalization script converts these to:

```text
data_norm/
├── SAMPLE1_SE.fastq.gz
├── SAMPLE2_SE.fastq.gz
└── SAMPLE3_SE.fastq.gz
```

#### FASTQ compression

Input FASTQ files must be gzip-compressed (`*.fastq.gz`).

If FASTQ files were generated using SRA Toolkit (`fasterq-dump`), compress them before running the pipeline:

```bash
gzip data/*.fastq
```

#### FASTQ normalization

The script:

```text
scripts/build_data_norm.sh
```

creates a normalized set of symbolic links in:

```text
data_norm/
```

Users should place raw FASTQ files in `data/` rather than manually creating files in `data_norm/`.

The normalization step converts supported input naming conventions into standardized names.

Paired-end:

```text
Sample1_R1.fastq.gz
Sample1_R2.fastq.gz
```

Single-end:

```text
Sample1_SE.fastq.gz
```

Run:

```bash
bash scripts/build_data_norm.sh
```

before generating the sample sheet.

---

### 3. Sample sheet

A sample sheet is required:

```text
samples.tsv
```

The file contains one sample ID per row:

| sample  |
| ------- |
| SAMPLE1 |
| SAMPLE2 |
| SAMPLE3 |

This file is generated automatically using:

```bash
bash scripts/generate_samples_tsv.sh
```

The script reads the normalized files in `data_norm/` and generates the sample IDs.

For paired-end samples, both R1 and R2 must be present.

For single-end samples, the normalized `_SE.fastq.gz` file must be present.

---

## Preprocessing steps required before Snakemake

Run the preprocessing scripts in order:

```bash
bash scripts/build_data_norm.sh
bash scripts/generate_samples_tsv.sh
```

These scripts:

1. Identify supported FASTQ naming conventions
2. Create normalized symbolic links in `data_norm/`
3. Identify sample IDs
4. Generate `samples.tsv`

Then run the Snakemake workflow:

```bash
bash scripts/run_snakemake.sh
```

---

## Automatic SE/PE detection

The Snakefile automatically classifies each sample by examining the files in `data_norm/`.

### Single-end sample

A sample is classified as single-end when:

```text
data_norm/SAMPLE_SE.fastq.gz
```

exists.

### Paired-end sample

A sample is classified as paired-end when both:

```text
data_norm/SAMPLE_R1.fastq.gz
data_norm/SAMPLE_R2.fastq.gz
```

exist.

The workflow performs validation to ensure that:

* A paired-end sample has both R1 and R2
* A sample does not contain both SE and PE files
* Every sample in `samples.tsv` has corresponding FASTQ data
* No duplicate sample IDs are present

The Snakefile reports the detected sample types at startup.

For example, for a single-end dataset:

```text
Detected 0 paired-end samples
Detected 12 single-end samples
Single-end: SAMPLE1, SAMPLE2, ...
```

For a paired-end dataset:

```text
Detected 12 paired-end samples
Detected 0 single-end samples
Paired-end: SAMPLE1, SAMPLE2, ...
```

For a mixed dataset:

```text
Detected 6 paired-end samples
Detected 6 single-end samples
Paired-end: SAMPLE1, SAMPLE2, ...
Single-end: SAMPLE7, SAMPLE8, ...
```

---

## Running the pipeline

Use the provided wrapper script:

```bash
bash scripts/run_snakemake.sh
```

This wrapper:

* Executes Snakemake from the repository root
* Uses the SLURM executor
* Loads the SLURM profile in `profiles/slurm/`
* Allows up to 100 concurrent jobs
* Passes additional Snakemake arguments through to the workflow

Additional Snakemake arguments can be supplied directly:

```bash
bash scripts/run_snakemake.sh -n
```

For example:

```bash
bash scripts/run_snakemake.sh --rerun-incomplete
```

or:

```bash
bash scripts/run_snakemake.sh --cores 1
```

Cluster-specific options should be configured through the SLURM profile rather than hard-coded into individual workflow rules.

---

## Dry run

A dry run is recommended before submitting the workflow:

```bash
bash scripts/run_snakemake.sh -n
```

The dry run should report the detected sample types and show the jobs that will be executed without submitting them to SLURM.

For single-end samples, the workflow should use:

* Single-end FastQC
* Single-end fastp
* Single-end STAR alignment

For paired-end samples, the workflow should use:

* Paired-end FastQC
* Paired-end fastp
* Paired-end STAR alignment

---

## Output structure

```text
results/
├── fastqc_raw/
├── fastp/
├── fastqc_fastp/
├── star/
│   ├── bam/
│   └── counts/
├── samtools/
├── rseqc/
│   ├── read_distribution/
│   └── infer_experiment/
└── multiqc/
    ├── raw/
    └── trimmed/
```

Key outputs include:

* Raw FastQC reports: `results/fastqc_raw/`
* Trimmed FASTQs: `results/fastp/`
* Trimmed FastQC reports: `results/fastqc_fastp/`
* STAR alignments: `results/star/bam/`
* STAR gene counts: `results/star/counts/`
* BAM indices: `results/samtools/`
* RSeQC reports: `results/rseqc/`
* MultiQC reports: `results/multiqc/`

---

## STAR reference index

The STAR genome index is generated automatically on the first run using the reference FASTA and GTF specified in `config.yaml`.

The generated index is stored in:

```text
reference/star_index/
```

This directory is large and computationally expensive to generate. It should not be deleted unless the STAR index needs to be rebuilt.

---

## Important design assumptions

* Supports single-end and paired-end RNA-seq
* Mouse reference genome (*Mus musculus*)
* FASTQ files must be gzip-compressed
* Input FASTQ files are placed in `data/`
* Normalized FASTQ files are generated in `data_norm/`
* Sample IDs must match the normalized FASTQ filenames
* Reference FASTA and GTF must be provided locally
* Workflow assumes a SLURM-based HPC system
* Software dependencies are isolated using conda environments
* STAR is used for genome alignment and gene-level quantification
* RSeQC is used for strandedness and read-distribution QC

---

## Troubleshooting

### Missing FASTQ errors

For a single-end sample, ensure:

```text
data_norm/SAMPLE_SE.fastq.gz
```

exists.

For a paired-end sample, ensure both:

```text
data_norm/SAMPLE_R1.fastq.gz
data_norm/SAMPLE_R2.fastq.gz
```

exist.

If files are missing, rerun:

```bash
bash scripts/build_data_norm.sh
bash scripts/generate_samples_tsv.sh
```

Then inspect:

```bash
ls -lh data_norm/
cat samples.tsv
```

---

### Sample detected as the wrong sequencing type

Inspect the normalized files:

```bash
ls -lh data_norm/
```

A single-end sample should have:

```text
SAMPLE_SE.fastq.gz
```

A paired-end sample should have:

```text
SAMPLE_R1.fastq.gz
SAMPLE_R2.fastq.gz
```

A sample cannot contain both SE and PE normalized files.

If the normalized files are incorrect, remove and rebuild the normalization directory:

```bash
bash scripts/build_data_norm.sh
bash scripts/generate_samples_tsv.sh
```

---

### `samples.tsv` contains no samples

Inspect the normalized directory:

```bash
ls -lh data_norm/
```

Then regenerate the sample sheet:

```bash
bash scripts/generate_samples_tsv.sh
```

Check the resulting file:

```bash
cat samples.tsv
```

It should contain:

```text
sample
SAMPLE1
SAMPLE2
SAMPLE3
```

---

### Snakemake executor issues

Verify that the SLURM executor plugin is installed:

```bash
conda list snakemake-executor-plugin-slurm
```

Also verify:

```bash
which snakemake
snakemake --version
```

If the Snakemake executable is not found, update `SNAKEMAKE_BIN` in:

```text
scripts/run_snakemake.sh
```

For example:

```bash
SNAKEMAKE_BIN="$HOME/miniforge3/envs/snakemake-pipeline/bin/snakemake"
```

---

### Environment issues

Ensure the Snakemake environment is active:

```bash
which conda
conda info
```

Then verify that Snakemake is available:

```bash
which snakemake
snakemake --version
```

---

### STAR index problems

The STAR index is generated from the FASTA and GTF specified in `config.yaml`.

Verify that both files exist:

```bash
ls -lh reference/genome/
ls -lh reference/annotation/
```

If the index needs to be rebuilt, remove:

```text
reference/star_index/
```

and rerun the workflow.

---

## Reproducibility notes

This pipeline is designed to be:

* Reproducible through conda-isolated software environments
* Scalable across HPC resources using SLURM
* Deterministic given identical inputs, reference files, software environments, and configuration
* Applicable to both single-end and paired-end RNA-seq datasets

The versions of individual bioinformatics tools are specified in the environment files under:

```text
envs/
```

The workflow configuration is specified in:

```text
config.yaml
```

Cluster-specific execution settings are specified under:

```text
profiles/slurm/
```

---

## Repository structure

The repository is organized approximately as follows:

```text
rnaseq-preprocessing-pipeline/
├── Snakefile
├── config.yaml
├── README.md
├── scripts/
│   ├── build_data_norm.sh
│   ├── generate_samples_tsv.sh
│   └── run_snakemake.sh
├── profiles/
│   └── slurm/
├── envs/
├── reference/
│   ├── genome/
│   ├── annotation/
│   └── star_index/
├── data/
├── data_norm/
├── results/
├── logs/
└── benchmarks/
```

Large input, output, and generated reference files should generally not be committed to GitHub.

---

## Author / maintenance

Olivia McDonald, PhD
[okfavor@outlook.com](mailto:okfavor@outlook.com)

---

## Citation

If you use this pipeline, please cite it as:

McDonald OF (2026). RNA-seq Preprocessing Pipeline - Conda (Mus musculus).
GitHub repository: [https://github.com/ofmcdonald/rnaseq-preprocessing-pipeline-conda](https://github.com/ofmcdonald/rnaseq-preprocessing-pipeline-conda)
Version: 2.0.0

---

## License

This software is released under the MIT License to facilitate reuse in academic and non-academic settings.
