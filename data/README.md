# Input FASTQ Files

Place raw FASTQ files in this directory.

The pipeline supports both **single-end (SE)** and **paired-end (PE)** RNA-seq data.

## Paired-end FASTQ files

Paired-end files can use common naming conventions such as:

```text
Sample1_R1.fastq.gz
Sample1_R2.fastq.gz
Sample2_R1.fastq.gz
Sample2_R2.fastq.gz
```

The pipeline also recognizes:

```text
Sample1_1.fastq.gz
Sample1_2.fastq.gz
Sample2_1.fastq.gz
Sample2_2.fastq.gz
```

and Illumina-style filenames:

```text
Sample1_R1_001.fastq.gz
Sample1_R2_001.fastq.gz
```

## Single-end FASTQ files

Single-end files should contain the sample name without a read-pair suffix:

```text
Sample1.fastq.gz
Sample2.fastq.gz
Sample3.fastq.gz
```

The `build_data_norm.sh` script automatically identifies these as single-end files and creates normalized symbolic links such as:

```text
data_norm/Sample1_SE.fastq.gz
data_norm/Sample2_SE.fastq.gz
```

## FASTQ requirements

* FASTQ files must be gzip-compressed (`*.fastq.gz`).
* Raw FASTQ files should be placed directly in this `data/` directory.
* Do not manually place files in `data_norm/`; that directory is generated automatically.
* Samples should have unique sample IDs.
* A paired-end sample must have both R1 and R2 files.
* A sample must not contain both single-end and paired-end FASTQ files.

## After adding FASTQ files

From the repository root, run:

```bash
bash scripts/build_data_norm.sh
bash scripts/generate_samples_tsv.sh
```

These scripts normalize the FASTQ filenames and generate `samples.tsv`.

Verify the normalized files:

```bash
ls -lh data_norm/
```

Then perform a Snakemake dry run before launching the full workflow:

```bash
bash scripts/run_snakemake.sh -n
```

The Snakefile will automatically determine whether each sample is single-end or paired-end and execute the appropriate workflow rules.

```
