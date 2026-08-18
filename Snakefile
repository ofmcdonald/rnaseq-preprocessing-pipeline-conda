import pandas as pd
import os

configfile: "config.yaml"

STAR_INDEX = "reference/star_index"

if not os.path.exists("samples.tsv"):
    raise FileNotFoundError(
        "samples.tsv not found. Run scripts/generate_samples_tsv.sh first."
    )

SAMPLES = pd.read_csv(
    "samples.tsv",
    sep="\t",
    dtype=str
)

if "sample" not in SAMPLES.columns:
    raise ValueError(
        "samples.tsv must contain a column named 'sample'"
    )

if not os.path.isdir("data_norm"):
    raise FileNotFoundError(
        "data_norm directory not found. Run scripts/build_data_norm.sh first."
    )

SAMPLE_IDS = SAMPLES["sample"].tolist()

if len(SAMPLE_IDS) == 0:
    raise ValueError(
        "No samples found in samples.tsv"
    )

dups = SAMPLES["sample"][SAMPLES["sample"].duplicated()]

if len(dups) > 0:
    raise ValueError(
        f"Duplicate sample IDs detected: {list(dups)}"
    )

# --------------------------------------------------
# IDENTIFY SINGLE-END VS PAIRED-END SAMPLES
# --------------------------------------------------

PE_SAMPLES = []
SE_SAMPLES = []

for sample in SAMPLE_IDS:

    r1 = f"data_norm/{sample}_R1.fastq.gz"
    r2 = f"data_norm/{sample}_R2.fastq.gz"
    se = f"data_norm/{sample}_SE.fastq.gz"

    has_r1 = os.path.exists(r1)
    has_r2 = os.path.exists(r2)
    has_se = os.path.exists(se)

    # Single-end
    if has_se:
        if has_r1 or has_r2:
            raise ValueError(
                f"Sample {sample} has both single-end and paired-end FASTQ files."
            )

        SE_SAMPLES.append(sample)

    # Paired-end
    elif has_r1 or has_r2:
        if not has_r1:
            raise FileNotFoundError(
                f"Missing R1 for paired-end sample: {sample}"
            )

        if not has_r2:
            raise FileNotFoundError(
                f"Missing R2 for paired-end sample: {sample}"
            )

        PE_SAMPLES.append(sample)

    # Nothing found
    else:
        raise FileNotFoundError(
            f"No FASTQ files found for sample: {sample}"
        )

print(f"Detected {len(PE_SAMPLES)} paired-end samples")
print(f"Detected {len(SE_SAMPLES)} single-end samples")

if PE_SAMPLES:
    print("Paired-end:", ", ".join(PE_SAMPLES))

if SE_SAMPLES:
    print("Single-end:", ", ".join(SE_SAMPLES))
    
def star_reads(wildcards):
    sample = wildcards.sample

    if sample in PE_SAMPLES:
        return [
            f"results/fastp/{sample}_R1_trimmed.fastq.gz",
            f"results/fastp/{sample}_R2_trimmed.fastq.gz"
        ]

    elif sample in SE_SAMPLES:
        return [
            f"results/fastp/{sample}_SE_trimmed.fastq.gz"
        ]

    else:
        raise ValueError(
            f"Sample {sample} not classified as PE or SE"
        )

# --------------------------------------------------
# FINAL TARGETS
# --------------------------------------------------

rule all:
    input:
        # FastQC raw - paired-end
        expand(
            "results/fastqc_raw/{sample}_R1_fastqc.zip",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastqc_raw/{sample}_R2_fastqc.zip",
            sample=PE_SAMPLES
        ),

        # FastQC raw - single-end
        expand(
            "results/fastqc_raw/{sample}_SE_fastqc.zip",
            sample=SE_SAMPLES
        ),

        # Trimmed FASTQ - paired-end
        expand(
            "results/fastp/{sample}_R1_trimmed.fastq.gz",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastp/{sample}_R2_trimmed.fastq.gz",
            sample=PE_SAMPLES
        ),

        # Trimmed FASTQ - single-end
        expand(
            "results/fastp/{sample}_SE_trimmed.fastq.gz",
            sample=SE_SAMPLES
        ),

        # FastQC trimmed - paired-end
        expand(
            "results/fastqc_fastp/{sample}_R1_trimmed_fastqc.zip",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastqc_fastp/{sample}_R2_trimmed_fastqc.zip",
            sample=PE_SAMPLES
        ),

        # FastQC trimmed - single-end
        expand(
            "results/fastqc_fastp/{sample}_SE_trimmed_fastqc.zip",
            sample=SE_SAMPLES
        ),

        # STAR outputs
        expand(
            "results/star/bam/{sample}.Aligned.sortedByCoord.out.bam",
            sample=SAMPLE_IDS
        ),
        expand(
            "results/samtools/{sample}.Aligned.sortedByCoord.out.bam.bai",
            sample=SAMPLE_IDS
        ),
        expand(
            "results/star/counts/{sample}.ReadsPerGene.out.tab",
            sample=SAMPLE_IDS
        ),

        # RSeQC
        expand(
            "results/rseqc/read_distribution/{sample}.read_distribution.txt",
            sample=SAMPLE_IDS
        ),
        expand(
            "results/rseqc/infer_experiment/{sample}.infer_experiment.txt",
            sample=SAMPLE_IDS
        ),

        # BED annotation
        "reference/annotation/Mus_musculus.GRCm39.bed",

        # MultiQC
        "results/multiqc/raw/multiqc_report.html",
        "results/multiqc/trimmed/multiqc_report.html"

# --------------------------------------------------
# STAR INDEX
# --------------------------------------------------

rule star_index:
    input:
        fasta=config["fasta"],
        gtf=config["gtf"]

    output:
        directory(STAR_INDEX)

    conda:
        "envs/star.yaml"

    threads: 8

    resources:
        mem_mb=64000,
        runtime=720

    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 149
        """

# --------------------------------------------------
# FASTQC RAW - PAIRED-END
# --------------------------------------------------

rule fastqc_raw_pe:
    input:
        r1="data_norm/{sample}_R1.fastq.gz",
        r2="data_norm/{sample}_R2.fastq.gz"

    output:
        zip1="results/fastqc_raw/{sample}_R1_fastqc.zip",
        zip2="results/fastqc_raw/{sample}_R2_fastqc.zip"

    log:
        "logs/fastqc_raw/{sample}.log"

    benchmark:
        "benchmarks/fastqc_raw/{sample}.txt"

    conda:
        "envs/fastqc.yaml"

    threads: 2

    resources:
        mem_mb=4000,
        runtime=60

    shell:
        """
        mkdir -p results/fastqc_raw
        mkdir -p logs/fastqc_raw

        fastqc \
            -t {threads} \
            -o results/fastqc_raw \
            {input.r1} \
            {input.r2} \
            &> {log}
        """

# --------------------------------------------------
# FASTQC RAW - SINGLE-END
# --------------------------------------------------

rule fastqc_raw_se:
    input:
        "data_norm/{sample}_SE.fastq.gz"

    output:
        "results/fastqc_raw/{sample}_SE_fastqc.zip"

    log:
        "logs/fastqc_raw/{sample}.log"

    benchmark:
        "benchmarks/fastqc_raw/{sample}.txt"

    conda:
        "envs/fastqc.yaml"

    threads: 2

    resources:
        mem_mb=4000,
        runtime=60

    shell:
        """
        mkdir -p results/fastqc_raw
        mkdir -p logs/fastqc_raw

        fastqc \
            -t {threads} \
            -o results/fastqc_raw \
            {input} \
            &> {log}
        """

# --------------------------------------------------
# FASTP - PAIRED-END
# --------------------------------------------------

rule fastp_pe:
    input:
        r1="data_norm/{sample}_R1.fastq.gz",
        r2="data_norm/{sample}_R2.fastq.gz"

    output:
        r1="results/fastp/{sample}_R1_trimmed.fastq.gz",
        r2="results/fastp/{sample}_R2_trimmed.fastq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json"

    log:
        "logs/fastp/{sample}.log"

    benchmark:
        "benchmarks/fastp/{sample}.txt"

    conda:
        "envs/fastp.yaml"

    threads: 2

    resources:
        mem_mb=8000,
        runtime=180

    shell:
        """
        mkdir -p results/fastp
        mkdir -p logs/fastp

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -w {threads} \
            -h {output.html} \
            -j {output.json} \
            &> {log}
        """

# --------------------------------------------------
# FASTP - SINGLE-END
# --------------------------------------------------

rule fastp_se:
    input:
        "data_norm/{sample}_SE.fastq.gz"

    output:
        trimmed="results/fastp/{sample}_SE_trimmed.fastq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json"

    log:
        "logs/fastp/{sample}.log"

    benchmark:
        "benchmarks/fastp/{sample}.txt"

    conda:
        "envs/fastp.yaml"

    threads: 2

    resources:
        mem_mb=8000,
        runtime=180

    shell:
        """
        mkdir -p results/fastp
        mkdir -p logs/fastp

        fastp \
            -i {input} \
            -o {output.trimmed} \
            -w {threads} \
            -h {output.html} \
            -j {output.json} \
            &> {log}
        """

# --------------------------------------------------
# FASTQC TRIMMED - PAIRED-END
# --------------------------------------------------

rule fastqc_trimmed_pe:
    input:
        r1="results/fastp/{sample}_R1_trimmed.fastq.gz",
        r2="results/fastp/{sample}_R2_trimmed.fastq.gz"

    output:
        zip1="results/fastqc_fastp/{sample}_R1_trimmed_fastqc.zip",
        zip2="results/fastqc_fastp/{sample}_R2_trimmed_fastqc.zip"

    log:
        "logs/fastqc_trimmed/{sample}.log"

    benchmark:
        "benchmarks/fastqc_trimmed/{sample}.txt"

    conda:
        "envs/fastqc.yaml"

    threads: 2

    resources:
        mem_mb=4000,
        runtime=60

    shell:
        """
        mkdir -p results/fastqc_fastp
        mkdir -p logs/fastqc_trimmed

        fastqc \
            -t {threads} \
            -o results/fastqc_fastp \
            {input.r1} \
            {input.r2} \
            &> {log}
        """

# --------------------------------------------------
# FASTQC TRIMMED - SINGLE-END
# --------------------------------------------------

rule fastqc_trimmed_se:
    input:
        "results/fastp/{sample}_SE_trimmed.fastq.gz"

    output:
        "results/fastqc_fastp/{sample}_SE_trimmed_fastqc.zip"

    log:
        "logs/fastqc_trimmed/{sample}.log"

    benchmark:
        "benchmarks/fastqc_trimmed/{sample}.txt"

    conda:
        "envs/fastqc.yaml"

    threads: 2

    resources:
        mem_mb=4000,
        runtime=60

    shell:
        """
        mkdir -p results/fastqc_fastp
        mkdir -p logs/fastqc_trimmed

        fastqc \
            -t {threads} \
            -o results/fastqc_fastp \
            {input} \
            &> {log}
        """

# --------------------------------------------------
# STAR ALIGNMENT
# --------------------------------------------------

rule star_align:
    input:
        reads=star_reads,
        index=STAR_INDEX

    output:
        bam="results/star/bam/{sample}.Aligned.sortedByCoord.out.bam",
        counts="results/star/counts/{sample}.ReadsPerGene.out.tab"

    log:
        "logs/star_align/{sample}.log"

    benchmark:
        "benchmarks/star_align/{sample}.txt"

    conda:
        "envs/star.yaml"

    threads: 8

    resources:
        mem_mb=48000,
        runtime=720

    shell:
        """
        mkdir -p results/star/bam
        mkdir -p results/star/counts
        mkdir -p results/star/tmp
        mkdir -p logs/star_align

        rm -rf results/star/tmp/{wildcards.sample}

        STAR \
            --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outTmpDir results/star/tmp/{wildcards.sample} \
            --outFileNamePrefix results/star/{wildcards.sample}. \
            &> {log}

        mv \
            results/star/{wildcards.sample}.Aligned.sortedByCoord.out.bam \
            {output.bam}

        mv \
            results/star/{wildcards.sample}.ReadsPerGene.out.tab \
            {output.counts}

        rm -rf results/star/tmp/{wildcards.sample}
        """

# --------------------------------------------------
# BAM INDEX
# --------------------------------------------------

rule samtools_index:
    input:
        bam="results/star/bam/{sample}.Aligned.sortedByCoord.out.bam"

    output:
        bai="results/samtools/{sample}.Aligned.sortedByCoord.out.bam.bai"

    log:
        "logs/samtools_index/{sample}.log"

    benchmark:
        "benchmarks/samtools_index/{sample}.txt"

    conda:
        "envs/samtools.yaml"

    threads: 2

    resources:
        mem_mb=4000,
        runtime=60

    shell:
        """
        mkdir -p results/samtools
        mkdir -p logs/samtools_index

        samtools index \
            -@ {threads} \
            -o {output.bai} \
            {input.bam} \
            &> {log}
        """

# --------------------------------------------------
# GTF -> BED12
# --------------------------------------------------

rule gtf_to_bed:
    input:
        gtf=config["gtf"]

    output:
        bed="reference/annotation/Mus_musculus.GRCm39.bed"

    conda:
        "envs/ucsc_tools.yaml"

    threads: 1

    resources:
        mem_mb=2000,
        runtime=30

    shell:
        """
        mkdir -p reference/annotation

        gtfToGenePred \
            {input.gtf} \
            reference/annotation/Mus_musculus.GRCm39.genepred

        genePredToBed \
            reference/annotation/Mus_musculus.GRCm39.genepred \
            {output.bed}

        rm reference/annotation/Mus_musculus.GRCm39.genepred
        """

# --------------------------------------------------
# RSeQC READ DISTRIBUTION
# --------------------------------------------------

rule rseqc_read_distribution:
    input:
        bam="results/star/bam/{sample}.Aligned.sortedByCoord.out.bam",
        bai="results/samtools/{sample}.Aligned.sortedByCoord.out.bam.bai",
        bed="reference/annotation/Mus_musculus.GRCm39.bed"

    output:
        "results/rseqc/read_distribution/{sample}.read_distribution.txt"

    log:
        "logs/rseqc/read_distribution/{sample}.log"

    benchmark:
        "benchmarks/rseqc/read_distribution/{sample}.txt"

    conda:
        "envs/rseqc.yaml"

    threads: 1

    resources:
        mem_mb=8000,
        runtime=120

    shell:
        """
        mkdir -p results/rseqc/read_distribution
        mkdir -p logs/rseqc/read_distribution

        read_distribution.py \
            -i {input.bam} \
            -r {input.bed} \
            > {output} 2> {log}
        """

# --------------------------------------------------
# RSeQC INFER EXPERIMENT (STRANDEDNESS)
# --------------------------------------------------

rule rseqc_infer_experiment:
    input:
        bam="results/star/bam/{sample}.Aligned.sortedByCoord.out.bam",
        bai="results/samtools/{sample}.Aligned.sortedByCoord.out.bam.bai",
        bed="reference/annotation/Mus_musculus.GRCm39.bed"

    output:
        "results/rseqc/infer_experiment/{sample}.infer_experiment.txt"

    log:
        "logs/rseqc/infer_experiment/{sample}.log"

    benchmark:
        "benchmarks/rseqc/infer_experiment/{sample}.txt"

    conda:
        "envs/rseqc.yaml"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=60

    shell:
        """
        mkdir -p results/rseqc/infer_experiment
        mkdir -p logs/rseqc/infer_experiment

        infer_experiment.py \
            -r {input.bed} \
            -i {input.bam} \
            > {output} 2> {log}
        """

# --------------------------------------------------
# MULTIQC RAW
# --------------------------------------------------

rule multiqc_raw:
    input:
        expand(
            "results/fastqc_raw/{sample}_R1_fastqc.zip",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastqc_raw/{sample}_R2_fastqc.zip",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastqc_raw/{sample}_SE_fastqc.zip",
            sample=SE_SAMPLES
        )

    output:
        "results/multiqc/raw/multiqc_report.html"

    log:
        "logs/multiqc/raw.log"

    benchmark:
        "benchmarks/multiqc/raw.txt"

    conda:
        "envs/multiqc.yaml"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    shell:
        """
        mkdir -p results/multiqc/raw
        mkdir -p logs/multiqc

        multiqc \
            results/fastqc_raw \
            --outdir results/multiqc/raw \
            --filename multiqc_report.html \
            --force \
            &> {log}
        """

# --------------------------------------------------
# MULTIQC TRIMMED
# --------------------------------------------------

rule multiqc_trimmed:
    input:
        expand(
            "results/fastqc_fastp/{sample}_R1_trimmed_fastqc.zip",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastqc_fastp/{sample}_R2_trimmed_fastqc.zip",
            sample=PE_SAMPLES
        ),
        expand(
            "results/fastqc_fastp/{sample}_SE_trimmed_fastqc.zip",
            sample=SE_SAMPLES
        )

    output:
        "results/multiqc/trimmed/multiqc_report.html"

    log:
        "logs/multiqc/trimmed.log"

    benchmark:
        "benchmarks/multiqc/trimmed.txt"

    conda:
        "envs/multiqc.yaml"

    threads: 1

    resources:
        mem_mb=4000,
        runtime=30

    shell:
        """
        mkdir -p results/multiqc/trimmed
        mkdir -p logs/multiqc

        multiqc \
            results/fastqc_fastp \
            --outdir results/multiqc/trimmed \
            --filename multiqc_report.html \
            --force \
            &> {log}
        """
