# ─── QA/QC rules ─────────────────────────────────────────────────────────────
# Steps:  FastQC (raw) → MultiQC → Fastp → FastQC (trimmed) → MultiQC
#
# Cluster modules (Compute Canada):
#   module load StdEnv/2020 fastqc/0.11.9
#   module load StdEnv/2020 fastp/0.23.1
#   Singularity image for MultiQC: multiqc-1.14.sif

QC_DIR    = config["out"]["qc"]
FASTP_DIR = config["out"]["fastp"]


# ─── FastQC — raw reads ───────────────────────────────────────────────────────
rule fastqc_raw:
    input:
        r1 = raw_r1,
        r2 = raw_r2,
    output:
        html_r1 = QC_DIR + "/raw/{sample}_1_fastqc.html",
        zip_r1  = QC_DIR + "/raw/{sample}_1_fastqc.zip",
        html_r2 = QC_DIR + "/raw/{sample}_2_fastqc.html",
        zip_r2  = QC_DIR + "/raw/{sample}_2_fastqc.zip",
        # Aliased outputs expected by rule all
        alias_r1 = QC_DIR + "/raw/{sample}_fastqc.html",
    params:
        outdir = QC_DIR + "/raw",
    conda:
        "../envs/fastqc.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 4000,
        runtime = 120,
    log:
        QC_DIR + "/logs/fastqc_raw/{sample}.log",
    shell:
        """
        fastqc {input.r1} {input.r2} --outdir {params.outdir} --threads {resources.cpus} \
            2> {log}
        # Create expected alias for rule all
        ln -sf $(basename {output.html_r1}) {output.alias_r1}
        """


# ─── MultiQC — raw reads ──────────────────────────────────────────────────────
rule multiqc_raw:
    input:
        expand(QC_DIR + "/raw/{sample}_1_fastqc.zip", sample=SAMPLES),
    output:
        html = QC_DIR + "/multiqc_raw/multiqc_report.html",
    params:
        indir  = QC_DIR + "/raw",
        outdir = QC_DIR + "/multiqc_raw",
        name   = "multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 60,
    log:
        QC_DIR + "/logs/multiqc_raw.log",
    shell:
        """
        multiqc {params.indir} \
            --outdir {params.outdir} \
            --filename {params.name} \
            --force \
            2> {log}
        """


# ─── Fastp — adapter trimming and quality filtering ───────────────────────────
rule fastp:
    input:
        r1 = raw_r1,
        r2 = raw_r2,
    output:
        r1       = FASTP_DIR + "/{sample}.R1.fq.gz",
        r2       = FASTP_DIR + "/{sample}.R2.fq.gz",
        unpaired = FASTP_DIR + "/{sample}_unpaired_combined.fq.gz",
        html     = FASTP_DIR + "/{sample}.html",
        json     = FASTP_DIR + "/{sample}.json",
    params:
        extra = config["params"]["fastp"].get("extra", ""),
    conda:
        "../envs/fastp.yaml"
    threads:
        config["resources"]["default"]["cpus"]
    resources:
        cpus    = config["resources"]["default"]["cpus"],
        mem_mb  = 20000,
        runtime = 240,
    log:
        FASTP_DIR + "/logs/{sample}.log",
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            --unpaired1 {output.unpaired} \
            --unpaired2 {output.unpaired} \
            -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} \
            --thread {threads} \
            {params.extra} \
            2> {log}
        """


# ─── FastQC — trimmed reads ───────────────────────────────────────────────────
rule fastqc_trimmed:
    input:
        r1 = FASTP_DIR + "/{sample}.R1.fq.gz",
        r2 = FASTP_DIR + "/{sample}.R2.fq.gz",
    output:
        html_r1 = QC_DIR + "/trimmed/{sample}.R1_fastqc.html",
        zip_r1  = QC_DIR + "/trimmed/{sample}.R1_fastqc.zip",
        html_r2 = QC_DIR + "/trimmed/{sample}.R2_fastqc.html",
        zip_r2  = QC_DIR + "/trimmed/{sample}.R2_fastqc.zip",
    params:
        outdir = QC_DIR + "/trimmed",
    conda:
        "../envs/fastqc.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 4000,
        runtime = 120,
    log:
        QC_DIR + "/logs/fastqc_trimmed/{sample}.log",
    shell:
        """
        fastqc {input.r1} {input.r2} --outdir {params.outdir} --threads {resources.cpus} \
            2> {log}
        """


# ─── MultiQC — trimmed reads ──────────────────────────────────────────────────
rule multiqc_trimmed:
    input:
        expand(QC_DIR + "/trimmed/{sample}.R1_fastqc.zip", sample=SAMPLES),
    output:
        html = QC_DIR + "/multiqc_trimmed/multiqc_report.html",
    params:
        indir  = QC_DIR + "/trimmed",
        outdir = QC_DIR + "/multiqc_trimmed",
        name   = "multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 60,
    log:
        QC_DIR + "/logs/multiqc_trimmed.log",
    shell:
        """
        multiqc {params.indir} \
            --outdir {params.outdir} \
            --filename {params.name} \
            --force \
            2> {log}
        """
