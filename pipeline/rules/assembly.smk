# ─── Assembly rules ───────────────────────────────────────────────────────────
# MEGAHIT co-assembly: one assembly per co-assembly group (site × type).
#
# Cluster modules (Compute Canada):
#   module load StdEnv/2020 megahit/1.2.9
#
# Note: individual-sample assemblies (used for exploration) are not wired into
# the binning workflow. Add rule megahit_individual (commented below) and
# request its output in rule all if needed.

COASSEMBLY_DIR = config["out"]["coassembly"]
FASTP_DIR      = config["out"]["fastp"]


# ─── MEGAHIT co-assembly ──────────────────────────────────────────────────────
rule megahit_coassembly:
    input:
        r1 = coassembly_r1,   # list of R1 files for all samples in the group
        r2 = coassembly_r2,
    output:
        # MEGAHIT writes directly into outdir; we expose the final contigs file
        contigs = COASSEMBLY_DIR + "/{site}/final.contigs.fa",
    params:
        outdir          = COASSEMBLY_DIR + "/{site}",
        min_contig_len  = config["params"]["megahit"]["min_contig_len"],
        presets         = config["params"]["megahit"]["presets"],
        extra           = config["params"]["megahit"]["extra"],
    conda:
        "../envs/megahit.yaml"
    threads:
        config["resources"]["assembly"]["cpus"]
    resources:
        cpus    = config["resources"]["assembly"]["cpus"],
        mem_mb  = config["resources"]["assembly"]["mem_mb"],
        runtime = config["resources"]["assembly"]["runtime"],
    log:
        COASSEMBLY_DIR + "/logs/{site}_megahit.log",
    shell:
        """
        # MEGAHIT requires the output directory to not already exist
        rm -rf {params.outdir}

        megahit \
            -1 $(IFS=,; echo "{input.r1}") \
            -2 $(IFS=,; echo "{input.r2}") \
            -o {params.outdir} \
            --presets {params.presets} \
            --min-contig-len {params.min_contig_len} \
            {params.extra} \
            -t {threads} \
            2> {log}
        """


# ─── (Optional) MEGAHIT individual assembly ───────────────────────────────────
# Uncomment and add to rule all if per-sample assemblies are needed.
#
# rule megahit_individual:
#     input:
#         r1 = FASTP_DIR + "/{sample}.R1.fq.gz",
#         r2 = FASTP_DIR + "/{sample}.R2.fq.gz",
#     output:
#         contigs = config["out"]["individual_assembly"] + "/{sample}/final.contigs.fa",
#     params:
#         outdir         = config["out"]["individual_assembly"] + "/{sample}",
#         min_contig_len = config["params"]["megahit"]["min_contig_len"],
#         # meta-sensitive is recommended for individual assemblies
#         presets        = "meta-sensitive",
#         extra          = "--no-mercy",
#     conda:
#         "../envs/megahit.yaml"
#     threads:
#         config["resources"]["assembly"]["cpus"]
#     resources:
#         cpus    = config["resources"]["assembly"]["cpus"],
#         mem_mb  = 125000,
#         runtime = 600,
#     log:
#         config["out"]["individual_assembly"] + "/logs/{sample}_megahit.log",
#     shell:
#         """
#         rm -rf {params.outdir}
#         megahit \
#             -1 {input.r1} -2 {input.r2} \
#             -o {params.outdir} \
#             --presets {params.presets} \
#             --min-contig-len {params.min_contig_len} \
#             {params.extra} \
#             -t {threads} \
#             2> {log}
#         """
