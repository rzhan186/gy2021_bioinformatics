# ─── Binning rules ────────────────────────────────────────────────────────────
# Pipeline:
#   1. filter_contigs      — anvi-script-reformat-fasta (min-len + simplify names)
#   2. bwa_index           — index filtered contigs
#   3. bwa_map             — map each sample back to its co-assembly
#   4. samtools_sort_index — sort + index BAM
#   5. metabat2_depth      — jgi_summarize_bam_contig_depths per site
#   6. metabat2            — MetaBAT2 binning
#   7. maxbin2             — MaxBin2 binning
#   8. vamb_bam_list       — build VAMB input TSV listing BAM paths
#   9. vamb                — VAMB binning (GPU)
#  10. concoct_cutup       — cut_up_fasta.py (CONCOCT preprocessing)
#  11. concoct_coverage    — concoct_coverage_table.py
#  12. concoct_run         — CONCOCT
#  13. concoct_merge       — merge_cutup_clustering.py
#  14. concoct_extract     — extract_fasta_bins.py
#  15. contig2bin_*        — Fasta_to_Contig2Bin.sh for each binner
#  16. dastool             — DAS_Tool refinement
#
# Cluster modules (Compute Canada):
#   Contig filtering : source anvio-7.1 environment
#   BWA + samtools   : module load bwa/0.7.17 samtools/1.17
#   MetaBAT2         : module load metabat/2.14
#   MaxBin2          : module load maxbin2/2.2.7 perl bowtie2
#   VAMB             : source ~/vamb4.1.1-ENV/bin/activate
#   CONCOCT          : apptainer run concoct-v1.1.0.sif
#   DAS_Tool         : module load r/4.0.2 prodigal diamond blast+/2.11.0 pullseq

import os

COASSEMBLY_DIR = config["out"]["coassembly"]
MAPPING_DIR    = config["out"]["mapping"]
BINS_DIR       = config["out"]["bins"]
FASTP_DIR      = config["out"]["fastp"]

BINNERS = ["metabat2", "maxbin2", "vamb", "concoct"]


# ─── 1. Filter co-assembly contigs ───────────────────────────────────────────
rule filter_contigs:
    input:
        fasta = COASSEMBLY_DIR + "/{site}/final.contigs.fa",
    output:
        fasta = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
    params:
        min_len = config["params"]["filter_contigs"]["min_len"],
    conda:
        "../envs/anvio.yaml"
    resources:
        cpus    = 4,
        mem_mb  = 16000,
        runtime = 60,
    log:
        COASSEMBLY_DIR + "/logs/filter_contigs/{site}.log",
    shell:
        """
        anvi-script-reformat-fasta {input.fasta} \
            --output-file {output.fasta} \
            --min-len {params.min_len} \
            --simplify-names \
            2> {log}
        """


# ─── 2. BWA index filtered contigs ───────────────────────────────────────────
rule bwa_index:
    input:
        fasta = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
    output:
        bwt   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa.bwt",
        pac   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa.pac",
        ann   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa.ann",
        amb   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa.amb",
        sa    = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa.sa",
    conda:
        "../envs/bwa_samtools.yaml"
    resources:
        cpus    = 4,
        mem_mb  = 16000,
        runtime = 60,
    log:
        COASSEMBLY_DIR + "/logs/bwa_index/{site}.log",
    shell:
        """
        bwa index -a bwtsw {input.fasta} 2> {log}
        """


# ─── 3 + 4. BWA mapping → sorted + indexed BAM ───────────────────────────────
# One job per (site, sample) pair — sample must belong to site's co-assembly group.
rule bwa_map:
    input:
        fasta = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
        bwt   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa.bwt",
        r1    = FASTP_DIR + "/{sample}.R1.fq.gz",
        r2    = FASTP_DIR + "/{sample}.R2.fq.gz",
    output:
        bam   = MAPPING_DIR + "/{site}/{sample}_sorted.bam",
        bai   = MAPPING_DIR + "/{site}/{sample}_sorted.bam.bai",
    params:
        # Discard unmapped reads (-F 0x4), output BAM (-bS)
        view_flags = config["params"]["bwa_mem"]["extra"],
        tmp_bam    = MAPPING_DIR + "/{site}/{sample}.unsorted.bam",
    conda:
        "../envs/bwa_samtools.yaml"
    threads:
        config["resources"]["mapping"]["cpus"]
    resources:
        cpus    = config["resources"]["mapping"]["cpus"],
        mem_mb  = config["resources"]["mapping"]["mem_mb"],
        runtime = config["resources"]["mapping"]["runtime"],
    log:
        MAPPING_DIR + "/logs/{site}/{sample}.log",
    shell:
        """
        bwa mem -t {threads} {input.fasta} {input.r1} {input.r2} 2>> {log} \
            | samtools view {params.view_flags} -bS --threads {threads} \
            > {params.tmp_bam}

        samtools sort -@ {threads} {params.tmp_bam} -o {output.bam} 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        rm -f {params.tmp_bam}
        """


# ─── 5. MetaBAT2 — depth calculation ─────────────────────────────────────────
rule metabat2_depth:
    input:
        bams    = site_sorted_bams,
        indices = site_bam_indices,
    output:
        depth = BINS_DIR + "/metabat2/{site}/depth.txt",
    conda:
        "../envs/metabat2.yaml"
    resources:
        cpus    = config["resources"]["binning"]["cpus"],
        mem_mb  = config["resources"]["binning"]["mem_mb"],
        runtime = 120,
    log:
        BINS_DIR + "/logs/metabat2_depth/{site}.log",
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bams} \
            2> {log}
        """


# ─── 6. MetaBAT2 binning ──────────────────────────────────────────────────────
rule metabat2:
    input:
        fasta = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
        depth = BINS_DIR + "/metabat2/{site}/depth.txt",
    output:
        # Sentinel directory; actual bins are {site}_bin_*.fa
        done  = BINS_DIR + "/metabat2/{site}/.done",
    params:
        outprefix      = BINS_DIR + "/metabat2/{site}/metabat2_{site}_bin_",
        min_contig_len = config["params"]["metabat2"]["min_contig_len"],
    conda:
        "../envs/metabat2.yaml"
    threads:
        config["resources"]["binning"]["cpus"]
    resources:
        cpus    = config["resources"]["binning"]["cpus"],
        mem_mb  = config["resources"]["binning"]["mem_mb"],
        runtime = config["resources"]["binning"]["runtime"],
    log:
        BINS_DIR + "/logs/metabat2/{site}.log",
    shell:
        """
        metabat2 \
            -i {input.fasta} \
            -a {input.depth} \
            -o {params.outprefix} \
            -m {params.min_contig_len} \
            -t {threads} \
            2> {log}
        touch {output.done}
        """


# ─── 7. MaxBin2 binning ───────────────────────────────────────────────────────
rule maxbin2:
    input:
        fasta = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
        depth = BINS_DIR + "/metabat2/{site}/depth.txt",   # reuse MetaBAT2 depth
    output:
        done  = BINS_DIR + "/maxbin2/{site}/.done",
    params:
        outprefix         = BINS_DIR + "/maxbin2/{site}/maxbin2_{site}_bin",
        min_contig_length = config["params"]["maxbin2"]["min_contig_length"],
    conda:
        "../envs/maxbin2.yaml"
    threads:
        config["resources"]["binning"]["cpus"]
    resources:
        cpus    = config["resources"]["binning"]["cpus"],
        mem_mb  = config["resources"]["binning"]["mem_mb"],
        runtime = config["resources"]["binning"]["runtime"],
    log:
        BINS_DIR + "/logs/maxbin2/{site}.log",
    shell:
        """
        run_MaxBin.pl \
            -contig {input.fasta} \
            -abund  {input.depth} \
            -out    {params.outprefix} \
            -min_contig_length {params.min_contig_length} \
            -thread {threads} \
            2> {log}
        touch {output.done}
        """


# ─── 8. VAMB — build BAM list TSV ────────────────────────────────────────────
rule vamb_bam_list:
    input:
        bams    = site_sorted_bams,
        indices = site_bam_indices,
    output:
        tsv = BINS_DIR + "/vamb/{site}/bam_list.tsv",
    run:
        with open(output.tsv, "w") as f:
            for bam in input.bams:
                sample = os.path.basename(bam).replace("_sorted.bam", "")
                f.write(f"{sample}\t{bam}\n")


# ─── 9. VAMB binning (GPU) ────────────────────────────────────────────────────
rule vamb:
    input:
        fasta   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
        bam_tsv = BINS_DIR + "/vamb/{site}/bam_list.tsv",
    output:
        bins_dir = directory(BINS_DIR + "/vamb/{site}/bins"),
        log_file = BINS_DIR + "/vamb/{site}/log.txt",
    params:
        outdir         = BINS_DIR + "/vamb/{site}",
        minfasta        = config["params"]["vamb"]["minfasta"],
        min_contig_len  = config["params"]["vamb"]["min_contig_len"],
        cuda_flag       = "--cuda" if config["params"]["vamb"]["use_cuda"] else "",
    conda:
        "../envs/vamb.yaml"
    threads:
        config["resources"]["vamb"]["cpus"]
    resources:
        cpus    = config["resources"]["vamb"]["cpus"],
        mem_mb  = config["resources"]["vamb"]["mem_mb"],
        runtime = config["resources"]["vamb"]["runtime"],
        gpu     = config["resources"]["vamb"]["gpu"],
    log:
        BINS_DIR + "/logs/vamb/{site}.log",
    shell:
        """
        # Remove output dir so VAMB can create it fresh
        rm -rf {params.outdir}

        vamb \
            --outdir {params.outdir} \
            --fasta  {input.fasta} \
            --bamfiles {input.bam_tsv} \
            --minfasta {params.minfasta} \
            -o C \
            -m {params.min_contig_len} \
            {params.cuda_flag} \
            -p {threads} \
            2> {log}
        """


# ─── 10. CONCOCT — cut up contigs ─────────────────────────────────────────────
rule concoct_cutup:
    input:
        fasta = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
    output:
        fa  = BINS_DIR + "/concoct/{site}/contigs_2000.fa",
        bed = BINS_DIR + "/concoct/{site}/contigs_2000.bed",
    params:
        chunk_size = config["params"]["concoct"]["chunk_size"],
        overlap    = config["params"]["concoct"]["overlap"],
    conda:
        "../envs/concoct.yaml"
    resources:
        cpus    = 4,
        mem_mb  = 16000,
        runtime = 60,
    log:
        BINS_DIR + "/logs/concoct/cutup_{site}.log",
    shell:
        """
        cut_up_fasta.py {input.fasta} \
            -c {params.chunk_size} \
            -o {params.overlap} \
            --merge_last \
            -b {output.bed} \
            > {output.fa} \
            2> {log}
        """


# ─── 11. CONCOCT — coverage table ────────────────────────────────────────────
rule concoct_coverage:
    input:
        bed     = BINS_DIR + "/concoct/{site}/contigs_2000.bed",
        bams    = site_sorted_bams,
        indices = site_bam_indices,
    output:
        tsv = BINS_DIR + "/concoct/{site}/coverage_table.tsv",
    conda:
        "../envs/concoct.yaml"
    resources:
        cpus    = config["resources"]["binning"]["cpus"],
        mem_mb  = config["resources"]["binning"]["mem_mb"],
        runtime = 120,
    log:
        BINS_DIR + "/logs/concoct/coverage_{site}.log",
    shell:
        """
        concoct_coverage_table.py {input.bed} {input.bams} \
            > {output.tsv} \
            2> {log}
        """


# ─── 12. CONCOCT — run ───────────────────────────────────────────────────────
rule concoct_run:
    input:
        fa       = BINS_DIR + "/concoct/{site}/contigs_2000.fa",
        coverage = BINS_DIR + "/concoct/{site}/coverage_table.tsv",
    output:
        clustering = BINS_DIR + "/concoct/{site}/clustering_gt1000.csv",
    params:
        outdir = BINS_DIR + "/concoct/{site}",
    conda:
        "../envs/concoct.yaml"
    threads:
        config["resources"]["binning"]["cpus"]
    resources:
        cpus    = config["resources"]["binning"]["cpus"],
        mem_mb  = config["resources"]["binning"]["mem_mb"],
        runtime = config["resources"]["binning"]["runtime"],
    log:
        BINS_DIR + "/logs/concoct/run_{site}.log",
    shell:
        """
        concoct \
            --composition_file {input.fa} \
            --coverage_file {input.coverage} \
            -b {params.outdir} \
            -t {threads} \
            2> {log}
        """


# ─── 13. CONCOCT — merge cut-up clustering back to original contigs ───────────
rule concoct_merge:
    input:
        clustering = BINS_DIR + "/concoct/{site}/clustering_gt1000.csv",
    output:
        merged = BINS_DIR + "/concoct/{site}/clustering_merged.csv",
    conda:
        "../envs/concoct.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 30,
    log:
        BINS_DIR + "/logs/concoct/merge_{site}.log",
    shell:
        """
        merge_cutup_clustering.py {input.clustering} > {output.merged} 2> {log}
        """


# ─── 14. CONCOCT — extract bin FASTA files ───────────────────────────────────
rule concoct_extract:
    input:
        fasta   = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
        merged  = BINS_DIR + "/concoct/{site}/clustering_merged.csv",
    output:
        bins_dir = directory(BINS_DIR + "/concoct/{site}/bins"),
    conda:
        "../envs/concoct.yaml"
    resources:
        cpus    = 4,
        mem_mb  = 16000,
        runtime = 60,
    log:
        BINS_DIR + "/logs/concoct/extract_{site}.log",
    shell:
        """
        mkdir -p {output.bins_dir}
        extract_fasta_bins.py {input.fasta} {input.merged} \
            --output_path {output.bins_dir} \
            2> {log}
        """


# ─── 15. Fasta_to_Contig2Bin — one rule per binner ───────────────────────────
# DAS_Tool requires a contig → bin TSV for each binner.

rule contig2bin_metabat2:
    input:
        done = BINS_DIR + "/metabat2/{site}/.done",
    output:
        tsv  = BINS_DIR + "/dastool/{site}/{site}_metabat2.contig2bin.tsv",
    params:
        bins_dir = BINS_DIR + "/metabat2/{site}",
        ext      = "fa",
    conda:
        "../envs/dastool.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 30,
    log:
        BINS_DIR + "/logs/contig2bin/metabat2_{site}.log",
    shell:
        """
        Fasta_to_Contig2Bin.sh -i {params.bins_dir} -e {params.ext} \
            > {output.tsv} 2> {log}
        """


rule contig2bin_maxbin2:
    input:
        done = BINS_DIR + "/maxbin2/{site}/.done",
    output:
        tsv  = BINS_DIR + "/dastool/{site}/{site}_maxbin2.contig2bin.tsv",
    params:
        bins_dir = BINS_DIR + "/maxbin2/{site}",
        ext      = "fasta",
    conda:
        "../envs/dastool.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 30,
    log:
        BINS_DIR + "/logs/contig2bin/maxbin2_{site}.log",
    shell:
        """
        Fasta_to_Contig2Bin.sh -i {params.bins_dir} -e {params.ext} \
            > {output.tsv} 2> {log}
        """


rule contig2bin_vamb:
    input:
        bins_dir = BINS_DIR + "/vamb/{site}/bins",
    output:
        tsv      = BINS_DIR + "/dastool/{site}/{site}_vamb.contig2bin.tsv",
    params:
        ext = "fna",
    conda:
        "../envs/dastool.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 30,
    log:
        BINS_DIR + "/logs/contig2bin/vamb_{site}.log",
    shell:
        """
        Fasta_to_Contig2Bin.sh -i {input.bins_dir} -e {params.ext} \
            > {output.tsv} 2> {log}
        """


rule contig2bin_concoct:
    input:
        bins_dir = BINS_DIR + "/concoct/{site}/bins",
    output:
        tsv      = BINS_DIR + "/dastool/{site}/{site}_concoct.contig2bin.tsv",
    params:
        ext = "fa",
    conda:
        "../envs/dastool.yaml"
    resources:
        cpus    = 2,
        mem_mb  = 8000,
        runtime = 30,
    log:
        BINS_DIR + "/logs/contig2bin/concoct_{site}.log",
    shell:
        """
        Fasta_to_Contig2Bin.sh -i {input.bins_dir} -e {params.ext} \
            > {output.tsv} 2> {log}
        """


# ─── 16. DAS_Tool — bin refinement ───────────────────────────────────────────
rule dastool:
    input:
        metabat2 = BINS_DIR + "/dastool/{site}/{site}_metabat2.contig2bin.tsv",
        maxbin2  = BINS_DIR + "/dastool/{site}/{site}_maxbin2.contig2bin.tsv",
        vamb     = BINS_DIR + "/dastool/{site}/{site}_vamb.contig2bin.tsv",
        concoct  = BINS_DIR + "/dastool/{site}/{site}_concoct.contig2bin.tsv",
        fasta    = COASSEMBLY_DIR + "/{site}/contigs.filtered.fa",
    output:
        summary  = BINS_DIR + "/dastool/{site}/dastool_DASTool_summary.txt",
        bins_dir = directory(BINS_DIR + "/dastool/{site}/dastool_DASTool_bins"),
    params:
        out_prefix     = BINS_DIR + "/dastool/{site}/dastool",
        search_engine  = config["params"]["dastool"]["search_engine"],
        score_threshold = config["params"]["dastool"]["score_threshold"],
        extra           = config["params"]["dastool"]["extra"],
    conda:
        "../envs/dastool.yaml"
    threads:
        config["resources"]["dastool"]["cpus"]
    resources:
        cpus    = config["resources"]["dastool"]["cpus"],
        mem_mb  = config["resources"]["dastool"]["mem_mb"],
        runtime = config["resources"]["dastool"]["runtime"],
    log:
        BINS_DIR + "/logs/dastool/{site}.log",
    shell:
        """
        DAS_Tool \
            -i {input.metabat2},{input.maxbin2},{input.vamb},{input.concoct} \
            -l metabat2,maxbin2,vamb,concoct \
            -c {input.fasta} \
            -o {params.out_prefix} \
            --search_engine {params.search_engine} \
            --score_threshold {params.score_threshold} \
            {params.extra} \
            --threads {threads} \
            2> {log}
        """
