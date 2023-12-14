from pathlib import Path
import re

# https://stackoverflow.com/a/16090640/715090
def natural_sort_key(s):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split("([0-9]+)", str(s))
    ]

R1_FASTQS = sorted(Path("raw-reads/").glob("*_R1_*.fastq.gz"), key=natural_sort_key)
NAMES = [r1.name.replace("_R1_", "_").replace(".fastq.gz", "") for r1 in R1_FASTQS]

if not R1_FASTQS:
    sys.exit("No FASTQ files found in raw-reads/ directory")


rule final:
    input: "report.html", "samples.pdf"


rule trim_ligation_index:
    output:
        fastq="round1/{name}_{extra}.fastq.gz",
        json="round1/{name}_{extra}.cutadapt.json",
    input:
        r1_fastq="raw-reads/{name}_R1_{extra}.fastq.gz",
        r2_fastq="raw-reads/{name}_R2_{extra}.fastq.gz",
        ligation_indices_with_linker_fasta="ligation-indices-with-linker.fasta",
    log:
        "round1/{name}_{extra}.cutadapt.log"
    threads: 8
    shell:
        "cutadapt"
        " -Z"
        " -j {threads}"
        " --json={output.json}"
        " --interleaved"
        " -o {output.fastq}"
        " -g ^file:{input.ligation_indices_with_linker_fasta}"
        " -e 2"
        " --no-indels"
        " --rename '{{header}} ligation_index={{r1.adapter_name}}'"
        " --discard-untrimmed"
        " {input.r1_fastq} {input.r2_fastq}"
        " > {log}"


rule trim_umi_and_rt_index:
    output:
        fastq="round2/{name}.fastq.gz",
        json="round2/{name}.cutadapt.json",
    input:
        fastq="round1/{name}.fastq.gz",
        rt_indices_fasta="rt-indices.fasta",
    log:
        "round2/{name}.cutadapt.log"
    threads: 8
    shell:
        "cutadapt"
        " -Z"
        " -j {threads}"
        " --json={output.json}"
        " --interleaved"
        " -o /dev/null"
        " -p {output.fastq}"
        " -u 8"
        " -g ^file:{input.rt_indices_fasta}"
        " -e 2"
        " --no-indels"
        " --rename '{{header}} umi={{r1.cut_prefix}} rt_index={{r1.adapter_name}}'"
        " --discard-untrimmed"
        #" --untrimmed-output untrimmed2/${name}.1.fastq.gz"
        #" --untrimmed-paired-output untrimmed2/${name}.2.fastq.gz"
        " {input.fastq}"
        " > {log}"


rule write_allowed_barcodes:
    output:
        allowed_barcodes_txt="allowed-barcodes.txt"
    input:
        ligation_indices_fasta="ligation-indices.fasta",
        rt_indices_fasta="rt-indices.fasta",
        p7_indices_fasta="p7-indices.fasta",
    params: script=Path(workflow.basedir) / "allowedbarcodes.py"
    shell:
        "python3 {params.script}"
        " {input.ligation_indices_fasta}"
        " {input.rt_indices_fasta}"
        " {input.p7_indices_fasta}"
        " > {output.allowed_barcodes_txt}.tmp"
        " &&"
        " mv {output.allowed_barcodes_txt}.tmp {output.allowed_barcodes_txt}"


rule simulate_cell_barcode:
    output:
        r1_fastq="cb-reads/{name}.1.fastq.gz",
        r2_fastq="cb-reads/{name}.2.fastq.gz",
        p7_mismatch_stats_tsv="cb-reads/{name}.p7mismatches.tsv"
    input:
        fastq="round2/{name}.fastq.gz",
        ligation_indices_fasta="ligation-indices.fasta",
        rt_indices_fasta="rt-indices.fasta",
        p7_indices_fasta="p7-indices.fasta",
    params: script=Path(workflow.basedir) / "simulatecb.py"
    shell:
        "python3 {params.script}"
        " --r1 {output.r1_fastq}"
        " --r2 {output.r2_fastq}"
        " --p7-stats {output.p7_mismatch_stats_tsv}"
        " {input.ligation_indices_fasta}"
        " {input.rt_indices_fasta}"
        " {input.p7_indices_fasta}"
        " {input.fastq}"


rule merge_p7_mismatch_stats:
    output:
        tsv="p7-mismatches.tsv",
    input:
        expand("cb-reads/{name}.p7mismatches.tsv", name=NAMES)
    shell:
        "head -n1 {input[0]} > {output}; "
        "for f in {input}; do sed 1d $f; done >> {output}"


rule star_solo:
    output:
        "star-out/Solo.out/Gene/Summary.csv",
        "star-out/Solo.out/Gene/UMIperCellSorted.txt",
        expand(
            "star-out/Solo.out/Gene/{which}/{name}",
            which=("filtered", "raw"),
            name=("matrix.mtx", "features.tsv", "barcodes.tsv")
        ),
    input:
        allowed_barcodes_txt="allowed-barcodes.txt",
        ref="ref/GRCm39/star/Genome",
        r1_fastqs=expand("cb-reads/{name}.1.fastq.gz", name=NAMES),
        r2_fastqs=expand("cb-reads/{name}.2.fastq.gz", name=NAMES),
    params:
        r1_fastqs=",".join(f"cb-reads/{name}.1.fastq.gz" for name in NAMES),
        r2_fastqs=",".join(f"cb-reads/{name}.2.fastq.gz" for name in NAMES),
    threads: 99
    shell:
        "STAR"
        " --genomeDir ref/GRCm39/star"
        " --readFilesCommand zcat"
        " --readFilesIn {params.r2_fastqs} {params.r1_fastqs}"
        " --runThreadN {threads}"
        " --outFileNamePrefix star-out/"
        " --outSAMtype BAM Unsorted"
        " --soloType CB_UMI_Simple"
        " --soloCBwhitelist {input.allowed_barcodes_txt}"
        " --soloCBstart 1"
        " --soloCBlen 30"
        " --soloUMIstart 31"
        " --soloUMIlen 8"


rule star_solo_cell_filter:
    output:
        expand("filtered/{name}", name=("matrix.mtx", "features.tsv", "barcodes.tsv")),
    input:
        expand("star-out/Solo.out/Gene/raw/{name}", name=("matrix.mtx", "features.tsv", "barcodes.tsv"))
    log: "filtered/Log.out"
    shell:
        "cd filtered"
        "; "
        "STAR --runMode soloCellFiltering ../star-out/Solo.out/Gene/raw/ ./ --soloCellFilter CellRanger2.2 200000 0.85 8"

rule report:
    output: "report.html"
    input:
        "p7-mismatches.tsv",
        "star-out/Solo.out/Gene/Summary.csv",
        "star-out/Solo.out/Gene/UMIperCellSorted.txt",
    params: qmd=Path(workflow.basedir) / "report.qmd"
    shell:
        "cp -n {params.qmd} . && "
        "quarto render report.qmd --to html --toc"


rule split_samples_and_plot:
    output:
        directory("per-sample/"),
        pdf="samples.pdf",
    input:
        expand("filtered/{name}", name=("matrix.mtx", "features.tsv", "barcodes.tsv")),
        samples_tsv="samples.tsv",
    params: script=Path(workflow.basedir) / "sample_split.R"
    shell:
        "Rscript {params.script}"
        " -i filtered/"
        " -m {input.samples_tsv}"
        " --pdf {output.pdf}"
        " -o per-sample/"
