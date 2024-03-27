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
    input: "out/report.html", "out/samples.pdf"


rule trim_ligation_index:
    output:
        fastq=temp("out/round1/{name}_{extra}.fastq.gz"),
        json="out/round1/{name}_{extra}.cutadapt.json",
    input:
        r1_fastq="raw-reads/{name}_R1_{extra}.fastq.gz",
        r2_fastq="raw-reads/{name}_R2_{extra}.fastq.gz",
        ligation_indices_with_linker_fasta="ligation-indices-with-linker.fasta",
    log:
        "out/round1/{name}_{extra}.cutadapt.log"
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
        fastq=temp("out/round2/{name}.fastq.gz"),
        json="out/round2/{name}.cutadapt.json",
    input:
        fastq="out/round1/{name}.fastq.gz",
        rt_indices_fasta="rt-indices.fasta",
    log:
        "out/round2/{name}.cutadapt.log"
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
        #" --untrimmed-output out/untrimmed2/${name}.1.fastq.gz"
        #" --untrimmed-paired-output out/untrimmed2/${name}.2.fastq.gz"
        " {input.fastq}"
        " > {log}"


rule write_allowed_barcodes:
    output:
        allowed_barcodes_txt="out/allowed-barcodes.txt"
    input:
        rt_indices_fasta="rt-indices.fasta",
        ligation_indices_fasta="ligation-indices.fasta",
        p7_indices_fasta="p7-indices.fasta",
    params: script=Path(workflow.basedir) / "allowedbarcodes.py"
    shell:
        "python3 {params.script}"
        " {input.rt_indices_fasta}"
        " {input.ligation_indices_fasta}"
        " {input.p7_indices_fasta}"
        " > {output.allowed_barcodes_txt}.tmp"
        " &&"
        " mv {output.allowed_barcodes_txt}.tmp {output.allowed_barcodes_txt}"


rule simulate_cell_barcode:
    output:
        r1_fastq=temp("out/cb-reads/{name}.1.fastq.gz"),
        r2_fastq=temp("out/cb-reads/{name}.2.fastq.gz"),
        p7_mismatch_stats_tsv="out/cb-reads/{name}.p7mismatches.tsv"
    input:
        fastq="out/round2/{name}.fastq.gz",
        rt_indices_fasta="rt-indices.fasta",
        ligation_indices_fasta="ligation-indices.fasta",
        p7_indices_fasta="p7-indices.fasta",
    params: script=Path(workflow.basedir) / "simulatecb.py"
    shell:
        "python3 {params.script}"
        " --r1 {output.r1_fastq}"
        " --r2 {output.r2_fastq}"
        " --p7-stats {output.p7_mismatch_stats_tsv}"
        " {input.rt_indices_fasta}"
        " {input.ligation_indices_fasta}"
        " {input.p7_indices_fasta}"
        " {input.fastq}"


rule merge_p7_mismatch_stats:
    output:
        tsv="out/p7-mismatches.tsv",
    input:
        expand("out/cb-reads/{name}.p7mismatches.tsv", name=NAMES)
    shell:
        "head -n1 {input[0]} > {output}; "
        "for f in {input}; do sed 1d $f; done >> {output}"


rule star_solo:
    output:
        "out/star/Solo.out/Gene/Summary.csv",
        "out/star/Solo.out/Gene/UMIperCellSorted.txt",
        expand(
            "out/star/Solo.out/Gene/{which}/{name}",
            which=("filtered", "raw"),
            name=("matrix.mtx", "features.tsv", "barcodes.tsv")
        ),
    input:
        allowed_barcodes_txt="out/allowed-barcodes.txt",
        ref="ref/GRCm39/star/Genome",
        r1_fastqs=expand("out/cb-reads/{name}.1.fastq.gz", name=NAMES),
        r2_fastqs=expand("out/cb-reads/{name}.2.fastq.gz", name=NAMES),
    params:
        r1_fastqs=",".join(f"out/cb-reads/{name}.1.fastq.gz" for name in NAMES),
        r2_fastqs=",".join(f"out/cb-reads/{name}.2.fastq.gz" for name in NAMES),
    threads: 99
    shell:
        "STAR"
        " --genomeDir ref/GRCm39/star"
        " --readFilesCommand zcat"
        " --readFilesIn {params.r2_fastqs} {params.r1_fastqs}"
        " --runThreadN {threads}"
        " --outFileNamePrefix out/star/"
        " --outSAMtype BAM Unsorted"
        " --soloType CB_UMI_Simple"
        " --soloCBwhitelist {input.allowed_barcodes_txt}"
        " --soloCBstart 1"
        " --soloCBlen 30"
        " --soloUMIstart 31"
        " --soloUMIlen 8"
        " --outSAMattributes NH HI nM AS CR UR GX GN"


rule filter_cells:
    output:
        expand("out/filtered/{name}", name=("matrix.mtx", "features.tsv", "barcodes.tsv")),
    input:
        expand("out/star/Solo.out/Gene/raw/{name}", name=("matrix.mtx", "features.tsv", "barcodes.tsv"))
    params: script=Path(workflow.basedir) / "filtercells.py"
    log: "out/filtered/filtercells.log"
    shell:
        "python3 {params.script}"
        " --umis 200"
        " out/star/Solo.out/Gene/raw/"
        " out/filtered/"
        " 2> {log}"


rule saturation:
    output:
        "out/saturation.tsv"
    input:
        bam="out/star/Aligned.out.bam",
        barcodes="out/filtered/barcodes.tsv"
    params: script=Path(workflow.basedir) / "saturation.py"
    shell:
        "python3 {params.script}"
        " --barcodes {input.barcodes}"
        " {input.bam}"
        " > {output}"


rule report:
    output: "out/report.html"
    input:
        "out/p7-mismatches.tsv",
        "out/star/Solo.out/Gene/Summary.csv",
        "out/star/Solo.out/Gene/UMIperCellSorted.txt",
        "out/saturation.tsv",
    params: qmd=Path(workflow.basedir) / "report.qmd"
    shell:
        "cp -n {params.qmd} out/ && "
        "quarto render out/report.qmd --execute-dir out/ --to html --toc && "
        "rm out/report.qmd"


rule split_samples_and_plot:
    output:
        directory("out/per-sample/"),
        pdf="out/samples.pdf",
    input:
        expand("out/filtered/{name}", name=("matrix.mtx", "features.tsv", "barcodes.tsv")),
        samples_tsv="samples.tsv",
    params: script=Path(workflow.basedir) / "sample_split.R"
    shell:
        "Rscript {params.script}"
        " -i out/filtered/"
        " -m {input.samples_tsv}"
        " --pdf {output.pdf}"
        " -o out/per-sample/"
