from pathlib import Path


R1_FASTQS = list(Path("reads/").glob("*_R1_*.fastq.gz"))
R2_FASTQS = [Path(str(r1).replace("_R1_", "_R2_")) for r1 in R1_FASTQS]
ROUND2_FASTQS = [Path("round2/" + str(r1.name).replace("_R1_", "_")) for r1 in R1_FASTQS]

if not R1_FASTQS:
    sys.exit("No FASTQ files found in reads/ directory")

print("Files:")
for r1, r2 in zip(R1_FASTQS, R2_FASTQS):
    print("-", r1, r2)

rule final:
    input: "report.html"


rule report:
    output: "report.html"
    input:
        "p7-mismatches.tsv",
        "star-out/Solo.out/Gene/Summary.csv",
        "star-out/Solo.out/Gene/UMIperCellSorted.txt",
    params: qmd=workflow.source_path("report.qmd")
    shell:
        "cp -n {params.qmd} . && "
        "quarto render report.qmd --to html --toc"


rule trim_ligation_index:
    output:
        fastq="round1/{name}_{extra}.fastq.gz",
        json="round1/{name}_{extra}.cutadapt.json",
    input:
        r1_fastq="reads/{name}_R1_{extra}.fastq.gz",
        r2_fastq="reads/{name}_R2_{extra}.fastq.gz",
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


rule simulate_cell_barcode:
    output:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        allowed_barcodes_txt="allowed-barcodes.txt",
        p7_mismatch_stats_txt="p7-mismatches.tsv",
    input:
        round2_fastqs=ROUND2_FASTQS,
        ligation_indices_fasta="ligation-indices.fasta",
        rt_indices_fasta="rt-indices.fasta",
        p7_indices_fasta="p7-indices.fasta",
    params: script=workflow.source_path("simulatecb.py")
    shell:
        "python3 {params.script}"
        " --list {output.allowed_barcodes_txt}"
        " --r1 {output.r1_fastq}"
        " --r2 {output.r2_fastq}"
        " --p7-stats {output.p7_mismatch_stats_txt}"
        " {input.ligation_indices_fasta}"
        " {input.rt_indices_fasta}"
        " {input.p7_indices_fasta}"
        " {input.round2_fastqs}"


rule star_solo:
    output:
        "star-out/Solo.out/Gene/Summary.csv",
        "star-out/Solo.out/Gene/UMIperCellSorted.txt",
        "star-out/Solo.out/Gene/filtered/matrix.mtx",
    input:
        allowed_barcodes_txt="allowed-barcodes.txt",
        ref="ref/GRCm39/star/Genome",
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
    threads: 99
    shell:
        "STAR"
        " --genomeDir ref/GRCm39/star"
        " --readFilesCommand zcat"
        " --readFilesIn {input.r2_fastq} {input.r1_fastq}"
        " --runThreadN {threads}"
        " --outFileNamePrefix star-out/"
        " --soloType CB_UMI_Simple"
        " --soloCBwhitelist {input.allowed_barcodes_txt}"
        " --soloCBstart 1"
        " --soloCBlen 30"
        " --soloUMIstart 31"
        " --soloUMIlen 8"
        " --outSAMtype BAM Unsorted"
