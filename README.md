# sciRNAseq3 pipeline

## Installation


1. Install Conda
2. [Configure the Bioconda channel](http://bioconda.github.io/#usage)
3. Create the Conda environment:

       conda env create -n scirnaseq -f environment.yml
4. Activate the environment:

       conda activate scirnaseq

The last step of activating the environment needs to be repeated when
starting a new shell (for example, after logging in again).


### On rackham (UPPMAX)

On the Swedish [UPPMAX](https://www.uppmax.uu.se/) cluster rackham,
Conda is preinstalled with Bioconda pre-configured,
so these commands suffice:

    module load conda
    conda env create -n scirnaseq -f environment.yml
    conda activate scirnaseq


## Preparing references

The reference and annotations need to be prepared.
This only needs to be done once.
These commands work for GRCm39:

    mkdir -p ref/GRCm39 && cd ref/GRCm39
    wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
    gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    gunzip Mus_musculus.GRCm39.110.gtf.gz
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm39.110.gtf


## Preparing the data

## Running the pipeline

## Read structure

See the Nature Protocols paper (https://doi.org/10.1038/s41596-022-00752-0)
and also
https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq3.html.


### R1

Structure:
```
#{9,10}CAGAGCN{8}#{10}
```
`#` is part of an index, `N` is part of UMI

1. 9-10 nt ligation index (some have 9 nt, some have 10 nt)
2. 6 nt linker
3. 8 nt UMI
4. 10 nt RT index

Total: 33 or 34 nt. Sequenced reads have 34 nt.


### Other reads

- R2 is the RNA read
- I1 is the third index
- I2 is the sample index (ignored by the pipeline)

The I1 and I2 reads do not need to exist as FASTQ files.


## Pipeline overview

* Input data is assumed to already have been demultiplexed by P7 index (I1)
* Cutadapt is run twice to move ligation index, UMI and RT index from R2 into the read header. R2 is empty afterwards and discarded.
* Script `simulatecb.py` simulates R2 reads similar to 10X Chromium reads.
  The sequence consists of:
  - RT index (10 nt)
  - Ligation index (10 nt, those that have only 9 nt have an `A` added to them)
  - P7 index (10 nt)
  - UMI
* STAR solo is run. It is configured to interpret the first 30 nt as cell barcode.


## Links

* An existing pipeline: https://github.com/JunyueC/sci-RNA-seq3_pipeline

Martin, B.K., Qiu, C., Nichols, E. et al. Optimized single-nucleus transcriptional profiling by combinatorial indexing. Nat Protoc 18, 188–207 (2023).
https://doi.org/10.1038/s41596-022-00752-0
