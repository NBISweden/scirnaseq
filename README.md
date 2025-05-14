# sciRNAseq3 pipeline

A [Snakemake](https://snakemake.readthedocs.io/) pipeline for processing
single-cell combinatorial indexing RNA sequencing (sci-RNA-seq) data.

sci-RNA-seq is described here:

> Martin, B.K., Qiu, C., Nichols, E. *et al.*
> Optimized single-nucleus transcriptional profiling by combinatorial indexing.
> *Nat Protoc* **18**, 188–207 (2023).
> https://doi.org/10.1038/s41596-022-00752-0

The pipeline presented here was developed independently.


## Installation


### On rackham (UPPMAX)

On the Swedish [UPPMAX](https://www.uppmax.uu.se/) cluster rackham,
Conda is preinstalled with Bioconda pre-configured,
so these commands suffice:

    module load conda
    conda env create -n scirnaseq -f environment.yml
    conda activate scirnaseq

The first and last step need to be repeated when
starting a new shell (for example, after logging in again).


### Everywhere else

1. Install Conda
2. [Configure the Bioconda channel](http://bioconda.github.io/#usage)
3. Create the Conda environment:

       conda env create -n scirnaseq -f environment.yml

4. Activate the environment:

       conda activate scirnaseq

The last step of activating the environment needs to be repeated when
starting a new shell (for example, after logging in again).


## Before running the pipeline

The reference and annotations need to be prepared.
This only needs to be done once.
These commands work for GRCm39:

    mkdir -p ref/GRCm39
    cd ref/GRCm39
    wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
    gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    gunzip Mus_musculus.GRCm39.110.gtf.gz
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm39.110.gtf
    cd ../..


## Running the pipeline

### Prepare the pipeline run directory

1. Create a new directory somewhere and change
   into it:

       mkdir run01
       cd run01

2. Make the compressed FASTQ files available in a `raw-reads/` subdirectory.
   The files must have the extension `.fastq.gz` and contain the text `_R1_`
   in their name. We recommended using symbolic links.
3. Make the indexed reference available, for example using a symbolic link:

       ln -s ../ref ref

4. Create these files:
   * `ligation-indices-with-linker.fasta`
   * `ligation-indices.fasta`
   * `p7-indices.fasta`
   * `rt-indices.fasta`

5. Create a samplesheet in TSV format named `samples.tsv`:
   It needs to have these columns:
   - `Barcode`: The 10-nt RT index
   - `Sample`: Sample name
   - `Replicate`: Replicate name (for example, use numbers 1, 2, 3, ...)


### Running Snakemake

Do a dry-run to test whether everything is set up correctly:

    snakemake --dry-run -p -s path/to/Snakefile

Replace `path/to/Snakefile` appropriately (for example, `../Snakefile`).
If there are no error messages, run the pipeline:

    snakemake -p --cores=all -s path/to/Snakefile


## Pipeline result files

All result files are placed into a newly created `out/` directory.
These are the relevant ones:

* `report.html`: QC report
* `p7-mismatches.tsv`: Statistics with P7 index mismatches.
  Also shown in `report.html`.
* `rds/`: Directory with Seurat objects (`.rds` files), one for each sample
  and one named ``allsamples.rds`` containing all samples.
* `star/Solo.out/Gene/filtered/`:
  Directory with Seurat-compatible counts matrix with overall counts.
  This is not split by sample and uses default STAR solo filtering.
* `filtered/`: Count matrix using adjusted STAR solo filtering that should
  result in more cells.
* `samples.pdf`: Various plots:
  - Number of cells per sample and replicate
  - Proportion of replicate per sample
  - nFeature_RNA
  - nCount_RNA
  - Top expressed genes


## Read structure

See the Nature Protocols paper (https://doi.org/10.1038/s41596-022-00752-0)
and also
https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq_family.html.

### Structure of R1

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
- I2 is the sample index

The I1 and I2 reads do not need to exist as FASTQ files.
I2 is ignored by the pipeline.


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

* The original pipeline published along with the Nature Protocols paper:
  https://github.com/JunyueC/sci-RNA-seq3_pipeline
