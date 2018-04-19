
# nextflow_toolbox

A toy Nextflow pipeline that contains most of the programming constructs you will need to
start writing your own pipelines.

nextflow: http://www.nextflow.io
nextflow-quickstart: http://www.nextflow.io/docs/latest/getstarted.html#get-started

### Overview

This Nextflow pipeline contains code to carry out the following steps:
 
 1) Initialize a Channel from a comma-separated text file.
 2) Create a Process using this input Channel.
 3) Group Channel elements/samples together.
 4) Merge two Channels into one.
 5) Group replicates into two input groups for a pairwise comparison.
 6) Ensure all process results have been generated before starting a process.

### The Pipeline Flow

![alt flow's directed acyclic graph][nextflow_toolbox_dag]

### Description

The pipeline takes the data/design.csv as input. It samples the fastq files for 1e3 reads. It runs fastqc and rsem on these sampled files. Samtools flagstat is run on the bam files. These bam files are converted to bigwig using the read counts extracted from flagstat for normalisation. The bam files are grouped by replicate group. This grouping is then used to merge the bam files. The file contrasts.csv is read in to define the contrasts to run. The merged bam files are then organised into their contrast groups for processing by the process compareBams. More detailed documentation is available in the script.

### Running

You can run the pipeline using the following command. By default itis configured for a slurm environment.
    git clone https://github.com/crickbabs/nextflow_toolbox
    cd rnaSeq_byBABS
    module load nextflow/0.27.2
    nextflow run nextflow_toolbox.nf

### Data

The test data provided was obtained from GEO submission 