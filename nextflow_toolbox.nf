
/*---------------------------------------------------------------------------------------
This Nextflow pipeline contains common programming constructs you are likely to 
need when constructing a Nextflow pipeline script.
 
 1) Initialize a Channel from a comma-separated text file.
 2) Create a Process using this input Channel.
 3) How to group Channel elements/samples together.
 4) How to merge two Channels into one.
 5) How to group replicates into two input groups for a pairwise comparison.
 6) How to ensure all process results have been generated before starting a process.

----------------------------------------------------------------------------------------*/

/* Define some variables used throughout the pipeline. These would orginarily be 
   declared in the parameters.yml file specficied on the command-line when running
   this script.
*/
designFile = "data/design.csv"
contrastsFile = "data/contrasts.csv"
publishMode = "copy"
publishOverwrite = true
RSEM_INDICES = "/camp/stp/babs/working/data/genomes/homo_sapiens/ensembl/GRCh38/release-86/genome_idx/rsem/star/75bp/genome"
FAI = "/camp/stp/babs/working/data/genomes/homo_sapiens/ensembl/GRCh38/release-86/genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai"
fastq_home = "/camp/stp/babs/working/software/nextflow/nextflow_toolbox/data"

//  1) To initialize a Channel from a comma-separated text file.
//     ---------------------------------------------------------

/* First we define our design file containing all the sample information.
   The input parameter 'input_file' is set to the path of this file. 
   This can be set on the command-line when running the Nextflow script or 
   defined in the parameters.yml file.
*/

// Define a new variable with the design file path.
File inputFile = new File( designFile )

// Check if the input file exists.
if( !inputFile.exists() ) throw new Exception( "input file " + inputFile.getPath() + " does not exist")

/* Initialize the Channel from which the Nextflow pipeline will run.

   The format of the input file is:
   sampleid,fastq file 1,groupid 

   This code creates a Channel where each iteration of the Channel will return 
   a Set containing the sample id, the group id and fastq file path.
   Sets provide a way of linking pieces of information together a bit like an array. 
   Here the sample id, group id and the fastq file are linked and will be returned
   together when the Channel is accessed by a Process.
*/

initialChannel = Channel 
  .fromPath( inputFile )
  .splitCsv( skip: 1 )
  .map { row ->
   [ row[0],
     row[2],
     new File( fastq_home + "/" + row[1] ) ]
   }

/* Here we use the Channel factory 'fromPath' to access the data contained within the file.
   There are other factory functions. 
   You can explicitly define the contents using 'from',
   read files directly using 'fromFilePairs' and
   create event driven processing using 'watchPath'

   https://www.nextflow.io/docs/latest/channel.html#

-------------------------------------------------------------------------*/

/*  2) Create a Process using this input Channel.
    --------------------------------------------*/

/*  Process 1. This Process takes the fastq files from
    initialChannel and randomly samples them for 1e3 reads.
    The sampled fastq files are used to define the output Channel.
  
    https://www.nextflow.io/docs/latest/process.html
*/

process sampleFastq {

  // This variable will be used to label the process on the system (slurm).
  tag { sampleid } 
  
  // Explicitly set the number of cpus required for this process.
  cpus 1
  // Set any other slurm command-line arguments.
  clusterOptions '--mem-per-cpu=10G'

  // Set which input Channel to use and define variables.
  input: set sampleid, groupid, fastqFile from initialChannel

  /* Define the output channels using the sampled fastq files.
     You need to create as many output channels as there are processes
     using these data. Here we have created two. 
     One as input into fastqc and the other into rsem */
  output: set sampleid, groupid, "*.sampled.fastq.gz" into \
          sampledFastqCh_fastqc, sampledFastqCh_rsem

  // The bash script that Nextflow will run given the input Channel parameters.
  // Nextflow variables are references like this ${variable}.
  // Be sure to escape any $ that you wish to be interpreted by bash. 
  script:

  // Defines the beginning of the bash script.
  """ 

  zcat ${fastqFile} | \
  awk '{ printf(\"%s\",\$0); n++; if(n%4==0) { printf(\"\\n\");} else { printf(\"\\t\\t\");} }' | \
  shuf | \
  head -n 1000 | \
  sed 's/\\t\\t/\\n/g' | \
  gzip > ${sampleid}.sampled.fastq.gz
  
  """
}

/*-------------------------------------------------------------
  Here we run fastqc on each of the sampled fastq files using
  the sampleFastqCh_fastqc Channel as input.
  The results are copied from the Nextflow output directory
  structure to the directory qc in the Nextflow run directory using
  the publishDir directive. This is useful for moving files you want to
  keep. Files and directories that you want published must be referenced
  in output channels.
*/
process fastqc {

  // This variable will be used to label the process on the system (slurm).
  tag { sampleid } 

  // Load fastqc module.
  module 'FastQC/0.11.7-Java-1.8.0_162'
  
  // Explicitly set the number of cpus required for this process along with any other.
  cpus 1
  // Slurm command-line arguments.
  clusterOptions '--mem-per-cpu=6G'

  // Write fastqc output to publish directory. These settings can be set in the parameters file
  publishDir "qc", mode: "copy", overwrite: true

  // Set which input Channel to use and define variables.
  input: set sampleid, groupid, fastqfile from sampledFastqCh_fastqc

  // Define the output Channel.
  output: set sampleid, groupid, "*" into fastqcCh

  // The bash script that Nextflow will run given the input Channel parameters.
  script:

  """ 
  fastqc -o ./ ${fastqfile}
  
  """
}

/*------------------------------------------------------------------
  Here we run rsem on each of the sampled fastq files. 32 cores are 
  used and gene/transcript quantification files and STAR genome bam files
  are used to create separate output Channels. These files are copied from 
  the Nextflow run tree to the alignment directory for saving.
*/
process rsem {

  tag { sampleid } 

  // Load rsem module.
  module 'RSEM/1.3.0-foss-2016b'
  module 'STAR/2.5.2a-foss-2016b'
  
  cpus 1
  clusterOptions '--mem-per-cpu=6G'
  
  // Publish STAR genome BAM and transcript/gene level quantifications.
  publishDir "alignment", mode: publishMode, overwrite: "true"

  // Set which input Channel to use and define variables.
  input: set sampleid, groupid, fastqfile from sampledFastqCh_rsem

  /* Define the output Channels. Here we have one for the gene-level quantifications,
     one for the transcript-level quantifications and one for the STAR aligned bam files */
  output:
    set sampleid, groupid, "*.genes.results" into rsem_genesCh
    set sampleid, groupid, "*.isoforms.results" into rsem_transcriptsCh
    set sampleid, groupid, "*.STAR.genome.bam" into \
        rsem_starbamCh_mergeBams, rsem_starbamCh_bam2bigwig, rsem_starbamCh_flagstat
  
  script:

  """
  rsem-calculate-expression \
  --temporary-folder tmp \
  --star \
  --num-threads 1 \
  --estimate-rspd \
  --seed 1 \
  --star-output-genome-bam \
  --star-gzipped-read-file \
  ${fastqfile}	\
  ${RSEM_INDICES} \
  ${sampleid}
  """
}

/* Runs flagstat on the rsem star genome files 
*/

process flagstat {

  tag{ sampleid }

  module 'SAMtools/1.5-foss-2016b'

  cpus 1
  clusterOptions '--mem-per-cpu=6G'

  publishDir "qc", mode: publishMode, overwrite: "true"

  input: set sampleid, groupid, bamfile from rsem_starbamCh_flagstat

  output: set sampleid, groupid, "*.flagstat" into flagstatCh

  script:

  """
  samtools flagstat ${bamfile} > ${sampleid}.flagstat
  """
}

/*------------------------------------------------------------------------------*/

//  3) How to group Channel elements/samples together.
//     -----------------------------------------------

/* In some instances it is desirable to group some of the samples together
   for analysis in replicate groups or merging of files. Here we use the map
   operator to map set elements contained in the rsem_starbamCh_mergeBams
   Channel to a new tuple. The first element is now the groupid. 
   The operator groupTuple merges the contents of rsem_starbamCh_mergeBams based 
   on this first element. 
*/

rsem_starbamCh_mergeBams
   .map{ s -> tuple( s[1], s[0], s[2] ) }
   .groupTuple()
   .set{ rsem_starbamCh_mergeBams }

// split groupid merged bam file channel for replicate processing below.
groupedBamsCh_1 = Channel.create()
groupedBamsCh_2 = Channel.create()
groupedBamsCh_3 = Channel.create()
rsem_starbamCh_mergeBams
  .into(
    groupedBamsCh_1,
    groupedBamsCh_2,
    groupedBamsCh_3 )

/* This process takes the newly merged star bam channel and merges the 
   files using samtools merge 
*/

process mergeBamFiles {

  tag { groupid }

  module 'SAMtools/1.5-foss-2016b'
  
  cpus 1
  clusterOptions '--mem-per-cpu=6G'
  
  input: set groupid, sampleids, file( bamfiles ) from groupedBamsCh_1

  output: set groupid, "*.bam" into mergedBamFilesCh
    
  script:

  """
  samtools merge ${groupid} ${bamfiles}
  """
}

/*-------------------------------------------------------------------*/

//  4) How to merge two Channels into one.
//     -----------------------------------

/* Sometimes the output from more than one process is required as input 
   into another. To achieve this and ensure synchronicity of the samples 
   we need to merge channels together. In this example we take the output 
   from samtools flagstat and the corresponding bam file as input to create a 
   bigwig file. The flagstat data are used to normalise the bigwig coverage
   vectors. Here one of the rsem starbam channels is merged with the flagstat
   channel using the first element in the sets from both channels using the join 
   operator.
*/

rsem_starbamCh_bam2bigwig
  .map{ s -> tuple( s[0], s[1], s[2] ) }
  .join( flagstatCh )
  .set{ rsem_starbamCh_bam2bigwig }

/* This process takes the rsem star genome bams and flag stat reports
   and creates normalised bigwig coverage files
*/
process bam2bigwig {

  tag{ s }

  module 'BEDTools/2.26.0-foss-2016b'
  module 'Kent_tools/20161115-linux.x86_64'
  
  cpus 1
  clusterOptions '--mem-per-cpu=6G'

  publishDir "visualisation", mode: publishMode, overwrite: true

  input:
    set s, groupid1, bamfile, groupid2, flagstatfile from rsem_starbamCh_bam2bigwig
    
  output:
    file "*.bigwig" into bam2bigwigCh

  script:

  // Total reads for scaling factor
  def totalreads = 1000000

  """
  # calculate scaling factor
  fragcount=`grep mapped ${flagstatfile} | head -n 1 | awk '{print \$1}'` 
  export scale=`bc <<< \"scale=6; ${totalreads}/\$fragcount\"`
  
  genomeCoverageBed \
   -ibam ${bamfile} \
   -bg \
   -pc \
   -scale \$scale \
   -g ${FAI} > ${s}.bedgraph
  
  wigToBigWig \
   -clip \
   ${s}.bedgraph \
   ${FAI} \
   ${s}.bigwig

  rm ${s}.bedgraph
  """
}

/*-------------------------------------------------------------------------*/

//  5) How to group replicates into two input groups for a pairwise comparison.
//     ------------------------------------------------------------------------


/* Merging replicate groups so both groups are available together is a bit more
   tricky. First we need to create an input channel from a contrasts file that defines 
   the pair-wise comparisons.

   This file defines a pairwise comparison per line using group ids as defined 
   above. The format used here is comparison_label,treatment_groupid,control_groupid. 

   The channel contains each comparison as a set entry, similar to how the sample 
   information is read from the design file above.
*/

contrastsCh = Channel
  .fromPath( contrastsFile )
  .splitCsv( skip: 1 )
  .map { row ->
   [ row[ 0 ], row[ 1 ], row[ 2 ] ] }

/* We then create separate treatment and control channels by first splitting this
   new contrasts channel into two
*/
contrastsCh_treatment = Channel.create()
contrastsCh_control = Channel.create()
contrastsCh
  .into( contrastsCh_treatment,
	 contrastsCh_control )

/* Replicate group samples are then joined to these channels from a data channel 
   created above using the groupid. In this case we are using bam files provided by 
   the groupedBamCh created above.
*/
contrastsCh_treatment
  .map{ s -> tuple( s[1], s[0] ) }
  .join( groupedBamsCh_2 )
  .map{ s -> tuple( s[1], s[0], s[3] ) }
  .set{ contrastsCh_treatment_groupedBams }

contrastsCh_control
  .map{ s -> tuple( s[2], s[0] ) }
  .join( groupedBamsCh_3 )
  .map{ s -> tuple( s[1], s[0], s[3] ) }
  .set{ contrastsCh_control_groupedBams }

/* We now have two channels containing grouped bam files indexed by groupid. 
   One channel for our treatment groups and one for our corresponding control 
   groups.
*/

/* We can now merge these two channels into one using the comparisonid.
   This will ensure the comparison groups are paired as defined in the 
   contrasts file. The channel contrastsWithBamsCh can now be used as 
   input into a process that compares the two groups.
*/
contrastsCh_treatment_groupedBams
  .map{ s -> tuple( s[0], s[1], s[2] ) }
  .join( contrastsCh_control_groupedBams )
  .set{ contrastsWithBamsCh }

/* Toy process to demonstrate how the replicate grouped 
   comparison channel can be used.
*/
process compareBams {

  errorStrategy 'ignore'
  
  tag{ contrast }

  cpus 1
  clusterOptions '--mem-per-cpu=6G'
 
  input:
    set contrast, t_target, file( t_bamfiles ), c_target, file( c_bamfiles ) from contrastsWithBamsCh

  output: file "*.output" into compareBams

  script:
  
  """
  compareBams ${t_bamfiles} ${c_bamfiles} > ${contrast}.output 
  """
}

/* A final user case takes care of processes that require all iterations of a 
   previous process to finish before running. A good example of this 
   is multiqc. Here we can use the collect operator to ensure all
   results are generated before the process starts. This multiqc example
   will wait until all fastqc processes have finished.
*/

process multiqc {
  
  module "multiqc/1.3-2016b-Python-2.7.12"
  
  cpus 1
  clusterOptions '--mem-per-cpu=6G'
  
  publishDir "qc", mode: publishMode, overwrite: true
  
  input:
    file( "*" ) from fastqcCh.collect()

  output:
    file "multiqc*" into multiqcQCCh 

  script:
  """
  multiqc ../../../qc
  """
}