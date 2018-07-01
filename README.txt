# polyDetect
polyDetect is used to detect and genotype polymorphic Alu insertions from short-read (<1KB) data. This program was designed to identify polymorphic Alu elements by evaluating a panel of primate samples. Thus, it should be used to process and evaluate multiple individuals, each with their own sequencing file. polyDetect must be run separately on each individual, producing its own unique output directory. The input can be either paired-end or single-end fastq files. However, it is recommended that if the data are publicly available on the NCBI-SRA database that you use fastQgenie to create your fastq file(s).
Date: June 29, 2018
Author: Vallmer Jordan 
Pre-requisites ============== polyDetect is a python script and requires Python 3 to be installed. It has been tested on a Linux based (Ubuntu 14.04) 64-bit machine. Additionally, it requires the following tools:
If you wish to create your fastq file using fastQgenie (recommended) you will need:
sratoolkit.2.6.2 or higher (make sure to add to your path) 
nesoni (make sure to add to your path)
If you wish to run polyDetect with a pre-existing fastq file, you will need to check the fastq sample.fq file located in this repository to make sure your fastq is formatted properly.
To run polyDetect, you will need the following tools:
bwa (make sure to add to your path)
bowTie2 (make sure to add to your path)
Running fastQgenie ===================
You must first construct a text file containing a list of SRA Run accessions. A sample SRR file is provided in this repository (sample.srr). From the same directory containing your custom SRR file, execute the fastQgenie script.

Running polyDetect =================== 
The polyDetect.py script should be executed from an empty directory followed by 8 mandatory commands:
python polydetect.py <alu> <ref> <pe.fq> <se.fq> <out> <id> chrom> <polyA>
<alu> - path consensus alu sequence (fasta format)
<ref> - path to reference assembly (fasta format)
<pe.fq> - paired-end interleaved fastq file
<se.fq> - single-end fastq file
<out> - path to output directory
<id> - Sample ID: custom string used to identify this sample. Must consist of letters and/or numbers (no spaces or special characters). This string will be the prefix for your output files.
<chrom> - total number of chromosomes in your reference genome
<polyA> - polyA consensus sequence. The standard Alu polyA consensus sequence recommended for this program can be found in this repository.

Running dataEval.py ===================
Make sure that for each sample, a directory exists containing only the files generated using polyDetect for that sample. The name of each directory should be the appropriate sample id <id>. Much match the prefix of the files contained within. Now place all of these directories in one master directory.
In this master directory, create a file containing a list of all of the sample IDs (See sample orgID list in this repository: sample.orgID)
From the master directory, run dataEval.py
dataEval.py <master> <idlist> <predict> <chrom> > <output>
<master> - path to master directory
<idlist> - path to orgID list
<predict> - path to file containing actual sequences spanning the point of insertion
<chrom> - total number of chromosomes in your reference genome
<output> - path to file contain output genotype datafile
