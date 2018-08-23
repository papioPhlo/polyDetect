# polyDetect
polyDetect is used to detect and genotype polymorphic Alu insertions from short-read (<1KB) data. This program was designed to identify polymorphic Alu elements by evaluating a panel of primate samples to a reference genome. Thus, it should be used to process and evaluate multiple individuals, each with their own sequencing file. polyDetect must be run separately on each individual, producing its own unique output directory. The input can be either paired-end or single-end fastq files. However, it is recommended that if the data are publicly available on the NCBI-SRA database that you use fastQgenie to create your fastq file(s).
Date: June 29, 2018
Author: Vallmer Jordan
 
Pre-requisites for polyDetect==============
-Python 3 installation
-Linux based 64-bit machine
-fastq file
	-If you wish to create your fastq file using fastQgenie* (in this repository-	see below for instructions) (recommended) you will need:
		-sratoolkit.2.6.2 or higher (make sure to add to your path)
		-nesoni (make sure to add to your path)
	-If you wish to run polyDetect with a pre-existing fastq file, you will need to check the two sample fastq files located in this repository to make sure your fastq is formatted properly. You're files match the format of these sample files identically, including the OP:i:N tag and the name of each sequence must be a number with 12 digits (first sequence should be named 000000000001).
-samtools (make sure to add to your path)
-bwa (make sure to add to your path)
-bowtie2 (make sure to add to your path)
-sam.py (in this repository)*
-sra.py (in this repository)*
-essentials.py (in this repository)*
-local_stats.py (in this repository)*
-meiDetect.py (in this repository)*
-polyAnna.py (in this repository)*

-index for all fasta files used (mobile element consensus sequence (a.k.a., Alu), polyA sequence, reference genome sequence) using all three of the following programs:
	-bwa
	Usage:
		$bwa index <sequence.fa>
		-generates .amb .ann .bwt .pac .sa files
	-samtools
	Usage:
		$samtools faidx <sequence.fa>
		-generates .fai file
	-bowtie2
	Usage:
		$bowtie2-build –f <sequence.fa> --large-index <outfile title>
		-generates .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2 files
		Example:
			$bowtie2-build Homo_sapiens.fa H_sap
			-output: H_sap.1.bt2 H_sap.2.bt2 H_sap.3.bt2 H_sap.4.bt2 H_sap.rev.1.bt2 H_sap.rev.2.bt2 files
			(This outfile title will be used to show how to run 						polyDetect.py)		
			NOTE: Your reference genome, prior to building the index, must be masked using the repeat you are looking for. This can be accomplished using the samtools function, maskfasta. You will need your FASTA reference genome and the .bed file of the location of your repeats in the reference genome. 
	

*It is recommended that these files be placed in their own directory, and that that directory is added to your path

Running fastQgenie ===================
You must first construct a text file containing a list of SRA Run accessions from a single biological sample. A sample SRR file is provided in this repository (sample.srr). From the same directory containing your custom SRR file, execute the fastQgenie script.
	Example:
		$ cd /directory_of_sra file
		$ fastQgenie.py (no commands)

Running polyDetect ===================
The polyDetect.py script should be executed from an empty directory followed by 8 mandatory commands:
python polydetect.py <alu> <ref> <pe.fq> <se.fq> <out> <id> chrom> <polyA>
<alu> - path consensus alu sequence (fasta format)**
<ref> - path to bowtie2 index (outfile title)**
<pe.fq> - paired-end interleaved fastq file***
<se.fq> - single-end fastq file***
<out> - path to output directory
<id> - Sample ID: custom string used to identify this sample. Must consist of letters and/or numbers (no spaces or special characters). This string will be the prefix for your output files.
<chrom> - total number of chromosomes in your reference genome
<polyA> - polyA consensus sequence.** The standard Alu polyA consensus sequence recommended for this program can be found in this repository.

	Example:
		$ polyDetect.py /directory_for_**/aluConsensus.fa directory_for_**/H_sap /directory_for_***/filename.PE.fq /directory_for_***/filename.SE.fq /outfile_directory /organism_ID chromosome_number  /directory_for_**/polyA.fa

**It is recommended that these files be placed in their own directory (with all of the other indexes generated for the fasta files)
***It is recommended that these files be placed in their own directory with the sample ID as the directory name

Running dataEval.py ===================
Make sure that for each sample, a directory exists containing only the files generated using polyDetect for that sample. The name of each directory should be the appropriate sample id <id>. Must match the prefix of the files contained within. Now place all of these directories in one master directory.
In this master directory, create a file containing a list of all of the sample IDs (See sample orgID list in this repository: sample.orgID)
From the master directory, run dataEval.py
dataEval.py <master> <idlist> <predict> <chrom> > <output>
<master> - path to master directory
<idlist> - path to orgID list
<predict> - path to file containing actual sequences spanning the point of insertion
<chrom> - total number of chromosomes in your reference genome
<output> - path to file contain output genotype datafile

