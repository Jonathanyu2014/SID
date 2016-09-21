Introduction
==================
Specific Insertions Detector (SID) is a Perl program to detect non-reference human transposon insertion. Through discordant reads detection and reads clustering, it could detect non-reference RIPs easily and quickly.

Memory requirement
==================
Its peak usage of memory is 28GB with one thread using 30-70X whole genome sequencing data of human.


Dependent softwares
==================
1. Samtools v0.1.18 (version 1.1 or later may result in some mistakes in the step of discordant reads detection) 
2. Perl Module: Bio::DB::Sam, Statistics::Descriptive, threads::shared, IO::File 
3. BLAST 


Input files preparation
==================
BAM file: paired-end sequencing data aligned by BWA 
FASTA file: The sequence of non-reference TEs

NOTEs
==================
1. The program allows putting 1 or more BAM files in a BAM_list file (plain text) as input. 
2. You must make a BLAST index for the TE sequence and put it in the same directory with TE FASTA file.


Publication
==================
We are preparing for a paper using this program.

Other tips
==================
1. When running this program, the input BAM file should not remove duplicates beforehand, or it may stop running suddenly. 
2. The parameter of '-run' cannot be used at present, and we will fix it soon.
