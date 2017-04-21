Introduction
==================
Specific Insertions Detector (SID) is a program to detect non-reference human transposon insertion. It is compiled in Perl and includes two steps, discordant reads detection and reads clustering. Generally, the first step collects informative reads and generates other necessary files, while the second step discovers the specific insertion sites and exports the final results into a plain text.

Memory requirement
==================
In the first step, SID usually consumes less than 1GB memory. While its peak usage of memory could up to 30GB with one thread using 30-70X whole genome sequencing data of human in the second step.


Dependent softwares
==================
1. Samtools v1.0 or later ( earlier versions may result in some mistakes in the step of discordant reads detection) 
2. Perl Module: Bio::DB::Sam, Statistics::Descriptive, threads::shared, IO::File 
3. BLAST (v2.2.25 or later)


Input files preparation
==================
1. BAM file: paired-end sequencing data aligned by BWA aln
2. FASTA file: The sequence of non-reference TEs and human reference genome. Both of them need BLAST indices. 


Publication
==================
We are preparing for a paper using this program.

Version
==================
1.   01discordant_v2.pl: v2.0
2.   02cluster.pl: v1.0

Tips
==================
1. The program allows putting 1 or more BAM files in a BAM_list file (plain text) as input. 
2. You must make a BLAST index for the TE sequence, and put it in the same directory with TE FASTA file.
3. When running this program, the input BAM file should not remove duplicates beforehand, or it may stop running accidentally.
4. The parameter of '-run' cannot be used at present, and we will fix it soon.
5. Please export the path of BLAST and Samtools to .bashrc before running SID.

Demo case
==================
An example of how to use SID and the demo input and output of SID. I uploaded the demo case to Google.
Please contact me if any questions: yuqichao@genomics.cn or yqc20101111@hotmail.com

https://drive.google.com/file/d/0B-5j9b_mSd_GaEI1eUxvamh2SWM/view?usp=sharing



