#!/usr/bin/perl

##################################################################
#            Reference guided assembly using BWA
#          Assumes reads have already been trimmed/QC
#  BWA, SAMtools, and Seqtk are required to be in the system path for this to work properly
##################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

# General variable declarations
my $refdir_ref;
my $reference_basename;
my $clean;
my $v;

# GetOpts required variables
my $dir;
my $pe1;
my $pe2;
my $reference;
my $output;

my $usage = "reference_based_assembly.pl\n\n
PURPOSE: Automating reference-based genome assembly.\n
This program assumes that BWA, Samtools, and Seqtk are all in the system path.\n\n
USAGE:	reference_based_assembly.pl\n
-1	paired-end forward reads in FASTQ format\n
-2	paired-end reverse reads in FASTQ format\n
-r	reference genome file in FASTA format\n
-d	name of the directory to store final output files\n
-o	final output file name\n\n
OPTIONAL PARAMETERS\n\n
-v	verbose\n
-clean	remove intermediate files.  **Keeps the sorted BAM file and BAM index**\n\n";

GetOptions(	'1=s' => \$pe1,
			'2=s' => \$pe2,
			'ref|r=s' => \$reference,
			'dir|d=s' => \$dir,
			'out|o=s' => \$output,
			'clean' => \$clean,
			'verbose|v' => \$v,
) or die "$usage";

# Parameter setup
$clean = defined($clean)? 1 : 0;
$v = defined($v)? 1 : 0;

# Make output dir(s)
print STDERR "Creating output directory $dir.\n" if ( $v == 1 );
system("mkdir -p $dir");
system("mkdir -p $dir\/ref");
system("mkdir -p $dir\/tmp");
system("mkdir -p $dir\/output");

# Copy the reference file to the new dir for indexing
$reference_basename = fileparse($reference);
$refdir_ref = "$dir\/ref\/$reference_basename";
system("cp $reference $refdir_ref");

# Index the reference file
print STDERR "Indexing the reference file provided: $reference\n" if ( $v == 1 );
system("bwa index -a is $refdir_ref");

# Separately align each paired end read file
print STDERR "Aligning the paired end reads to the reference.\n" if ( $v == 1 );
system("bwa aln $refdir_ref $pe1 > $dir\/tmp\/pe1-1.sai");
system("bwa aln $refdir_ref $pe2 > $dir\/tmp\/pe1-2.sai");

# Combine the alignments into a single SAM file
print STDERR "Generating SAM file from aligned reads.\n" if ( $v == 1 );
system("bwa sampe -f $dir\/tmp\/pe1-12.sam $refdir_ref $dir\/tmp\/pe1-1.sai $dir\/tmp\/pe1-2.sai $pe1 $pe2");

# Convert SAM to BAM format
print STDERR "Converting SAM file to BAM format\n" if ( $v == 1 );
system("samtools view -@ 8 -S $dir\/tmp\/pe1-12.sam -b -o $dir\/tmp\/pe1-12.bam");

# Sort and index the BAM file
print STDERR "Sorting and indexing the BAM file.\n" if ( $v == 1 );
system("samtools sort -@ 8 -o $dir\/tmp\/pe1-12_sorted.bam $dir\/tmp\/pe1-12.bam");
system("samtools index -@ 8 $dir\/tmp\/pe1-12_sorted.bam");

# Create consensus file with mpileup
print STDERR "Creating the consensus file.\n" if ( $v == 1 );
system("samtools mpileup -v --no-BAQ -f $refdir_ref $dir\/tmp\/pe1-12_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A - > $dir\/output\/$output");
#--------------This command may take some time to run
print STDERR "File converstion complete.  Output Fasta is saved as $dir\/output\/$output.\n" if ( $v == 1 );

# Run Quast to QC the output assembly
print STDERR "Running assembly QC with Quast.\n" if ( $v == 1 );
system("quast.py -t 8 -m 1000 -e -o $dir\/output\/ref_assembly_quast $dir\/output\/$output");
system("google-chrome $dir\/output\/ref_assembly_quast\/report.html");

# Run alnstat.py script to print some basic stats about the BAM file we mapped
print STDERR "Generating alignment statistics from the BAM file.\n\n";
system("alnstat.py $dir\/tmp\/pe1-12_sorted.bam");
print STDERR "\n\n";

# Clean up the index and intermediate files if option is turned on
print STDERR "Cleaning up temporary files\n" if ( $clean == 1 && $v == 1 );
if ( $clean == 1 )	{
	system("rm -Rf $dir\/ref");
	system("mv $dir\/tmp\/pe1-12_sorted.bam $dir\/output\/pe1-12_sorted.bam");
	system("mv $dir\/tmp\/pe1-12_sorted.bam.bai $dir\/output\/pe1-12_sorted.bam.bai");
	system ("rm -Rf $dir\/tmp");
}

print STDERR "Reference-based assembly completed and Quast results saved.  Exiting.\n";
exit;
