#!/usr/bin/perl

##################################################################
#            Trim and QC Raw Sequence Reads
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
my $open_html;
my $clean;
my $v;

# GetOpts required variables
my $dir;
my $pe1;
my $pe2;

my $usage = "trim_qc_reads.pl\n\n
PURPOSE: Pre-assembly quality control of raw sequence data.\n
This program assumes that FastQC, Seqtk, Scythe, and SolexaQA++ are all in the system path.\n\n
USAGE:	reference_based_assembly.pl\n
-1	paired-end forward reads in FASTQ format\n
-2	paired-end reverse reads in FASTQ format\n
-d	name of the directory to store final output files\n
OPTIONAL PARAMETERS\n\n
-open	automatically opens the FastQC HTML files in Google Chrome
-v	verbose\n
-clean	remove intermediate files.  **Keeps the sorted BAM file and BAM index**\n\n";

GetOptions(	'1=s' => \$pe1,
			'2=s' => \$pe2,
			'dir|d=s' => \$dir,
			'open' => \$open_html,
			'clean' => \$clean,
			'verbose|v' => \$v,
) or die "$usage";

# Parameter setup
$clean = defined($clean)? 1 : 0;
$v = defined($v)? 1 : 0;
$open_html = defined($open_html)? 0 : 1;

# Create directory structure and copy the raw read data into it for processing.
print STDERR "Creating the directory structure for data handling.\n" if ( $v == 1 );
system("mkdir -p $dir\/reads\/{raw,tmp,trimmed}");
system("cp $pe1 $dir\/reads\/raw\/read1.fastq");
system("cp $pe2 $dir\/reads\/raw\/read2.fastq");

# Run initial FastQC on the raw reads
print STDERR "QC-ing raw reads with FastQC.\n" if ( $v == 1 );
system("fastqc $dir\/reads\/raw\/*.fastq");
system("google-chrome *.html") if ( $open_html == 1 );

# Force trim library bias with Seqtk
print STDERR "Force trimming read data.\n" if ( $v == 1 );
system("seqtk trimfq -b 15 -e 5 $dir\/reads\/raw\/read1.fastq > $dir\/reads\/tmp\/read1.fastq");
system("seqtk trimfq -b 15 -e 5 $dir\/reads\/raw\/read2.fastq > $dir\/reads\/tmp\/read2.fastq");

# Quality trim with Seqtk  (need to add in functionality for different Q scores)
print STDERR "Quality trimming read data at Q20.\n" if ( $v == 1 );
system("seqtk trimfq -q 0.01 $dir\/reads\/tmp\/read1.fastq > $dir\/reads\/tmp\/trimmed1.fastq");
system("seqtk trimfq -q 0.01 $dir\/reads\/tmp\/read2.fastq > $dir\/reads\/tmp\/trimmed2.fastq");

# Remove Illumina Nextera adapter files (add in feature to allow input of different adapter fasta file)
print STDERR "Removing Nextera adapters.\n" if ( $v == 1 );
system("scythe -a /media/hunter/Data/Bioinformatics/scythe/nextera_transposase.fasta -o $dir\/reads\/tmp\/adapters1.fastq $dir\/reads\/tmp\/trimmed1.fastq");
system("scythe -a /media/hunter/Data/Bioinformatics/scythe/nextera_transposase.fasta -o $dir\/reads\/tmp\/adatpers2.fastq $dir\/reads\/tmp\/trimmed2.fastq");

# Length sort the trimmed reads (add in feature to sort by user-specified length)
print STDERR "Sorting the trimmed reads by length.\n" if ( $v == 1 );
system("SolexaQA++ lengthsort -l 25 -d $dir\/reads\/trimmed $dir\/reads\/tmp\/adapters*.fastq");

# Move and rename the length sorted files for ease of processing
system("mv $dir\/reads\/tmp\/adapters1.fastq.paired $dir\/reads\/trimmed\/adapters1.paired.fastq");
system("mv $dir\/reads\/tmp\/adapters1.fastq.single $dir\/reads\/trimmed\/adapters1.single.fastq");
system("mv $dir\/reads\/tmp\/adapters2.fastq.paired $dir\/reads\/trimmed\/adapters2.paired.fastq");
system("mv $dir\/reads\/tmp\/adapters1.fastq.discard $dir\/reads\/trimmed\/adapters1.discard.fastq");
system("mv $dir\/reads\/tmp\/adapters1.fastq.summary.txt $dir\/reads\/trimmed\/lengthsort_summary.txt");
system("mv $dir\/reads\/tmp\/adapters1.fastq.summary.txt.pdf $dir\/reads\/trimmed\/lengthsort_summary.pdf");

# Run the final FastQC on the clean data
print STDERR "Running final FastQC on cleaned data.\n" if ( $v == 1 );
system("fastqc $dir\/reads\/trimmed\/*.paired.fastq");
system("google-chrome *.html") if ( $open_html == 1 );

# Clean the raw and intermediate data files
print STDERR "Cleaning up the raw and intermediate data files.\n" 0f ( $v == 1 );
system("rm -Rf reads/raw reads/tmp");

print STDERR "Data has been cleaned.  Exiting.\n\n";
exit;
