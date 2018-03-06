#!/usr/bin/perl

# Run Spades
# Runs multiple instances of Spades back to back and then parses the Quast outputs to help determine the best assembly
# Requires Spades, Quast and optionally, SolexaQA++ in the system path to function correctly.

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

# General variable declarations
my $pe1;
my $pe2;
my $pes1;
my $pes2;
my $dir;
my $read_len;
my $v;
my $clean;

my $usage = "run_spades.pl\n
PURPOSE: Denovo assembly of paired-end sequencing reads.
This program assumes that Spades and Quast (optionally, SolexaQA++) in the system path.\n
USAGE:	run_spades.pl -1 reads1.fastq -2 reads2.fastq -d output_directory [-s singlereads.fastq -rl INT -v -clean]
-1	paired-end forward reads in FASTQ format
-2	paired-end reverse reads in FASTQ format
-s1	single-end reads 1 in FASTQ format (enter 'NA' if not using)
-s2	single-end reads 2 in FASTQ format (enter 'NA' if not using)
-d	name of the directory to store final output files\n
OPTIONAL PARAMETERS\n
-rl	read length (must be a positive integer)
-v	verbose
-clean	remove intermediate files.  **Keeps the FASTG graph files and the Quast reports**\n\n";

GetOptions(	'1=s' => \$pe1,
			'2=s' => \$pe2,
			'single1|s1=s' => \$pes1,
			'single2|s2=s' => \$pes2,
			'dir|d=s' => \$dir,
			'rl' => \$read_len,
			'clean' => \$clean,
			'verbose|v' => \$v,
) or die "$usage";

# Parameter checks
#--------------- Do read files exist?
die "One or more of the paired-end sequence read files don't exist!\n$!\n" unless ( -e $pe1 && -e $pe2 );

if ( defined($pes1) ) 	{
	die "Singleton sequence reads file 1 doesn't exist!\n$!\n" unless ( -e $pes1 );
}

if ( defined($pes2) ) 	{
	die "Singleton sequence reads file 2 doesn't exist!\n$!\n" unless ( -e $pes2 );
}

# Set up optional parameters
$clean = defined($clean)? 1 : 0;
$v = defined($v)? 1 : 0;
$pes1 = defined($pes1)? $pes1 : 'NA';
$pes2 = defined($pes2)? $pes2 : 'NA';

# Create the directory structure
print STDERR "Creating the directory structure.\n" if ( $v == 1 );
system("mkdir -p $dir\/reads $dir\/tmp $dir\/output");
system("mkdir -p $dir\/output\/fastg $dir\/output\/fasta $dir\/output\/quast");

# Copy the read files over to the reads dir for now (delete these later to save space)
system("cp $pe1 $dir\/reads\/pe1.fastq");
system("cp $pe2 $dir\/reads\/pe2.fastq");
system("cp $pes1 $dir\/reads\/pes1.fastq") unless ( $pes1 eq 'NA' );
system("cp $pes2 $dir\/reads\/pes2.fastq") unless ( $pes2 eq 'NA' );

# Estimate read lenght if not defined by the user
#---------------- Is read legnth a positive integer less than 602?
unless ( $read_len )	{
	print STDERR "Read length was undefined.  Estimating read length with SolexaQA++.\n" if ( $v == 1 );	
	system("SolexaQA++ analysis $dir\/reads\/pe1.fastq -d $dir\/tmp -h 20 --torrent");
	 
	$read_len = qx/tail -n 1 $dir\/tmp\/pe1.fastq.segments | awk '{ print \$1 }'/;
	print STDERR "Estimated read length: $read_len\n" if ( $v == 1 );
}
die "Error: unspecified or invalid read length!\n$usage" unless ( defined($read_len) && $read_len <= 602 && $read_len > 0 );

# Run Spades with multiple kmers and save output to denovo/tmp
print STDERR "Running spades with multiple kmer values.\n"  if ( $v == 1 );

#--------------------- If using optional singleton reads:
if ( defined($pes1) && defined($pes2) ) 	{
	system("spades.py --careful -k 41 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq --pe1-s $dir\/reads\/pes1.fastq --pe1-s $dir\/reads\/pes2.fastq -o $dir\/tmp\/spades41") if ( $read_len >= 100 );
	system("spades.py --careful -k 77 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq --pe1-s $dir\/reads\/pes1.fastq --pe1-s $dir\/reads\/pes2.fastq -o $dir\/tmp\/spades77") if ( $read_len >= 150 );
	system("spades.py --careful -k 99 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq --pe1-s $dir\/reads\/pes1.fastq --pe1-s $dir\/reads\/pes2.fastq -o $dir\/tmp\/spades99") if ( $read_len >= 150 );
	system("spades.py --careful -k 127 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq --pe1-s $dir\/reads\/pes1.fastq --pe1-s $dir\/reads\/pes2.fastq -o $dir\/tmp\/spades127") if ( $read_len >= 200 );
}

#--------------------- If NOT using singleton reads:
elsif ( not defined($pes1) )	{
	system("spades.py --careful -k 41 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq -o $dir\/tmp\/spades41") if ( $read_len >= 100 );
	system("spades.py --careful -k 77 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq -o $dir\/tmp\/spades77") if ( $read_len >= 150 );
	system("spades.py --careful -k 99 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq -o $dir\/tmp\/spades99") if ( $read_len >= 150 );
	system("spades.py --careful -k 127 --pe1-1 $dir\/reads\/pe1.fastq --pe1-2 $dir\/reads\/pe2.fastq -o $dir\/tmp\/spades127") if ( $read_len >= 200 );
}

# Run Quast for each assembly
print STDERR "Running Quast on the assemblies.\n" if ( $v == 1 );
system("quast.py -m 1000 -s -o $dir\/output\/quast\/q41 $dir\/tmp\/spades41\/scaffolds.fasta") if ( -e "$dir\/tmp\/spades41\/scaffolds.fasta" );
system("quast.py -m 1000 -s -o $dir\/output\/quast\/q77 $dir\/tmp\/spades77\/scaffolds.fasta") if ( -e "$dir\/tmp\/spades77\/scaffolds.fasta" );
system("quast.py -m 1000 -s -o $dir\/output\/quast\/q99 $dir\/tmp\/spades99\/scaffolds.fasta") if ( -e "$dir\/tmp\/spades99\/scaffolds.fasta" );
system("quast.py -m 1000 -s -o $dir\/output\/quast\/q127 $dir\/tmp\/spades127\/scaffolds.fasta") if ( -e "$dir\/tmp\/spades127\/scaffolds.fasta" );

# Parse the Quast output and create the results table
print STDERR "Creating the Quast matrix file.\n" if ( $v == 1 );
my $matrix = "$dir\/output\/quast\/assembly_metrics.tsv";
open MATRIX, ">", $matrix or die "Cant create matrix file. $!\n";

#------------------------- Print the column headers
print MATRIX "Assembly\t# contigs (>= 0 bp)\t# contigs (>= 1000 bp)\t# contigs (>= 5000 bp)\t# contigs (>= 10000 bp)\t# contigs (>= 25000 bp)\t# contigs (>= 50000 bp)\tTotal length (>= 0 bp)\tTotal length (>= 1000 bp)\tTotal length (>= 5000 bp)\tTotal length (>= 10000 bp)\tTotal length (>= 25000 bp)\tTotal length (>= 50000 bp)\t# contigs\tLargest contig\tTotal length\tGC (%)\tN50\tN75\tL50\tL75\t# N's per 100 kbp\n";

#------------------------- Append the actual data
my @quasts = qw(q41 q77 q99 q127);
foreach my $quast ( @quasts )	{

	my $report = "$dir\/output\/quast\/$quast\/report.tsv";
	my @vars;
	my $var = qx/awk 'BEGIN { FS="\t" } ; { print \$2 }' $report/;
	$var =~ s/\n/\t/g;
	push @vars, $var;
	
	print MATRIX (join"\t", @vars), "\n";
}
close MATRIX;

# Move the important files to the output folders from the tmp directory
print STDERR "Moving useful files to output directory.\n" if ( $v == 1 );
#----------------- FastG files
system("mv $dir\/tmp\/spades41\/assembly_graph.fastg $dir\/output\/fastg\/spades41.fastg") if ( -e "$dir\/tmp\/spades41\/*.fastg" );
system("mv $dir\/tmp\/spades77\/assembly_graph.fastg $dir\/output\/fastg\/spades41.fastg") if ( -e "$dir\/tmp\/spades77\/*.fastg" );
system("mv $dir\/tmp\/spades99\/assembly_graph.fastg $dir\/output\/fastg\/spades41.fastg") if ( -e "$dir\/tmp\/spades99\/*.fastg" );
system("mv $dir\/tmp\/spades127\/assembly_graph.fastg $dir\/output\/fastg\/spades41.fastg") if ( -e "$dir\/tmp\/spades127\/*.fastg" );
#----------------- FastA files
system("mv $dir\/tmp\/spades41\/scaffolds.fasta $dir\/output\/fasta\/spades41_scaffolds.fasta") if ( -e "$dir\/tmp\/spades41\/scaffolds.fasta" );
system("mv $dir\/tmp\/spades77\/scaffolds.fasta $dir\/output\/fasta\/spades77_scaffolds.fasta") if ( -e "$dir\/tmp\/spades77\/scaffolds.fasta" );
system("mv $dir\/tmp\/spades99\/scaffolds.fasta $dir\/output\/fasta\/spades99_scaffolds.fasta") if ( -e "$dir\/tmp\/spades99\/scaffolds.fasta" );
system("mv $dir\/tmp\/spades127\/scaffolds.fasta $dir\/output\/fasta\/spades127_scaffolds.fasta") if ( -e "$dir\/tmp\/spades127\/scaffolds.fasta" );

# Clean up the tmp files if --clean is activated
print STDERR "Cleaning up temporary files\n" if ( $clean == 1 && $v == 1 );
if ( $clean == 1 )	{
	system("rm -Rf $dir\/reads $dir\/tmp");
}

print STDERR "De novo assembly with Spades has been completed.  Exiting.\n";
exit;
