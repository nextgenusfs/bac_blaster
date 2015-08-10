#!/usr/bin/env perl

# AUTHOR: Jon Palmer
# LAST REVISED: March 2015
# palmer3 at wisc dot edu
# All rights reserved.

# script to process folder of ABI trace files and BLAST against target
# Dependencies include: Python 2.7, Biopython, BLAST, tracetuner in scripts directory

use strict;
use warnings;
use Getopt::Long;

#set command line arguments
my ($folder, $db, $out) = @ARGV;
my $version="bac_blaster.pl  v0.2.1";
my $helpAsked;
my $qual = 0;
my $cpus = 1;
my $length = 60;
my $osname = $^O;
GetOptions(
	'f|folder:s'=>\$folder, #folder containing ABI trace files
	'd|database:s'=>\$db, #name of Blast database, must be in DB folder, i.e. Nectria
  	'o|output:s'=>\$out, #name of file output
	'q|qual:i'=>\$qual, #quality score to trim reads, default is 15
	'c|cpus:i'=>\$cpus, #number of cpus to use for blast search, default is 1
	'l|length:i'=>\$length, #min length to keep quality trimmed reads
	'h|help:s'=>\$helpAsked,
	'v|version'=>sub{print $version."\n"; exit;},
);

if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(!defined($folder)) {
	prtUsage();
	exit;
}

if($qual == 0){
	$qual = 15;
}

#read ABI file names into array to capture all files
opendir(DIR, $folder) or die "couldn't open $folder: $!\n";
my @filenames = grep { /\.ab1$/ } readdir DIR;
closedir DIR;

#get short name from output
my $shortname = $out;
$shortname =~ s{\.[^.]+$}{};

#Sub routine for running tracetuner windows
sub wintracetuner{
 print "-------------------------------------------------------------\n";
 print "Now Running TraceTuner to extract reads from ABI files.\n";
my $tracer = `scripts/tracetuner_3.0.6beta/rel/x86-win32/ttuner.exe -V -id $folder -sa output.fasta -qa output.qual`;
}
#Sub routine for running tracetuner Unix
sub unixtracetuner{
 print "-------------------------------------------------------------\n";
 print "Now Running TraceTuner to extract reads from ABI files.\n";
my $tracer = `scripts/tracetuner_3.0.6beta/rel/Linux/ttuner -V -id $folder -sa output.fasta -qa output.qual`;
}
#Sub routine for running python conversion script
sub winfasta2fastq{
 print "-------------------------------------------------------------\n";
 print "Converting to FASTQ.\n";
my $fasta2fastq = `C:/Python27/python scripts/fasta2fastq.py output.fasta output.qual > output.fastq`;
}
#Sub routine for running python conversion script
sub unixfasta2fastq{
 print "-------------------------------------------------------------\n";
 print "Converting to FASTQ.\n";
my $fasta2fastq = `python scripts/fasta2fastq.py output.fasta output.qual > output.fastq`;
}
#Sub routine for running quality trimming
sub wintoolkit_trim{
 print "-------------------------------------------------------------\n";
 print "Quality trimming sequences to values > $qual and > $length length.\n";
my $quality_trim = `type output.fastq | perl scripts/trim_fastq.pl -q $qual -l $length > trim.fastq`
}
#Sub routine for running quality trimming
sub unixtoolkit_trim{
 print "-------------------------------------------------------------\n";
 print "Quality trimming sequences to values > $qual and > $length length.\n";
my $quality_trim = `cat output.fastq | perl scripts/trim_fastq.pl -q $qual -l $length > trim.fastq`
}
#Sub routine for converting back to fasta
sub winfastq2fasta{
my $conversion = `C:/Python27/python -c "from Bio import SeqIO;SeqIO.convert('trim.fastq','fastq','$shortname.fasta','fasta')" `
}
#Sub routine for converting back to fasta
sub unixfastq2fasta{
my $conversion = `python -c "from Bio import SeqIO;SeqIO.convert('trim.fastq','fastq','$shortname.fasta','fasta')" `
}
#Sub routine for running blast
sub blaster{
 print "-------------------------------------------------------------\n";
 print "Running BLAST search against $db database.";
my $blast_results = `blastn -query $shortname.fasta -db DB/$db -outfmt "6 qseqid sseqid score evalue length nident positive sstart send qlen" -out bac_blaster_results.tmp -max_target_seqs 1`;
my $tbl_results = 'bac_blaster_results.tmp';
my $output = $out;
open(OUTPUT, '>'.$output) or die;
print OUTPUT "BAC Clone\tSeq_length\tTarget found\tScore (Bits)\tExpect(E-value)\tAlign-length\tIdentities\tPositives\tChr/supercontig\tStart\tEnd\n";
close $output;
open(my $fh, "<:encoding(UTF-8)", $tbl_results) or die;
open(OUTPUT, '>>'.$output) or die;
while (<$fh>) {
	my @row = split(/\s+/,$_);
	print OUTPUT "$row[0]\t$row[9]\t$row[1]:$row[7]-$row[8]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[1]\t$row[7]\t$row[8]\n";
}
close $fh;
close $output;
}

#sub routine for help
sub prtHelp {
	print "\n$0 options: \n\n";
	print "### Input a folder name containing ABI trace files (Required)\n";
	print "  Example:  -f sample_data\n";
	print "\n";
	print "### Input a species to BLAST against (Required)\n";
	print "  Example:  -d Nectria \n";
	print "\n";
	print "Pre-installed Options are:\n";
	print "	Nectria: Nectria haematococca\n";
	print "	A_wentii: Aspergillus wentii\n";
	print "	A_aculeatus: Aspergillus aculeatus\n";
	print "	P_expansum: Penicillium expansum\n";
	print "	P_marneffei: Penicillium marneffei\n";
	print "\n";
	print "### Name of output text file (Required)\n";
	print "  Example:  -o my_blast_search.txt\n";
	print "\n";
	print "### Optional flags\n";
	print " -c number of cpus to use for BLAST search\n";
	print " -q quality score to trim reads, default is 15\n";
	print " -l minimum length read to keep, default is 60 bp\n";
	print "\n";
	print "### Get help on how to use this script:\n";
	print "   -h  Prints this help message\n";
	print "#############################################\n";
	print "Complete example of script:\n";
	print "perl bac_blaster.pl -f sample_data -d Nectria -o nectria_blast.txt";
	print "\n";
}

#print usage subroutine
sub prtUsage {
	print "OS version: $osname\n";
	print "Usage: perl $0 <options>\n";
	prtHelp();
}

#Run sub routines if windows
if( $osname !~ 'MSWin32' ){
	unixtracetuner();
	unixfasta2fastq();
	unixtoolkit_trim();
	unixfastq2fasta();
	blaster();
}
	else {
	wintracetuner();
	winfasta2fastq();
	wintoolkit_trim();
	winfastq2fasta();
	blaster();
}

print "\n";
print "-------------------------------------------------------------\n";
print "Script has successfully run.  Your BLAST results are located in $out\n";
print "and the quality filtered multi-fasta file is called: $shortname.fasta\n";
print "-------------------------------------------------------------\n";
unlink 'bac_blaster_results.tmp'; #deletes the original blast output file
unlink 'trim.fastq';
unlink 'output.fastq';
unlink 'output.qual';
unlink 'output.fasta';
