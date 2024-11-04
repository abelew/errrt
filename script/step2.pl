#!/usr/bin/env perl
use strict;
use Modern::Perl;
use autodie qw":all";
use warnings;
use diagnostics;
use Cwd;
use Getopt::Long;
use errrt;

use Bio::SeqIO;
use File::Basename;
use FileHandle;
use Storable qw"store retrieve nstore";

## Our expected template sequence (including the PCR primers)

## These are the sequences for the first(in vitro) experiment.
##my $template = 'GTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTC';
##my $final = 'CTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC';

## These are the sequences used in the second experiment.
my $template = 'TATCCACTGGCTACATGAACCGCCACCAGGATAATTTTTCCTTCTAGATGTGTGCAATCTAGTTGCCATATTCCTGGACTACAGTCTACTTGTCCATGCATGGCTTCTCCTTTTAGCTGACATTTATCACAGCTGGCTACTATTTCTTTTGCTACTACAGGTGGTAAGTTAAAATCACTAGCCATGGCTCTCCAATTACTGTGATATTTCTCATGTTCTTCTTGGGCCTTATCTATTCCACTGTCTCTTATACACATCTCCGAGCCCACGAGAC';
my $final = 'CGCCACCAGGATAATTTTTCCTTCTAGATGTGTGCAATCTAGTTGCCATATTCCTGGACTACAGTCTACTTGTCCATGCATGGCTTCTCCTTTTAGCTGACATTTATCACAGCTGGCTACTATTTCTTTTGCTACTACAGGTGGTAAGTTAAAATCACTAGCCATGGCTCTCCAATTACTGTGATATTTCTCATGTTCTT';

my $length = 9;
my $padding = 6;
my $prefix = 14;
my $library_name = 'template.fasta';
my $input = 'step1.extendedFrags.fastq.xz';
my $output = 'step2.fasta';
my $verbose = 1;
GetOptions("length=i" => \$length, # numeric
           "prefix=i" => \$prefix,
           "input=s" => \$input, # string
           "output=s" => \$output,
           "padding=s" => \$padding,
           "template=s" => \$template,
           "final=s" => \$final,
           "verbose" => \$verbose) # flag
    or die("Error in command line arguments\n");
my $terminal_string = substr($template, $padding, $length);
print "Using ${terminal_string} to extract reads from ${input}.\n";
## This terminal string is from the 5' of the forward template sequence above.  We will use it as a handle to pull off the 14 nt before each template read.
my $read_indices = errrt::Extract_Reads(
    input => $input,
    output => $output,
    final => $final,
    terminal => $terminal_string,
    prefix => $prefix,
);
