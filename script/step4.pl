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
##my $template = 'GTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTC';

my $template = 'TATCCACTGGCTACATGAACCGCCACCAGGATAATTTTTCCTTCTAGATGTGTGCAATCTAGTTGCCATATTCCTGGACTACAGTCTACTTGTCCATGCATGGCTTCTCCTTTTAGCTGACATTTATCACAGCTGGCTACTATTTCTTTTGCTACTACAGGTGGTAAGTTAAAATCACTAGCCATGGCTCTCCAATTACTGTGATATTTCTCATGTTCTTCTTGGGCCTTATCTATTCCACTGTCTCTTATACACATCTCCGAGCCCACGAGAC';
my $input = 'step3.txt.xz';
my $output = 'step4.txt';
my $verbose;
my $debug = undef;
GetOptions(
    "input=s" => \$input, # string
    "output=s" => \$output,
    "debug" => \$debug,
    "template" => \$template,
    "verbose" => \$verbose) # flag
    or die("Error in command line arguments\n");
print "Parsing ${input} and writing to ${output}\n";
my $indices = errrt::Parse_Fasta(
    input => $input,
    output => $output,
    template => $template);
