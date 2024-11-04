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
## These are the sequences in the first experiment
##my $template = 'GTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTC';
## my $final = 'CTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC';

## These are the sequences used in the second experiment.
my $template = 'TATCCACTGGCTACATGAACCGCCACCAGGATAATTTTTCCTTCTAGATGTGTGCAATCTAGTTGCCATATTCCTGGACTACAGTCTACTTGTCCATGCATGGCTTCTCCTTTTAGCTGACATTTATCACAGCTGGCTACTATTTCTTTTGCTACTACAGGTGGTAAGTTAAAATCACTAGCCATGGCTCTCCAATTACTGTGATATTTCTCATGTTCTTCTTGGGCCTTATCTATTCCACTGTCTCTTATACACATCTCCGAGCCCACGAGAC';
my $final = 'CGCCACCAGGATAATTTTTCCTTCTAGATGTGTGCAATCTAGTTGCCATATTCCTGGACTACAGTCTACTTGTCCATGCATGGCTTCTCCTTTTAGCTGACATTTATCACAGCTGGCTACTATTTCTTTTGCTACTACAGGTGGTAAGTTAAAATCACTAGCCATGGCTCTCCAATTACTGTGATATTTCTCATGTTCTT';

my $library_name = 'template.fasta';
my $fasta_cmd = 'fasta36';
my $fasta_args = '-3 -d 1 -b 1 -Z 1 -a -w 250 -n ';
my $input = 'step2.fasta';
my $output = 'step3.txt.xz';
my $verbose;
GetOptions(
    "input=s" => \$input,
    "output=s" => \$output,
    "template=s" => \$template,
    "library_name=s" => \$library_name,
    "final=s" => \$final,
    "fasta_cmd=s" => \$fasta_cmd,
    "fasta_args=s" => \$fasta_args);

print "Running ${fasta_cmd} with ${output} and ${library_name}.\n";
my $fasta_run = errrt::Run_Fasta(
    input => $input,
    output => $output,
    library_name => $library_name,
    library_seq => $template,
    fasta_cmd => $fasta_cmd,
    fasta_args => $fasta_args,
);
