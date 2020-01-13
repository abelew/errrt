#!/usr/bin/env perl
use lib '/mnt/sshfs/cbcbsub/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/rt_erates_destefano_2019/errrt/lib';
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
my $template = 'GTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTC';
my $final = 'CTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC';
my $length = 9;
my $prefix = 14;
my $library_name = 'template.fasta';
my $input = 'step1.extendedFrags.fastq.xz';
my $output = 'step2.fasta';
my $verbose;
GetOptions("length=i" => \$length, # numeric
           "prefix=i" => \$prefix,
           "input=s" => \$input, # string
           "output=s" => \$output,
           "template=s" => \$template,
           "final=s" => \$final,
           "verbose" => \$verbose) # flag
    or die("Error in command line arguments\n");
my $terminal_string = substr($template, 0, $length);
print "Using ${terminal_string} to extract reads from ${input}.\n";
## This terminal string is from the 5' of the forward template sequence above.  We will use it as a handle to pull off the 14 nt before each template read.
my $read_indices = errrt::Extract_Reads(
    input => $input,
    output => $output,
    final => $final,
    terminal => $terminal_string,
    prefix => $prefix,
);
