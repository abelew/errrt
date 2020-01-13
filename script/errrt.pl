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

## Step 1 options: merge read 1 and read 2 into a single pile of sequences.
my $read1 = '';
my $read2 = '';
my $flash = '/cbcb/sw/RedHat-7-x86_64/common/local/flash/1.2.11/bin/flash';
my $flash_args = ' -d 1 -b 1 --compress-prog=xz --compress-prog-args=-9e ';
my $step1_out = 'step1';
## Step 2 options: find reads containing at least the first n($length) nucleotides of the template.
my $template = 'GTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTC';
my $final = 'CTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC';
my $length = 9;
my $prefix = 14;
my $step2_out = 'step2';
## Step 3 options: Use fasta36 to look for mutant reads.
my $fasta_cmd = '/cbcb/sw/RedHat-7-x86_64/common/local/fasta/36.3.5/bin/fasta36';
my $fasta_args = '-3 -d 1 -b 1 -Z 1 -a -w 250 -n ';
my $step3_out = 'step3';
## Step 4 options: Parse the fasta36 output.
my $step4_out = 'step4';
my $verbose;

GetOptions(
    ## Step 1 args
    "read1=s" => \$read1,
    "read2=s" => \$read2,
    "flash=s" => \$flash,
    "flash_args=s" => \$flash_args,
    "step1_out=s" => \$step1_out,
    ## Step 2 args
    "length=i" => \$length, # numeric
    "prefix=i" => \$prefix,
    "template=s" => \$template,
    "final=s" => \$final,
    "step2_out=s" => \$step2_out,
    ## Step 3 args
    "fasta_cmd=s" => \$fasta_cmd,
    "fasta_args=s" => \$fasta_args,
    "step3_out=s" => \$step3_out,
    ## Step 4 args
    "step4_out=s" => \$step4_out,
    "verbose" => \$verbose);


## Perform Step 1
print "Running ${flash} with ${flash_args} on ${read1} and ${read2}.\n";
my $merged = errrt::Merge_Reads(
    read1 => $read1,
    read2 => $read2,
    flash => $flash,
    flash_args => $flash_args,
    output => $step1_out,
);

## Perform Step 2
my $terminal_string = substr($template, 0, $length);
print "Using ${terminal_string} to extract reads.\n";
## This terminal string is from the 5' of the forward template sequence above.
## We will use it as a handle to pull off the 14 nt before each template read.
my $read_indices = errrt::Extract_Reads(
    input => $step1_out,
    final => $final,
    terminal => $terminal_string,
    prefix => $prefix,
    output => $step2_out,
);

## Perform Step 3
print "Running ${fasta_cmd}.\n";
my $fasta_run = errrt::Run_Fasta(
    input => $step2_out,
    library_seq => $template,
    fasta_cmd => $fasta_cmd,
    fasta_args => $fasta_args,
    output => $step3_out,
);

## Perform Step 4
print "Parsing the fasta36 output.\n";
my $indices = errrt::Parse_Fasta(
    input => $step3_out,
    template => $template,
    output => $step4_out,
);
