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
my $fasta_cmd = '/cbcb/sw/RedHat-7-x86_64/common/local/fasta/36.3.5/bin/fasta36';
my $fasta_args = '-d 1 -b 1';
my $input = 'merged.extendedFrags.fastq.xz';
my $verbose;
GetOptions("length=i" => \$length, # numeric
           "prefix=i" => \$prefix,
           "input=s" => \$input, # string
           "template=s" => \$template,
           "library_name=s" => \$library_name,
           "final=s" => \$final,
           "fasta_cmd=s" => \$fasta_cmd,
           "fasta_args=s" => \$fasta_args,
           "verbose" => \$verbose) # flag
    or die("Error in command line arguments\n");
my $output = basename($input, ('.xz', '.gz', '.bz2'));
$output = basename($output, ('.fastq'));
$output = basename($output, ('.extendedFrags'));
$output = qq"${output}.fasta";
my $terminal_string = substr($template, 0, $length);
print "Using ${terminal_string} to extract reads from ${input}.\n";
## This terminal string is from the 5' of the forward template sequence above.  We will use it as a handle to pull off the 14 nt before each template read.
my $read_indices = errrt::Extract_Reads(
    input => $input,
    output => $output,
    terminal => $terminal_string,
    prefix => $prefix,
);

print "Running ${fasta_cmd} with ${output} and ${library_name}.\n";
my $fasta_run = errrt::Run_Fasta(
    input => $output,
    library_name => $library_name,
    library_seq => $template,
    fasta_cmd => $fasta_cmd,
    fasta_args => $fasta_args,
);
my $fasta = basename($output, ('.fasta'));
my $fasta_result = "${fasta}_fasta.txt";

print "Parsing ${fasta_result}\n";
my $indices = errrt::Parse_Fasta(input => $fasta_result);

print "Counting results\n";
my $counts = errrt::Count_Parsed(input => $indices);
