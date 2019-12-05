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
my $input = 'merged_fasta.txt.gz';
my $output = 'merged_fasta_parsed.txt';
my $verbose;
my $debug = undef;
GetOptions(
    "input=s" => \$input, # string
    "output=s" => \$output,
    "debug" => \$debug,
    "verbose" => \$verbose) # flag
    or die("Error in command line arguments\n");
print "Parsing ${input} and writing to ${output}\n";
my $indices = errrt::Parse_Fasta(input => $input, output => $output);
