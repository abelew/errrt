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

my $flash = 'flash';
my $output = 'step1';
my $flash_args = ' --compress-prog=xz --compress-prog-args=-9e';
my $correct = 0;
my $read1 = '';
my $read2 = '';
GetOptions("read1=s" => \$read1,
           "read2=s" => \$read2,
           "correct" => \$correct,
           "flash=s" => \$flash,
           "output=s" => \$output,
           "flash_args=s" => \$flash_args,
       );

print "Running ${flash} with ${flash_args} on ${read1} and ${read2} correction:$correct.\n";
my $merged = errrt::Merge_Reads(
    correct => $correct,
    read1 => $read1,
    read2 => $read2,
    flash => $flash,
    output => $output,
    flash_args => $flash_args);
