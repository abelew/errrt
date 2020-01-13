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

my $flash = '/cbcb/sw/RedHat-7-x86_64/common/local/flash/1.2.11/bin/flash';
my $output = 'step1';
my $flash_args = ' --compress-prog=xz --compress-prog-args=-9e';
my $read1 = '';
my $read2 = '';
GetOptions("read1=s" => \$read1,
           "read2=s" => \$read2,
           "flash=s" => \$flash,
           "output=s" => \$output,
           "flash_args=s" => \$flash_args,
       );

print "Running ${flash} with ${flash_args} on ${read1} and ${read2}.\n";
my $merged = errrt::Merge_Reads(
    read1 => $read1,
    read2 => $read2,
    flash => $flash,
    output => $output,
    flash_args => $flash_args);
