use strict;
use warnings;

use Module::Build;

my $builder = Module::Build->new(
    module_name => 'errrt',
    license => 'perl',
    dist_author => 'atb <abelew@umd.edu>',
    dist_abstract => 'Search for reverse transcriptase mutations.',
    create_makefile_pl => 0,   #'traditional',
    script_files => ['script/errrt.pl', 'script/step1.pl', 'script/step2.pl',
                     'script/step3.pl', 'script/step4.pl'],
    create_readme => 0,
    requires => {
        'File::Basename' => '0',
        'Bio::SeqIO' => '0',
        'Modern::Perl' => '0',
        'Getopt::Long' => '0',
    },
    build_requires => {
        'Test::More' => '0',
    },
);

$builder->create_build_script();
