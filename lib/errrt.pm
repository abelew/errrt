package errrt;
use File::Basename;
use FileHandle;
use IO::Handle;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Bio::SearchIO::fasta;
use Bio::Seq;
use Cwd;
use FileHandle;
use File::Basename qw"dirname basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use POSIX qw"ceil";
use Storable qw"store retrieve nstore";

use strict;
use 5.008_005;
our $VERSION = '0.01';

sub Extract_Reads {
    my (%args) = @_;
    my $inputted = FileHandle->new("less $args{input} |");
    my $in = Bio::SeqIO->new(-fh => $inputted, -format => 'Fastq');
    my $hist = FileHandle->new(">idx_count.txt");
    my $fasta = FileHandle->new(">$args{output}");
    my $data_file = "$args{output}.pdata";

    my $indices = {};
    my $count = 0;
    my $success = 0;
  READS: while (my $in_seq = $in->next_dataset()) {
        $count++;
        my $id = $in_seq->{'-descriptor'};
        my $sequence = $in_seq->{'-seq'};
        my $qual = $in_seq->{'-raw_quality'};
        ## Assume the $length long terminal portion of the sequence is in the forward direction.
        my $direction = "f";
        ## If we do not observe the terminal sequence in a read,
        ## Then set the direction to 'reverse', reverse complement
        ## the read, and look again.
        unless ($sequence =~ /$args{terminal}/) {
            $direction = "r";
            $sequence = reverse($sequence);
            $sequence =~ tr/ATGC/TACG/;
        }
        ## Let us now pull the ${prefix} long fragment as the index along
        ## with all the sequence which follows and print it to a fasta file
        ## which will be aligned against the template with fasta36.
        if ($sequence =~ m/^(.+)([ATGC]{$args{prefix}})($args{terminal}.*)$/) {
            my $index = $2;
            my $data = $3;
            my $success++;
            print $fasta ">${count}_${direction}_${index}
$data
";
            if (defined($indices->{$index})) {
                $indices->{$index} = $indices->{$index} + 1;
            } else {
                $indices->{$index} = 1;
            }
        } ## End checking for a read with a valid terminus.
    } ## End iterating over every read.
    foreach my $k (keys %{$indices}) {
        print $hist "$k $indices->{$k}\n";
    }
    $hist->close();
    $inputted->close();
    $fasta->close();
    print "Wrote ${success} index containing reads.\n";
    store($indices, $data_file);
    return($indices)
}

sub Run_Fasta {
    my (%args) = @_;
    my $library = FileHandle->new(">$args{library_name}");
    print $library ">template_sequence
$args{library_seq}
";
    my $output = basename($args{input}, ('.fasta'));
    my $error = qq"${output}_fasta.err";
    $output = qq"${output}_fasta.txt";

    my $runner = qq"$args{fasta_cmd} $args{fasta_args} $args{input} $args{library_name} 2>${error} | xz -9e > ${output}.xz";
    my $handle = IO::Handle->new;
    print "Starting $runner\n";
    open($handle, "$runner |");
    while (my $line = <$handle>) {
        print "$line\n";
    }
    close($handle);
    return($output);
}

sub Parse_Fasta {
    my (%args) = @_;
    my $f = FileHandle->new("less $args{input} |");
    my $output = "parsed.txt";
    my $out = FileHandle->new(">${output}");
    my $searchio = Bio::SearchIO->new(-format => 'fasta', -fh => $f);
    my $data_file = "${output}.pdata";
    my $results = 0;
    my $indices = {};
    while (my $result = $searchio->next_result()) {
        $results++;
        while (my $hit = $result->next_hit) {
            my $query_name = $result->query_name();
            my $query_length = $result->query_length();
            my $accession = $hit->accession();
            my $library_id = $hit->name();
            my $length = $hit->length();
            my $score = $hit->raw_score();
            my $sig = $hit->significance();
            my $ident = $hit->frac_identical();
            my $lstrand = $hit->strand('hit');
            my ($count, $direction, $index) = split(/_/, $query_name);
            my @hitlst = ();
            if ($indices->{$index}) {
                @hitlst = @{$indices->{$index}};
            } else {
                $indices->{$index} = [];
            }

            my $hit_len;
            my $hit_matches;
            my $first_qstart = -1;
            my $first_qend = -1;
            my $first_lstart = -1;
            my $first_lend = -1;
            my $hsp_count = 0;
            my $hits = {};

            while (my $hsp = $hit->next_hsp) {
                my @mismatch_indices = $hsp->seq_inds('query', 'mismatch');
                if (scalar(@mismatch_indices) > 0) {
                    my $q_seq = $hsp->query_string;
                    my $h_seq = $hsp->hit_string;
                    for my $ind (@mismatch_indices) {
                        my $q = substr($q_seq, $ind - 1, 1);
                        my $h = substr($h_seq, $ind - 1, 1);
                        $hits->{$ind}->{type} = 'mismatch';
                        $hits->{$ind}->{ref} = $q;
                        $hits->{$ind}->{hit} = $h;
                        ##print "MISMATCH: $query_name: $ind $q $h\n";
                    }
                }
                my @del_indices = $hsp->seq_inds('query', 'gap');
                if (scalar(@del_indices) > 0) {
                    if (scalar(@del_indices) > 0) {
                        my $q_seq = $hsp->query_string;
                        my $h_seq = $hsp->hit_string;
                        for my $ind (@del_indices) {
                            my $q = substr($q_seq, $ind - 2, 6);
                            my $h = substr($h_seq, $ind - 2, 6);
                            $hits->{$ind}->{type} = 'del';
                            $hits->{$ind}->{ref} = $q;
                            ##print "DEL: $query_name: $ind $q $h\n";
                        }
                    }
                }
                my @ins_indices = $hsp->seq_inds('hit', 'gap');
                if (scalar(@ins_indices) > 0) {
                    if (scalar(@ins_indices) > 0) {
                        my $q_seq = $hsp->query_string;
                        my $h_seq = $hsp->hit_string;
                        for my $ind (@ins_indices) {
                            my $q = substr($q_seq, $ind - 2, 6);
                            my $h = substr($h_seq, $ind - 2, 6);
                            $hits->{$ind}->{type} = 'ins';
                            $hits->{$ind}->{ref} = $h;
                            ##print "INS: $query_name: $ind $q $h\n";
                        }
                    }
                }
            }
            push(@hitlst, $hits);
            $indices->{$index} = \@hitlst;
        }   ## End of each hit of a result.
    }       ## End of each result.
    $f->close();
    store($indices, $data_file);
    return($indices);
}

sub Count_Parsed {
    my (%args) = @_;
    print "Starting counters.\n";
    my $input = $args{input};
    my @count_by_hits = ();
    foreach my $i (keys %{$input}) {
        my @hitlst = @{$input->{$i}};
        my $num_hits = scalar(@hitlst);
        if ($count_by_hits[$num_hits]) {
            $count_by_hits[$num_hits]++;
        } else {
            $count_by_hits[$num_hits] = 1;
        }
    }
    use Data::Dumper;
    my $c = -1;
    for my $count (@count_by_hits) {
        $c++;
        if ($count) {
            print "There are $count reads with $c hits.\n";
        }
    }
}

1;
__END__

=encoding utf-8

=head1 NAME

errrt - Blah blah blah

=head1 SYNOPSIS

  use errrt;

=head1 DESCRIPTION

errrt is

=head1 AUTHOR

Ashton Trey Belew E<lt>abelew@gmail.comE<gt>

=head1 COPYRIGHT

Copyright 2019- Ashton Trey Belew

=head1 LICENSE

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 SEE ALSO

=cut
