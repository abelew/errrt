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

use Data::Dumper;

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
        }                      ## End checking for a read with a valid terminus.
    }                          ## End iterating over every read.
    foreach my $k (keys %{$indices}) {
        print $hist "$k\t$indices->{$k}\n";
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
    my $output = $args{output};
    my $out = FileHandle->new(">${output}");
    my $searchio = Bio::SearchIO->new(-format => 'fasta', -fh => $f);
    my $data_file = basename($output, ('.txt', '.xz', '.gz'));
    $data_file = basename($output, ('.txt'));
    $data_file .= ".pdata";
    my $results = 0;
    my $indices = {};
  SEARCHLOOP: while (my $result = $searchio->next_result()) {
        $results++;
        ## There is only ever 1 hit.
        my $hit = $result->next_hit;
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
        }

        my $hit_len;
        my $hit_matches;
        my $first_qstart = -1;
        my $first_qend = -1;
        my $first_lstart = -1;
        my $first_lend = -1;
        my $hsp_count = 0;
        my $hits = {};

        my $hsp = $hit->next_hsp;
        my @query_mismatches = $hsp->seq_inds('query', 'mismatch');
        my @hit_mismatches = $hsp->seq_inds('hit', 'mismatch');
        my @product_gap = $hsp->seq_inds('query', 'gap');
        my @template_gap = $hsp->seq_inds('hit', 'gap');
        my $q_seq = $hsp->query_string;
        my $h_seq = $hsp->hit_string;
        ## There is an important caveat here: The indices returned by seq_inds()
        ## do not take into account the positions taken up by gap characters.
        ## As a result, if one wishes to get the changed nucleotide _after_ an
        ## indel, then any attempt to find those nucleotides is doomed to
        ## failure until one rewrites the sequence without the '-' characters.
        my $q_raw = $q_seq;
        $q_raw =~ s/\-//g;
        my $h_raw = $h_seq;
        $h_raw =~ s/\-//g;
        ## I expect the length of query_mismatches and hit_mismatches to be identical.
        if (scalar(@query_mismatches) > 0) {
            ##print "HMMMMM: 1       9        19        29        39        49        59        69        79        89        99\n";
            ##print "QESTME: $q_seq\n";
            ##print "HESTME: $h_seq\n";
            for my $i (0 .. $#query_mismatches) {
                my $q_miss = $query_mismatches[$i] - 1;
                my $h_miss = $hit_mismatches[$i] - 1;
                my $qseq = substr($q_raw, $q_miss, 1);
                my $hseq = substr($h_raw, $h_miss, 1);
                ##print "TESTME1: $qseq $q_miss $hseq $h_miss\n";
                $hits->{$q_miss}->{type} = 'mismatch';
                $hits->{$q_miss}->{ref} = $qseq;
                $hits->{$q_miss}->{hit} = $hseq;
                ##print "MISMATCH: ${query_name}: ${q_miss} ${qseq} ${hseq}\n";
            }

            if (scalar(@template_gap) > 0) {
                if (scalar(@template_gap) > 0) {
                    for my $ind (@template_gap) {
                        my $q = substr($q_seq, $ind - 2, 6);
                        my $h = substr($h_seq, $ind - 2, 6);
                        $hits->{$ind}->{type} = 'ins';
                        $hits->{$ind}->{ref} = $q;
                        ##print "INS (template gap): ${query_name}: ${ind} ${q} ${h}\n";
                    }
                }
            }

            if (scalar(@product_gap) > 0) {
                if (scalar(@product_gap) > 0) {
                    for my $ind (@product_gap) {
                        my $q = substr($q_seq, $ind - 2, 6);
                        my $h = substr($h_seq, $ind - 2, 6);
                        $hits->{$ind}->{type} = 'del';
                        $hits->{$ind}->{ref} = $h;
                        ##print "DEL (product gap): ${query_name}: ${ind} ${q} ${h}\n";
                    }
                }
            }

            push(@hitlst, $hits);
            my $hitlength = scalar(@hitlst);
            $indices->{$index} = \@hitlst;
        }                       ## End of each hit of a result.
    }                           ## End of each result.
    $f->close();
    ##my $stored = store($indices, $data_file);
    print "Stored indices to $data_file.\n";
    ## I think this is a good place to write some information about what was extracted.
    ## I need a reminder of this data structure.
    ## $data->{some index}->[ { 10 }->{type}='ins' ->{ref}='A'}, { 20 }->... ];
    ## $data->{same index}->[ { 12 }->{type}='mis' ] ...
    ## So hash of arrays of hashes of hashes.

    foreach my $idx (keys %{$indices}) {
        ## Each element of @hitlst is one read.
        my @hitlst = @{$indices->{$idx}};
        my $num = scalar(@hitlst);
        ## Each hit in hitlst is either an empty or an interesting hash.
        my $c = 0;
      POSLOOP: foreach my $datum (@hitlst) {
            $c++;
          HITLOOP: foreach my $pos (keys %{$datum}) {
                my $posref = $datum->{$pos};
                my $string = qq"";
                if ($posref->{type} eq 'ins') {
                    $string = qq"${idx}\t${c}\t${num}\t${pos}\tins\t$posref->{ref}\t \n";
                    print $string;
                } elsif ($posref->{type} eq 'del') {
                    $string = qq"${idx}\t${c}\t${num}\t${pos}\tdel\t \t$posref->{ref}\n";
                    print $string;
                } else {
                    $string = qq"${idx}\t${c}\t${num}\t${pos}\tmis\t$posref->{ref}\t$posref->{hit}\n";
                    print $string;
                }
                print $out $string;
            }  ## End looking at each position inside the list of positions.
        }      ## End looking at the list of positions in a given read.
    }          ## End looking at the reads in an index.
    print "\n";
    $out->close();
    return($indices);
}

sub Count_Parsed {
    my (%args) = @_;
    print "Starting counters.\n";
    my $input = $args{input};
    my $type = ref($input);
    my $data;
    if (!$type) {
        ## This is not a reference, and therefore a filename.
        $data = retrieve($input);
    } else {
        $data = $input;
    }
    use Data::Dumper;
    print Dumper $data;
    my @count_by_hits = ();
    foreach my $i (keys %{$data}) {
        my @hitlst = @{$data->{$i}};
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
