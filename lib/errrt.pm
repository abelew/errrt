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
use IPC::Open3;
use POSIX qw"ceil";
use Storable qw"store retrieve nstore";

use Data::Dumper;

use strict;
use 5.008_005;
our $VERSION = '0.01';

=head1 NAME

errrt - A few methods for boiling down RT error rate experiments.

=head1 SYNOPSIS

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

=head2 Methods

=over 4

=item C<Merge_Reads>

Invoke flash to merge the read1/read2 pairs into single reads.

=cut
sub Merge_Reads {
    my (%args) = @_;
    my $out = FileHandle->new(">flash.out");
    ##$|=1;
    ##my $pid = open3(\*IN,\*OUT, \*ERR , '/bin/bash');
    my $runner = qq"bash -c '$args{flash} \\
  -o $args{output} \\
  $args{flash_args} \\
  <(less $args{read1}) \\
  <(less $args{read2})'";
    ##print IN "${runner}\n";
    ##my $result = <OUT>;
    ##print $out $result;
    ##$out->close();
    my $handle = IO::Handle->new;
    print "Starting $runner\n";
    open($handle, "$runner |");
    while (my $line = <$handle>) {
        print $out "$line\n";
    }
    close($handle);
    $out->close();
}

=item C<Extract_Reads>

Look for a prefix sequence at the beginning of every read to identify reads suitable
for future analysis.  Write them to a new fasta file.

=cut
sub Extract_Reads {
    my (%args) = @_;
    my $final_seq = $args{final};
    my $inputted = FileHandle->new("less $args{input} |");
    my $in = Bio::SeqIO->new(-fh => $inputted, -format => 'Fastq');
    my $hist = FileHandle->new(">idx_count.txt");
    my $fasta = FileHandle->new(">$args{output}");
    my $ident_file = basename($args{output}, ('.fasta'));
    my $ident = FileHandle->new(">${ident_file}_identical_reads.txt");
    print $ident "read_num\tdirection\tindex\n";
    my $data_file = "$args{output}.pdata";

    my $indices = {};
    my $count = 0;
    my $identical = 0;
    my $non_identical = 0;
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
        unless ($sequence =~ m/$args{terminal}/) {
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
            $success++;
            if ($data =~ m/$final_seq/) {
                $identical++;
                print $ident "${count}\t${direction}\t${index}\n";
            } else {
                $non_identical++;
                print $fasta ">${count}_${direction}_${index}
$data
";
            }
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
    $ident->close();
    print "${success} indexed reads; identical: ${identical}, non-identical: ${non_identical}\n";
    store($indices, $data_file);
    return($indices)
}

=head2 C<Run_Fasta>

Invoke the fasta36 suite of programs with some reasonable options for learning the errors
from a fasta file of suitable reads.

=cut
sub Run_Fasta {
    my (%args) = @_;
    my $library = FileHandle->new(">$args{library_name}");
    print $library ">template_sequence
$args{library_seq}
";
    my $output = $args{output};
    if (!$output) {
        $output = basename($args{input}, ('.fasta'));
        $output = qq"${output}.txt.xz";
    }
    my $error = qq"${output}.err";

    my $runner = qq"$args{fasta_cmd} $args{fasta_args} $args{input} $args{library_name} 2>${error} |\
   xz -9e > ${output}";
    my $handle = IO::Handle->new;
    print "Starting $runner\n";
    open($handle, "$runner |");
    while (my $line = <$handle>) {
        print "$line\n";
    }
    close($handle);
    return($output);
}

=head2 C<Parse_Fasta>

Read the output from fasta36 and extract a data structure of mismatches, insertions, and deletions.

=cut
sub Parse_Fasta {
    my (%args) = @_;
    my $template_string = $args{template};
    my @template = split(//, $template_string);
    my $f = FileHandle->new("less $args{input} |");
    my $output = $args{output};
    my $out_base = basename($output, ('.txt'));
    my $out = FileHandle->new(">${output}");
    my $numbers = FileHandle->new(">${out_base}_numbers.txt");
    my $bads = FileHandle->new(">${out_base}_bads.txt");
    ## Write a header for the output tsv.
    print $out "readid\tindex\tdirection\tposition\ttype\treference\thit\n";

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
        ##my $score = $hit->raw_score();
        ##my $sig = $hit->significance();
        ##my $ident = $hit->frac_identical();
        my $hstrand = $hit->strand('hit');
        my $qstrand = $hit->strand('query');
        my ($readid, $direction, $index) = split(/_/, $query_name);
        if ($hstrand ne $qstrand) {
            ## I found a couple of reads which are messed up and reverse complementing
            ## I am not completely certain this will make them go away, but I think it will.
            print "SKIPPING: ${query_name} hit strand: $hstrand query strand: $qstrand\n";
            next SEARCHLOOP;
        }
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
        ## IMPORTANT POINT OF CONFUSION:
        ## The hit is the template sequence!
        ## The query is my provided sequence.
        my @template_mismatches = $hsp->seq_inds('hit', 'mismatch');
        my @product_mismatches = $hsp->seq_inds('query', 'mismatch');
        my @template_gaps = $hsp->seq_inds('hit', 'gap');
        my @product_gaps = $hsp->seq_inds('query', 'gap');
        my $num_mis = scalar(@template_mismatches);
        my $tem_gap = scalar(@template_gaps);
        my $pro_gap = scalar(@product_gaps);
        my $num_muts = $num_mis + $tem_gap + $pro_gap;
        print $numbers "${query_name}\t${num_muts}\n";
        my $t_seq = $hsp->hit_string;
        my @tseq = split(//, $t_seq);
        my $p_seq = $hsp->query_string;
        my @pseq = split(//, $p_seq);
        my $t_raw = $t_seq;
        $t_raw =~ s/\-//g;
        $t_raw =~ s/^\s+//g;
        my @traw = split(//, $t_raw);
        my $p_raw = $p_seq;
        $p_raw =~ s/\-//g;
        $p_raw =~ s/^\s+//g;
        my @praw = split(//, $p_raw);
        ## I am increasingly confident that we cannot trust indels at the end of alignments.
        ## No matter how I mess with the gap score/penalties, I get dumb alignments at the end.
        ## The tests which look at the end of the aligned sequence attempt to address this problem
        ## but I think a more aggressive approach will be required.  I think I will have to check
        ## the position of the indel and just arbitrarily drop it if it is >= 200.

        ## FIXME: I think there is an inconsistency in how I calculate the
        ## product nucleotide.  In the case of insertions and mismatches I use
        ## @praw, but in the case of deletions I use @pseq.  I would also like
        ## to figure out a way to consistently understand why/when removing the
        ## gap characters is necessary.  Presumably this is related to the weird
        ## requirement of adding $i to $index_position.  One would assume that
        ## we should use the raw data without the addition, but that does not
        ## result in sensible results.
        if (scalar(@template_gaps) > 0) {
            for my $i (0 .. $#template_gaps) {
                my $inc = $i + 1;
                my $total = scalar(@template_gaps);
                my $t_pos = $template_gaps[$i];
                my $index_position = $t_pos;
                ## I think the next two lines are a particularly nasty hack, and
                ## I do not understand why they seem to be necessary to get
                ## correct results.
                my $product_nt = $praw[$index_position + $i];
                my $template_nt = $tseq[$index_position + $i];
                my $template_length = scalar(@traw);
                if ($index_position >= $template_length) {
                    print "We passed the end of the alignment.\n";
                    next SEARCHLOOP;
                } elsif ($index_position >= 200) {
                    next SEARCHLOOP;
                }
                ##print "INS (template gap): ${index} ${index_position} $template_nt <$product_nt>\n"; #
                $hits->{$index_position}->{type} = 'ins';
                $hits->{$index_position}->{to} = $product_nt;
                print $out qq"${readid}\t${index}\t${direction}\t${index_position}\tins\t \t${product_nt}\n";
            }
        }
        if (scalar(@product_gaps) > 0) {
            for my $i (0 .. $#product_gaps) {
                my $p_pos = $product_gaps[$i];
                my $inc = $i + 1;
                my $total = scalar(@product_gaps);
                my $index_position = $p_pos;
                my $product_nt = $pseq[$index_position + $i];
                my $template_nt = $tseq[$index_position + $i];
                ## Now check that the alignment did not end.
                my $product_length = scalar(@pseq);
                if ($index_position >= $product_length) {
                    print "We passed the end of the alignment.\n";
                    next SEARCHLOOP;
                } elsif ($index_position >= 200) {
                    next SEARCHLOOP;
                }
                ##print "DEL (product gap): ${index} ${index_position} <${template_nt}> $product_nt\n";
                $hits->{$index_position}->{type} = 'del';
                $hits->{$index_position}->{from} = $template_nt;
                print $out qq"${readid}\t${index}\t${direction}\t${index_position}\tdel\t${template_nt}\t \n";
            }
        }
        if (scalar(@template_mismatches) > 0) {
            ## I expect the length of query_mismatches and hit_mismatches to be identical.
            for my $i (0 .. $#template_mismatches) {
                ## Get the index position of the mismatch.
                ## This can be confusing because of indels.
                ## I want them 0 indexed, but fasta36 gives them as 1 indexed, so subtract 1.
                my $t_pos = $template_mismatches[$i] - 1;
                my $p_pos = $product_mismatches[$i] - 1;
                ## $template_nt appears to be the correct one.
                my $template_nt = $template[$t_pos];
                my $inc = $i + 1;
                my $total = scalar(@template_mismatches);
                ## The product nucleotide is correct, except if there are indels
                ## before the mismatches...  In that case, shenanigans occur,
                ## and I absolutely cannot seem to figure out why.
                my $product_nt = $praw[$p_pos];
                my $index_position = $t_pos;
                $hits->{$index_position}->{type} = 'mis';
                $hits->{$index_position}->{from} = $template_nt;
                $hits->{$index_position}->{to} = $product_nt;
                ##print qq"MIS ${index}\t${inc}\t${total}\t${index_position}\tmis\t${template_nt}\t${product_nt}\n";
                print $out qq"${readid}\t${index}\t${direction}\t${index_position}\tmis\t${template_nt}\t${product_nt}\n";
            }
        }
        push(@hitlst, $hits);
        my $hitlength = scalar(@hitlst);
        $indices->{$index} = \@hitlst;
    } ## End of each hit of a result.
    $f->close();
    my $stored = store($indices, $data_file);
    print "Stored indices to $data_file.\n";
    ## I think this is a good place to write some information about what was extracted.
    ## I need a reminder of this data structure.
    ## $data->{some index}->[ { 10 }->{type}='ins' ->{ref}='A'}, { 20 }->... ];
    ## $data->{same index}->[ { 12 }->{type}='mis' ] ...
    ## So hash of arrays of hashes of hashes.
    $out->close();
    $numbers->close();
    $bads->close();
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


=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=cut
