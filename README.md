# NAME

errrt - A Module and script to examine Reverse Transcriptase error rates.

# SYNOPSIS

    use errrt;

# DESCRIPTION

errrt is intended to help process and analyze sequencing data generated as per
the paper: "Accurate fidelity analysis of the reverse transcriptase by a modified
next-generation sequencing."  As a simple synopsis of the experiment, a short
template sequence is used to generate cDNA with primers that have a random
terminus.

The resulting DNA sequences are then used to make illumina compatible sequencing
libraries.  Upon sequencing and merging the forward and reverse reads, we expect
to get a relatively homogenous set of sequences which look generally like this:

<pre>
NNNNNNN{14 nt. of random oligo} Amplified template sequence
       ^^^^ I will call this 14 nt. the RT index from now on.
</pre>

The hope is that we will have >= 3 reads for many 14 nucleotide random oligos.
If that is indeed true, then if we observe a mismatch between 1 sequenced read
out of a group of >= 3 and the template, then, we can safely assume that the
error was generated during the creation of the sequencing libraries.  In
contrast, if we observe a mismatch for a single base across all 3 reads of a
group, then we may assume that the mismatch was generated during reverse
transcription.

This script is intended to read the merged input library, count the reads
observed per index, perform an alignment against the template, categorize each
class of mismatch (insertion, deletion, transition, transversion), and count how
many are observed which span groups >= 3 reads vs. ones which do not.

## Step 1: Merge the forward/reverse reads

I use the tool 'flash' for this.  For example, using the first 40,000 reads from
one sample.

```{bash flash}
/cbcb/sw/RedHat-7-x86_64/common/local/flash/1.2.11/bin/flash \
  -o merged --compress-prog=xz --compress-prog-args=-9e \
  <(less share/r1_test.fastq.gz) <(less share/r2_test.fastq.gz)
```

The above command generates the file 'merged.extendedFrags.fastq.xz'.

## Step 2: Extract the template sequences by RT index.

Here is the first merged read in the data.

<pre>
AAAGGTAGTGCTGAATTCGATCACGCTGAGGGTGATTTTGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTCCCCGTCAAGCTCTAAATC
                         xxxxxxxxxxxxxx^^^^^^^^ Beginning of the non-random
                         ^^^^                   primer region
                         14 nt. RT index
</pre>

This should find all reads containing the non-random region, extract the
14 nt. random 'index' primer, categorize them by 'index,' and quantify the
mismatches between the rest of the sequence and the relevant portion of the
template.

With that in mind, here is the entire template region:

<pre>
GTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTC
</pre>

The beginning and end of this template will not be used for quantification,
leaving the following:

<pre>
CTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC
</pre>

The function Extract_Reads() accepts the following arguments:

* input:  Arbitrarily compressed merged forward/reverse reads (defaulting to
  'merged.extendedFrags.fastq.xz').
* output: Filename to write the extracted reads (defaulting to the input prefix
  .fasta, compressed I think).
* terminal: The first n nucleotides of the template, where n is defined by the
  command-line argument 'length' (defaulting to 9).

This function is pretty simple, it just uses a regex to look for the 'terminal'
sequence in each read.  If it does not find that, it reverse complements the
read and tries again.  If that terminal sequence was observed, it writes an
entry to a fasta file where the sequence ID is the
{read ID}_{read direction}_{terminal sequence} and the sequence is the remainder
of the read.

Upon completion, it wries a tab delimited file showing how many times every
RT index was observed along.

A shortcut script for this is provided as 'extract.pl'.

## Step 3: Compare the extracted reads to the template

The function Run_Fasta() runs my favorite aligner, fasta36, in order to compare
the set of extracted reads to our template sequence.  This is almost certainly
not the most efficient way of handling this step; however it has a significant
advantage coming from the fact that I already had an implementation handy for
parsing the fasta36 output in fashion which was basically perfect for the
questions inherent in this experiment.

It takes the following arguments:

* input: The output from the previous step.
* library_name: template.fasta by default.
* library_seq: as printed up above by default.
* fasta_cmd: fasta36 by default.
* fasta_args: -b 1 -d 1, e.g. grab the first hit only by default.

## Step 4: Parse the fasta36 output

Here is where I copy/pasted an old function I wrote for comparing gene sequences
from various trypanosomes years ago.  It takes the following arguments:

* input: The output from the previous step.
* output: parsed.txt by default.

This just uses Bio::SearchIO to read the output from fasta36 and give back a
data structure 'indices' which is basically a blown-up version of the what was
created in Step 2, thus:

1. each key of indices is one of the RT indices.
2. The value of each key is a list comprising a series of hits, each hit has:
  * A key which is the position within the template sequence.  If this is blank,
    then the element of the list tells us that the read it came from was
    identical to the template.  Therefore, if the lengths of these lists do not
    match the information found in the tab delimited file above, then we know
    something went wrong, the test suite will be definition contain this.
  * Each keyed template position in turn contains a 'type': insertion, deletion,
    and mismatch; a 'ref': the base found in the template at that position; and
    'hit': in the case of mismatches this records the new sequence observed.

This data structure is stored on disk for future analysis.

## Step 5 Poke at the set of observed differences.

I have not yet implemented this because I am not sure what I want to do here
yet.  There are a few obvious things to do:

1. Filter that data structure looking for only the changes observed which are
   found in > n reads with an indentical RT index.  The relatively small number
   of reads expected in this category define the errors introduced by the
   reverse transcriptase.
2. Simply record a matrix where columns are the positions in the template and
   the rows are differences.  This will be interesting to me to just get a
   distribution of error rates over all.
3. Do #2, but for the #1 filtered data.  If done intelligently, this should
   provide a way to do some interesting statistics with the data.

# AUTHOR

Ashton Trey Belew <abelew@gmail.com>

# COPYRIGHT

Copyright 2019- Ashton Trey Belew

# LICENSE

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

# SEE ALSO
