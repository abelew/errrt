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

## Step 2: Note the template sequences.

Here is the first merged read in the data.

<pre>
AAAGGTAGTGCTGAATTCGATCACGCTGAGGGTGATTTTGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAATAAGATATCATCGGCTTTCCCCGTCAAGCTCTAAATC
                         xxxxxxxxxxxxxx^^^^^^^^ Beginning of the non-random
                         ^^^^                   primer region
                         14 nt random region
</pre>

This script should find all reads containing the non-random region, extract the
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


# AUTHOR

Ashton Trey Belew <abelew@gmail.com>

# COPYRIGHT

Copyright 2019- Ashton Trey Belew

# LICENSE

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

# SEE ALSO
