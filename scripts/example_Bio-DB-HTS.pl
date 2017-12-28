#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::DB::HTS;


my $referenceFasta = '/data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa';
my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => 'testCase.cram');

my $alignment_iterator = $sam->features(-iterator => 1);	
while(my $alignment = $alignment_iterator->next_seq)
{
	my ($ref,$matches,$query) = $alignment->padded_alignment;
	unless(length($ref) == length($query))
	{
		warn Dumper("Alignment length mismatch", length($ref), length($query), $alignment->query->name);
		next;
	}	
}