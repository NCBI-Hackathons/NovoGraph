#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;
use File::Path;

$| = 1;

my $BAM;
GetOptions (
	'BAM:s' => \$BAM, 
);

die "Please specify --BAM" unless($BAM);
die "--BAM $BAM not existing" unless(-e $BAM);

my $sam = Bio::DB::Sam->new(-bam => $BAM);

my %total_contig_IDs;
my %contigID_atLeastOneAlignment;

my @sequence_ids = $sam->seq_ids();
foreach my $referenceSequenceID (@sequence_ids)
{	
	my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);
	while(my $alignment = $alignment_iterator->next_seq)
	{
		my $readID = $alignment->query->name;
		$total_contig_IDs{$readID}++;
				
		# We wilfully ignore all the reference sequence existence / length conditions present in BAM2MAFFT.pl, as all
		# chr* entries should be there and sufficiently long
		
		next unless($referenceSequenceID =~ /chr[XY\d]+/);
		$contigID_atLeastOneAlignment{$readID}++;
	}
}

my $n_ignored = scalar(keys %total_contig_IDs) - scalar(keys %contigID_atLeastOneAlignment);

print "\nTotal contig entries:", scalar(keys %total_contig_IDs), "\n";
print "\t", "Skipped: ", $n_ignored, "\n";
print "\t", "Produce alignments: ", scalar(keys %contigID_atLeastOneAlignment), "\n";



