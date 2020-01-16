#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use Bio::DB::HTS;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;
use File::Path;

## Usage: 
## countExpectedGlobalAlignments.pl --BAM <path to BAM output of FIND_GLOBAL_ALIGNMENTS.pl>
##
## Example command:
## 	./countExpectedGlobalAlignments.pl --BAM .../intermediate_files/forMAFFT.bam

$| = 1;

my $BAM;
GetOptions (
	'BAM:s' => \$BAM, 
);

die "Please specify --BAM" unless($BAM);
die "--BAM $BAM not existing" unless(-e $BAM);

my $sam = Bio::DB::HTS->new(-bam => $BAM);

my %total_contig_IDs;
my %contigID_atLeastOneAlignment;
my %skipped_alignments_from;
my %lengths_by_category;
my %lengths_by_readID;

my @sequence_ids = $sam->seq_ids();
my $nbases_used = 0;
my $nbases_skipped = 0;
foreach my $referenceSequenceID (@sequence_ids)
{	
	my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);
	while(my $alignment = $alignment_iterator->next_seq)
	{
		my $readID = $alignment->query->name;
		my $sequenceLength = length($alignment->query->dna);
		die "Duplicate contig ID $readID ?" if($total_contig_IDs{$readID});
		$total_contig_IDs{$readID}++;
				
		$lengths_by_readID{$readID} = $sequenceLength;
		# We wilfully ignore all the reference sequence existence / length conditions present in BAM2MAFFT.pl, as all
		# chr* entries should be there and sufficiently long
		
		unless($referenceSequenceID =~ /chr[XY\d]+/)
		{
			warn $referenceSequenceID;
			$skipped_alignments_from{$referenceSequenceID}{$readID}++;
			push(@{$lengths_by_category{'skipped'}}, $sequenceLength);
			$nbases_skipped += $sequenceLength;			
			next;
		}
		$contigID_atLeastOneAlignment{$readID}++;
		push(@{$lengths_by_category{'taken'}}, $sequenceLength);		
		$nbases_used += $sequenceLength;		
	}
}

open(SKIPPED, '>', '_skippedAlignmentsFrom') or die;


print "\nNon-chr contig sources (details in file _skippedAlignmentsFrom:\n";
my @nonChrSources = sort {scalar(keys %{$skipped_alignments_from{$b}}) <=> scalar(keys %{$skipped_alignments_from{$a}})} keys %skipped_alignments_from;
my $sourceI = 0;
foreach my $key (@nonChrSources)
{
	my @readIDs = keys %{$skipped_alignments_from{$key}};
	if($sourceI <= 10)
	{
		print "\t - $key ", scalar(@readIDs), " \n";
	}
	
	foreach my $readID (@readIDs)
	{
		print SKIPPED join("\t", $key, $readID), "\n";
	}
	
	$sourceI++;
}
my $n_ignored = scalar(keys %total_contig_IDs) - scalar(keys %contigID_atLeastOneAlignment);
close(SKIPPED);

print "\nRead length details are in file _readLength_bySkippedOrNot\n";

open(RLBYCAT, '>', '_readLength_bySkippedOrNot') or die;
foreach my $key (keys %lengths_by_category)
{
	foreach my $rL (@{$lengths_by_category{$key}})
	{
		print RLBYCAT join("\t", $key, $rL), "\n";
	}	
}
close(RLBYCAT);

print "\nTotal contig entries:", scalar(keys %total_contig_IDs), "\n";
print "\t", "Skipped: ", $n_ignored, " / ", sprintf("%.2f", $nbases_skipped/1e6), "Mb\n";
print "\t", "Produce alignments: ", scalar(keys %contigID_atLeastOneAlignment), " / ", sprintf("%.2f", $nbases_used/1e6), "Mb\n";



