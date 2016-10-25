#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
$| = 1;

my $referenceFasta;
my $BAM;

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

unless(-e 'forVG')
{
	mkdir('forVG') or die "Cannot mkdir forVG directory";
}

my $reference_href = readFASTA($referenceFasta);
my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my @sequence_ids = $sam->seq_ids();
foreach my $referenceSequenceID (@sequence_ids)
{	
	print "Processing $referenceSequenceID", "\n";
	
	next unless(defined $reference_href->{$referenceSequenceID});

	my @alignments = $sam->get_features_by_location(-seq_id => $referenceSequenceID);
	
	print "Region $referenceSequenceID, have ", scalar(@alignments), " alignments.\n";

}

sub readFASTA
{
	my $file = shift;	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000000) == 0)
		{
		# 	print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			$currentSequence = substr($line, 1);
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}
