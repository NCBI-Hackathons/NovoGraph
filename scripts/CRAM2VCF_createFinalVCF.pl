#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

my $CRAM;
my $referenceFasta;
my $output;


 
GetOptions (
	'CRAM:s' => \$CRAM, 
	'referenceFasta:s' => \$referenceFasta, 
	'output:s' => \$output,
);
	
die "Please specify --CRAM" unless($CRAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --output" unless($output);

die "--CRAM $CRAM not existing" unless(-e $CRAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $CRAM);

my @sequence_ids = $sam->seq_ids();

open(OUT, ">", $output) or die "Cannot open $output";
print OUT qq(##fileformat=VCFv4.2
##fileDate=20161026
##source=BAM2VCF.pl
##reference=file://$referenceFasta), "\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";
#foreach my $referenceSequenceID (@sequence_ids)
my @referenceSequenceIDs = @sequence_ids;
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	my $fn_for_BAM2VCF = $output . '.part_'. $referenceSequenceID;
	my $fn_VCF = $fn_for_BAM2VCF . '.VCF';
	
	die "File $fn_for_BAM2VCF not present? Have you run CRAM2VCF.pl?" unless(-e $fn_for_BAM2VCF);
	unless(-e $fn_VCF)
	{
		warn "File $fn_VCF not existing - skip, but generate big VCF anyway.";
		next;
	}
	
	my $last_update_time_inputforVCF = (stat($fn_for_BAM2VCF))[9];
	my $last_update_time_VCF = (stat($fn_VCF))[9];	
	if($last_update_time_VCF < $last_update_time_inputforVCF)
	{
		warn "File $fn_VCF is older than $fn_for_BAM2VCF - skip, but generate big VCF anyway.";
		next;	
	}
	
	open(VCF, '<', $fn_VCF) or die "Cannot open $fn_VCF";
	while(<VCF>)
	{
		print OUT $_;
	}
	close(VCF);
	
}	

close(OUT);

print "\n\nGenerated file $output";
