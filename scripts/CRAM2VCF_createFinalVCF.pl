#!/usr/bin/env perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

## Usage: 
## CRAM2VCF_createFinalVCF.pl --CRAM <path to CRAM> 
##                            --referenceFasta <path to reference FASTA>  
##                            --output <path to VCF created by CRAM2VCF.pl>
##
## Example command:
##     ./CRAM2VCF_createFinalVCF.pl --CRAM /intermediate_files/combined.cram
##                                  --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa 
##                                  --output VCF/graph.vcf


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
##source=CRAM2VCF.pl
##reference=file://$referenceFasta), "\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";
#foreach my $referenceSequenceID (@sequence_ids)
my @referenceSequenceIDs = @sequence_ids;
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	my $fn_for_CRAM2VCF = $output . '.part_'. $referenceSequenceID;
	my $fn_VCF = $fn_for_CRAM2VCF . '.VCF';
	
	die "File $fn_for_CRAM2VCF not present? Have you run CRAM2VCF.pl?" unless(-e $fn_for_CRAM2VCF);
	
	my $fn_VCF_done = $fn_VCF . '.done';
	unless(get_done($fn_VCF_done))
	{
		die "File $fn_VCF_done not indicating completion - skip.";		
		next;
	}
	
	unless(-e $fn_VCF)
	{
		warn "File $fn_VCF not existing - skip, but generate big VCF anyway.";
		next;
	}
	
	my $last_update_time_inputforVCF = (stat($fn_for_CRAM2VCF))[9];
	my $last_update_time_VCF = (stat($fn_VCF))[9];	
	if($last_update_time_VCF < $last_update_time_inputforVCF)
	{
		warn "File $fn_VCF is older than $fn_for_CRAM2VCF - skip, but generate big VCF anyway.";
		next;	
	}
	
	open(VCF, '<', $fn_VCF) or die "Cannot open $fn_VCF";
	while(<VCF>)
	{
		my @fields = split(/\t/, $_);
		unless(scalar(@fields) == 8)
		{
			warn "Weird number of fields in line $. of $fn_VCF -- is $#fields + 1, but want 8";
			next;
		}
		print OUT $_;
	}
	close(VCF);
	
}	

close(OUT);

print "\n\nGenerated file $output\n\n";

sub get_done
{
	my $fn = shift;
	return 0 unless($fn);
	open(F, '<', $fn) or die "Cannot open $fn";
	my $done = <F>;
	chomp($done);
	$done = (length($done) > 0) ? substr($done, 0, 1) : 0;	
	close(F);
	return $done;
}
