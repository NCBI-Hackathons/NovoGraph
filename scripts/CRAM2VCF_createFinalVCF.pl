#!/usr/bin/env perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;

## Usage: 
## CRAM2VCF_createFinalVCF.pl --CRAM <path to CRAM> 
##                            --referenceFasta <path to reference FASTA>  
##                            --prefix <path to VCF created by CRAM2VCF.pl>
##
## Example command:
## ./CRAM2VCF_createFinalVCF.pl --CRAM /intermediate_files/combined.cram
##                              --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa 
##                              --prefix VCF/graph.vcf


$| = 1;

my $CRAM;
my $referenceFasta;
my $output;
 
GetOptions (
	'CRAM:s' => \$CRAM, 
	'referenceFasta:s' => \$referenceFasta, 
	'prefix:s' => \$output,
);
	
die "Please specify --CRAM" unless($CRAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --prefix" unless($output);

die "--CRAM $CRAM not existing" unless(-e $CRAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $fn_CRAM_index = $CRAM . '.extract.index';
die "Expected file $fn_CRAM_index not found" unless(-e $fn_CRAM_index);

my @sequence_ids;
open(SAMHEADER, '<', $fn_CRAM_index) or die "Cannot open $fn_CRAM_index";
while(<SAMHEADER>)
{
	chomp;
	if($_ =~ /\@SQ\s+SN:(\S+)\s+LN:/)
	{
		push(@sequence_ids, $1);
	}
}
close(SAMHEADER);

my $VCF1 = $output . '.separated.VCF';
my $VCF2 = $output . '.overlapping.VCF';

open(OUT1, ">", $VCF1) or die "Cannot open $VCF1";
open(OUT2, ">", $VCF2) or die "Cannot open $VCF2";
print OUT1 qq(##fileformat=VCFv4.2
##fileDate=20161026
##source=CRAM2VCF.pl
##reference=file://$referenceFasta), "\n";
print OUT1 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";
print OUT2 qq(##fileformat=VCFv4.2
##fileDate=20161026
##source=CRAM2VCF.pl
##reference=file://$referenceFasta
##INFO=<ID=FROM,Number=A,Type=Integer,Description="Origin">), "\n";
print OUT2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";
#foreach my $referenceSequenceID (@sequence_ids)
my @referenceSequenceIDs = @sequence_ids;
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	my $fn_for_CRAM2VCF = $output . '.part_'. $referenceSequenceID;
	
	my $fn_VCF_1 = $fn_for_CRAM2VCF . '.separated.VCF';
	my $fn_VCF_2 = $fn_for_CRAM2VCF . '.overlapping.VCF';
	
	die "File $fn_for_CRAM2VCF not present? Have you run CRAM2VCF.pl?" unless(-e $fn_for_CRAM2VCF);
	
	my $fn_VCF_done = $fn_for_CRAM2VCF . '.done';
	unless(get_done($fn_VCF_done))
	{
		die "File $fn_VCF_done not indicating completion - skip.";		
		next;
	}
	
	unless((-e $fn_VCF_1) and (-e $fn_VCF_2))
	{
		warn "Either output file ( $fn_VCF_1 or $fn_VCF_2) not existing - skip, but generate big VCF anyway.";
		next;
	}
	
	my $last_update_time_inputforVCF = (stat($fn_for_CRAM2VCF))[9];
	my $last_update_time_VCF_1 = (stat($fn_VCF_1))[9];	
	my $last_update_time_VCF_2 = (stat($fn_VCF_2))[9];	
	if(($last_update_time_VCF_1 < $last_update_time_inputforVCF) or ($last_update_time_VCF_2 < $last_update_time_inputforVCF))
	{
		warn "File $fn_VCF_1 or $fn_VCF_2 is older than $fn_for_CRAM2VCF - skip, but generate big VCF anyway.";
		next;	
	}
	
	open(VCF, '<', $fn_VCF_1) or die "Cannot open $fn_VCF_1";
	while(<VCF>)
	{
		my @fields = split(/\t/, $_);
		unless(scalar(@fields) == 8)
		{
			warn "Weird number of fields in line $. of $fn_VCF_1 -- is $#fields + 1, but want 8";
			next;
		}
		print OUT1 $_;
	}
	close(VCF);
	
	open(VCF, '<', $fn_VCF_2) or die "Cannot open $fn_VCF_2";
	while(<VCF>)
	{
		my @fields = split(/\t/, $_);
		unless(scalar(@fields) == 8)
		{
			warn "Weird number of fields in line $. of $fn_VCF_2 -- is $#fields + 1, but want 8";
			next;
		}
		print OUT2 $_;
	}
	close(VCF);
	
}	

close(OUT1);
close(OUT2);

print "\n\nGenerated file $VCF1 and $VCF2\n\n(The .overlapping.VCF file is a straighforward representation of the input multiple sequence alignment in VCF format, with potentially overlapping variant alleles; the .separated.VCF file maintains a strict horizontal separation between variant alleles in the VCF. Which file you want to use depends on your downstream inference pipeline.)\n\n";

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
