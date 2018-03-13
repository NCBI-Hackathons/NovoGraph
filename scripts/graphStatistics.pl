#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   

$| = 1;

my $VCF;
my $referenceFasta;
my $output;
 
GetOptions (
	'VCF:s' => \$VCF, 
	'referenceFasta:s' => \$referenceFasta, 
);

die "Please specify --VCF" unless($VCF);
die "File --VCF $VCF not existing" unless(-e $VCF);
die "Please specify --referenceFasta" unless(-e $referenceFasta);
die "File --referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $sawHeaderLine = 0;
my $lastEntryLastReferencePos;
my $lastEntryChromosome;
my $variantPositions_total = 0;
my %variantPositionsPerChromosome;
my $variantAlleles_total = 0;
my %variantAllelesPerChromosome;

my %variantAllelesLengthDifference_perChromosome;
my %variantAlleles_lengthDiff0_basesDiff_perChromosome;

open(VCF, '<', $VCF) or die "Cannot open $VCF";
while(<VCF>)
{
	chomp;
	next unless($_);
	next if(substr($_, 0, 2) eq '##');
	if(substr($_, 0, 1) eq '#')
	{
		die if($sawHeaderLine);
		$sawHeaderLine = 1;
	}
	else
	{
		die unless($sawHeaderLine);
		my @fields = split(/\t/, $_);
		
		my $chromosome = $fields[0];
		my $position = $fields[1];
		my $referenceAllele = $fields[3];
		my $variantAlleles = $fields[4];
		
		my @variantAlleles = split(/,/, $variantAlleles);
		if((defined $lastEntryLastReferencePos) and ($chromosome eq $lastEntryChromosome))
		{
			warn "It seems that there are overlapping vasriants in $VCF line $. - this allele starting position $chromosome:$position, last allele until $lastEntryChromosome:$lastEntryLastReferencePos" unless($position > $lastEntryLastReferencePos);
		}	
		
		$variantPositions_total++;
		$variantPositionsPerChromosome{$chromosome}++;
		
		$variantAlleles_total += scalar(@variantAlleles);
		$variantAllelesPerChromosome{$chromosome} += scalar(@variantAlleles);
		
		foreach my $vA (@variantAlleles)
		{
			my $lD = length($referenceAllele) - length($vA);
			$variantAllelesLengthDifference_perChromosome{'ALL'}{$lD}++;
			$variantAllelesLengthDifference_perChromosome{$chromosome}{$lD}++;
			if($lD == 0)
			{
				my $basesDifferent = 0;
				for(my $pI = 0; $pI < length($referenceAllele); $pI++)
				{
					if(substr($referenceAllele, $pI, 1) ne substr($vA, $pI, 1))
					{
						$basesDifferent++;
					}	
				}
				$variantAlleles_lengthDiff0_basesDiff_perChromosome{'ALL'}{$basesDifferent}++;
				$variantAlleles_lengthDiff0_basesDiff_perChromosome{$chromosome}{$basesDifferent}++;
			}
		}
		
		$lastEntryLastReferencePos = $position + length($referenceAllele) - 1;
		$lastEntryChromosome = $chromosome;
	}
}
close(VCF); 

my $fn_out = $VCF . '.statistics';

open(OUT, '>', $fn_out.'.lengthDiff') or die "Cannot open $fn_out.lengthDiff";
foreach my $chromosome (keys %variantAllelesLengthDifference_perChromosome)
{
	foreach my $value (keys %{$variantAllelesLengthDifference_perChromosome{$chromosome}})
	{
		print OUT join("\t", $chromosome, $value, $variantAllelesLengthDifference_perChromosome{$chromosome}{$value}), "\n";
	}
}
close(OUT);

open(OUT, '>', $fn_out.'.lengthDiff0_basesDifferent') or die "Cannot open $fn_out.lengthDiff0_basesDifferent";
foreach my $chromosome (keys %variantAlleles_lengthDiff0_basesDiff_perChromosome)
{
	foreach my $value (keys %{$variantAlleles_lengthDiff0_basesDiff_perChromosome{$chromosome}})
	{
		print OUT join("\t", $chromosome, $value, $variantAlleles_lengthDiff0_basesDiff_perChromosome{$chromosome}{$value}), "\n";
	}
}
close(OUT);

print "\nSummary of graph $VCF\n";
print "\tSize of the graph file: ", sprintf("%.2f", (-s $VCF)/(1024**3)), "Gb\n";
print "\tNumber of variant positions: $variantPositions_total\n";
print "\tNumber of variant allleles : $variantAlleles_total\n";
print "\nFurther statistics are in $fn_out.lengthDiff and $fn_out.lengthDiff0_basesDifferent\n\n";
