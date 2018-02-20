use strict;
use Data::Dumper;
$| = 1;

my $referenceFile = '/data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa';

my $expectedSNPs = 'pileup21_cram.SNPs';
#my $expectedSNPs = 'VCF/uber_vcf.vcf.part_chr21.VCF.expectedSNPs';
#my $expectedSNPs = 'pileup21_original_primary.SNPs';
my $VCF = 'VCF/uber_vcf.vcf.part_chr21.VCF';

my $n_lookingFor_total = 0;
my %expected;
open(EXPECTED, '<', $expectedSNPs) or die "Cannot open $expectedSNPs";
while(<EXPECTED>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @fields = split(/\t/, $line);
	die unless(scalar(@fields) == 3);
	die if($expected{$fields[0]}{$fields[1]}{$fields[2]});
	$expected{$fields[0]}{$fields[1]}{$fields[2]} = 1;
	$n_lookingFor_total++;
}
close(EXPECTED);

print "\nNow scanning for $n_lookingFor_total alleles in $VCF\n";

my %positionsCoveredByINDEL;
open(VCF, '<', $VCF) or die "Cannot open $VCF";
while(<VCF>)
{ 
	my $line = $_;
	chomp($line);
	
	next if(substr($line, 0, 2) eq '##');
	
	if(substr($line, 0, 1) eq '#')
	{
		# my @header_fields = split(/\t/, $line);			
		# die Dumper("Wrong genotyping index", \@header_fields) unless($header_fields[$VCF_gt_field] eq 'NA12878');
		next;
	}
			
	my @fields = split(/\t/, $line);
	my $chr = $fields[0];
	my $position = $fields[1];
	my $refAllele = $fields[3];
	my @alternativeAlleles = split(/,/, $fields[4]);
	die "No alternative allele in line $. of $VCF? $line" unless(scalar(@alternativeAlleles));

	next unless($expected{$chr});

	my @alternativeAlleles_length = map {length($_)} @alternativeAlleles;
	my $alternativeAlleles_sameLength_all = 1;
	foreach my $aA (@alternativeAlleles)
	{
		$alternativeAlleles_sameLength_all = ($alternativeAlleles_sameLength_all and (length($aA) == length($refAllele)));
	}
	
	for(my $i = 0; $i < length($refAllele); $i++)
	{
		my $p = $position + $i;
		$positionsCoveredByINDEL{$chr}{$p} = 1;
	}
		
	foreach my $alternativeAllele (@alternativeAlleles)
	{
		my $minL = (length($refAllele) < length($alternativeAllele)) ? length($refAllele) : length($alternativeAllele);
		die unless($minL > 0);
		for(my $i = 0; $i < $minL; $i++)
		{
			my $refPos = $position + $i;
			my $refC = substr($refAllele, $i, 1);
			my $alleleC = substr($alternativeAllele, $i, 1);
			die unless($refC and $alleleC);
			if($expected{$chr}{$refPos}{$alleleC})
			{
				$expected{$chr}{$refPos}{$alleleC} = 0;
			}
		}
	}
}

my $found_alleles = 0;
my $notfound_alleles_potentially = 0;
my $notfound_alleles_forSure = 0;
foreach my $chr (sort keys %expected)
{
	foreach my $pos (sort keys %{$expected{$chr}})
	{
		foreach my $allele (sort keys %{$expected{$chr}{$pos}})
		{
			if($expected{$chr}{$pos}{$allele} == 0)
			{
				$found_alleles++;
			}
			else
			{
				if($positionsCoveredByINDEL{$chr}{$pos})
				{
					$notfound_alleles_potentially++;
				}
				else
				{
					#warn "Not found: $chr / $pos / $allele";
					#exit if($notfound_alleles_forSure >= 5);				
					$notfound_alleles_forSure++; 
				}
			}
		}
	}
}

print "\n\nFound $found_alleles / $n_lookingFor_total alleles from $expectedSNPs in $VCF -- $notfound_alleles_forSure certainly missing, $notfound_alleles_potentially not sure\n\n";

sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= uc($line);
		}
	}	
	close(F);
		
	return \%R;
}
