use strict;
use Data::Dumper;
$| = 1;

my $referenceFile = '/data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa';

# my $file = 'pileup21_cram';
my $file = 'pileup21_original_primary';
my $outfile = $file . '.SNPs';

print "Read reference $referenceFile\n";
my $reference_href = readFASTA($referenceFile, 1);

print "\nProcessing $file..\n";

my %haveNonRefAlleles;
open(F, '<', $file) or die "Cannot open $file";
open(FOUT, '>', $outfile) or die "Cannot open $outfile";
my %printed;
while(<F>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	
	my @fields = split(/\t/, $line);
	
	die "Reference contig $fields[0] not defined in $referenceFile" unless(exists $reference_href->{$fields[0]});
	
	my $position = $fields[1];
	my $coverage = $fields[3];
	my $alleles = $fields[4];	

	my %next_alleles_nextColumn;
	my @alleles;
	my $ignored_ref_skips = 0;
	for(my $i = 0; $i < length($alleles); $i++)
	{
		my $C = substr($alleles, $i, 1);
		if($C eq '^')
		{
			$i++;
			next;
		}
		next if($C eq '$');
		
		if(($C eq '+') or ($C eq '-'))
		{
			my $j = 1;
			while(substr($alleles, $i + $j + 1, 1) =~ /^\d$/)
			{
				$j++;
			}
			my $lengthInsertion = substr($alleles, $i+1, $j);
			die unless($lengthInsertion =~ /^\d+$/);
			my $allele = substr($alleles, $i + 1 + $j, $lengthInsertion);
			if(substr($alleles, $i, 1) eq '+')
			{
				$alleles[$#alleles] .= $allele;
			}
			else
			{
				# push(@alleles, '-');
				my $deletionAllele = '-' x $lengthInsertion;
				die unless(length($deletionAllele) == $lengthInsertion);
				$next_alleles_nextColumn{$deletionAllele}++;
			}
			$i = ($i + 1 + $j + $lengthInsertion) - 1;
		}
		elsif(($C eq '>') or ($C eq '<') or ($C eq '*'))
		{
			if(($C eq '>') or ($C eq '<'))
			{
				$ignored_ref_skips++;
				# die Dumper("This is unexpected, position $position of BAM $BAM");
			}
			else
			{
				push(@alleles, '-');		
			}
		}
		else
		{
			push(@alleles, substr($alleles, $i, 1));
		}	
	}	
		
	@alleles = map {uc($_)} @alleles;
	
	my @offendingAlleles = grep {$_ !~ /^[ACGTN\-]+$/} @alleles; 
	
	if(scalar(@offendingAlleles))
	{
		die Dumper("Offending alleles", \@offendingAlleles, $position, $alleles, $line);
	}
	  
	my $totalAlleles = scalar(@alleles);
	
	unless(($totalAlleles + $ignored_ref_skips) == $coverage)
	{ 
		warn Dumper("Mismatching coverage $totalAlleles vs $ignored_ref_skips", \@offendingAlleles, $position, $alleles, $line);
	}
	
	my $refC = substr($reference_href->{$fields[0]}, $position-1, 1);
	foreach my $allele (@alleles)
	{
		next unless($allele =~ /^[ACGT]$/i);
		if($refC ne $allele)
		{
			my $info_key = join('_', $fields[0], $position, $allele);
			next if($printed{$info_key});
			print FOUT join("\t", $fields[0], $position, $allele), "\n";
			$printed{$info_key} = 1;
		}
	}
}	
close(F);
close(FOUT);

print "\nProduced file $outfile - total alleles: ", scalar(keys %printed), "\n";

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
