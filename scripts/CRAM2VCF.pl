#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

# Example command:
# To check correctness of INPUT for mafft:
# 	./BAM2VCF.pl --BAM /data/projects/phillippy/projects/hackathon/Graph_Genomes_CSHL/scripts/uber_sorted.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output uber_vcf.vcf
# ./BAM2VCF.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.fa --output uber_vcf.vcf
my $CRAM;
my $referenceFasta;
my $output;
my $bin_BAM2VCF = '../BAM2VCF/BAM2VCF';
my $contigLengths;
unless(-e $bin_BAM2VCF)
{
	die "BAM2VCF binary $bin_BAM2VCF not present - run 'make all' in the directory.";
}
#printHaplotypesAroundPosition(2, {
	# 0 => [
		# ['AAAC', 'TTTT', 'r0', 0, 3],
	# ],
	# 2 => [
		# ['ACGT', 'TTTT', 'r1', 2, 5],
		# ['AC-GT', 'TTATT', 'r2', 2, 5],	
		# ['ACGT', 'AA_A', 'r3', 2, 5],		
	# ],
	#2 => [
		# ['-CGT', 'TTTT', 'r4', 3, 5],	
	#	['ACG-', 'TTTT', 'r5', 2, 4],	
	#],
	#0 => [
	#	 ['AAAC', 'TTT-', 'r7', 0, 3],
	#],	
#});
 

GetOptions (
	'CRAM:s' => \$CRAM, 
	'referenceFasta:s' => \$referenceFasta, 
	'output:s' => \$output,
	'contigLengths:s' => \$contigLengths, 	
);
	
die "Please specify --CRAM" unless($CRAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --output" unless($output);

die "--CRAM $CRAM not existing" unless(-e $CRAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my %expectedLengts;
if($contigLengths)
{
	open(L, '<', $contigLengths) or die "Cannot open --contigLengths $contigLengths";
	while(<L>)
	{
		my $l = $_; 
		chomp($l);
		my @f = split(/\t/, $l);
		$expectedLengts{$f[0]} = $f[1];
	}
	close(L);
}

my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $CRAM);

my @sequence_ids = $sam->seq_ids();

my $testing = 0;

my $reference_href;
if(not $testing)
{
	print "Reading $referenceFasta\n";
	$reference_href = readFASTA($referenceFasta);
	print "\t...done.\n";
}

open(OUT, ">", $output) or die "Cannpt open $output";
print OUT qq(##fileformat=VCFv4.2
##fileDate=20161026
##source=BAM2VCF.pl
##reference=file://$referenceFasta), "\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";
#foreach my $referenceSequenceID (@sequence_ids)
my @referenceSequenceIDs = @sequence_ids;
if($CRAM =~ /TRY2/)
{
	# @referenceSequenceIDs = ("chr1");
}
my %alignments_starting_at_test;
if($testing)
{
	# ?CGT-A--CGT ref
	# -CTT-A--CGT r0
	#      AAACGT r1
	#   GTAT      r2
	#
		 
	@referenceSequenceIDs = qw/test1/;
	$reference_href = {
		test1 => "ACGTACGT",
	};
	
	%alignments_starting_at_test = (
	 1 => [
		 ['CGT-A--CGT', 'CTT-A--CGT', 'r0', 1, 7],
	 ],
	 4 => [
		 ['A--CGT', 'AAACGT', 'r1', 4, 7],
	 ],	 
	 2 => [
		 ['GT-A', 'GTAT', 'r2', 2, 4],
	 ],		 
	);
}

mkdir('forVCF');
die unless(-e 'forVCF');
# @referenceSequenceIDs = qw/chr21/;
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	print "Processing $referenceSequenceID .. \n";
	die "Sequence ID $referenceSequenceID not defined in $referenceFasta" unless(exists $reference_href->{$referenceSequenceID});
	
	my $l_ref_sequence = length($reference_href->{$referenceSequenceID});
	my @gap_structure;
	$#gap_structure = ($l_ref_sequence - 1);
	print "Set...";
	#for(my $i = 0; $i < $l_ref_sequence; $i++)
	#{
	#	#$gap_structure[$i] = -1;
	#}
	print " .. done.\n";
	die unless(scalar(@gap_structure) == $l_ref_sequence);
	
	my $fn_for_BAM2VCF = $output . '.part_'. $referenceSequenceID;
	open(D, '>', $fn_for_BAM2VCF) or die "Cannot open $fn_for_BAM2VCF";
	print D $reference_href->{$referenceSequenceID}, "\n";
	my $n_alignments = 0;
	my %alignments_starting_at;
	if(not $testing)
	{
		my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);	
		while(my $alignment = $alignment_iterator->next_seq)
		{
			$n_alignments++;		

			my $alignment_start_pos = $alignment->start - 1;
			my ($ref,$matches,$query) = $alignment->padded_alignment;
			unless(length($ref) == length($query))
			{
				warn Dumper("Alignment length mismatch", length($ref), length($query), $alignment->query->name);
				next;
			}	
			
			if($contigLengths)
			{
				my $query_noGaps = $query;
				$query_noGaps =~ s/[\-_\*]//g;
				my $expectedLength = $expectedLengts{$alignment->query->name};
				die "No length for " . $alignment->query->name unless(defined $expectedLength);
				unless($expectedLength == length($query_noGaps))
				{
					die Dumper("Sequence length mismatch", $query_noGaps, length($query_noGaps), $expectedLength, length($alignment->query->dna), $alignment->query->name);
				}
			}

			
			
			my $ref_preAll = $ref;
			my $query_preAll = $query;
			
			my $firstMatch;
			my $lastMatch;
			for(my $i = 0; $i < length($ref); $i++)
			{
				if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*') and (substr($query, $i, 1) ne '-') and (substr($query, $i, 1) ne '*'))
				{
					$firstMatch = $i unless(defined $firstMatch);
					$lastMatch = $i;
				}
				
				my $isDeletion = (((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*')) and ((substr($query, $i, 1) eq '-') or (substr($query, $i, 1) eq '*')));
				my $isInsertion = (((substr($query, $i, 1) ne '-') and (substr($query, $i, 1) ne '*')) and ((substr($ref, $i, 1) eq '-') or (substr($ref, $i, 1) eq '*')));
				die if($isDeletion and $isInsertion);
				
				if($isDeletion)
				{
					# warn Dumper("Deletions before start?", $ref, $query, $alignment->query->name, $alignment->cigar_str);
				}
				if($isInsertion)
				{
					# warn Dumper("One of the insertions!", $ref, $query, $alignment->query->name, $alignment->cigar_str);
				}	
			}
			die unless(defined $firstMatch);
			
			$ref = substr($ref, $firstMatch, $lastMatch - $firstMatch + 1);
			$query = substr($query, $firstMatch, $lastMatch - $firstMatch + 1);		
			
			if($alignment->query->name eq 'Korean.gi|1078261939|gb|LPVO02001249.1|')
			{
				# warn Dumper($alignment->query->name, $ref_preAll, $query_preAll, $ref, $query);
			}
			my $gaps_left_side = 0;
			my $gaps_right_side = 0;
			
			for(my $i = 0; $i < length($ref); $i++)
			{
				if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*'))
				{
					last;
				}
				else
				{
					$gaps_left_side++;
				}
			}
			
			for(my $i = length($ref) - 1; $i >= 0; $i--)
			{
				if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*'))
				{
					last;
				}
				else
				{
					$gaps_right_side++;
				}
			}
			
			my $starstar_left_side = 0;
			my $starstar_right_side = 0;
			for(my $i = 0; $i < length($ref); $i++)
			{
				if((substr($ref, $i, 1) eq '*') and (substr($query, $i, 1) eq '*'))
				{
					$starstar_left_side++;
				}
				else
				{
					last;
				}
			}
			
			for(my $i = length($ref) - 1; $i >= 0; $i--)
			{
				if((substr($ref, $i, 1) eq '*') and (substr($query, $i, 1) eq '*'))
				{
					$starstar_right_side++;
				}
				else
				{
					last;
				}
			}
			
			
			die unless($gaps_left_side == 0);
			die unless($gaps_right_side == 0);
			
			$ref = substr($ref, $gaps_left_side, length($ref) - $gaps_left_side - $gaps_right_side);
			$query = substr($query, $gaps_left_side, length($query) - $gaps_left_side - $gaps_right_side);		
			#$ref = substr($ref, 0, length($ref) - $gaps_right_side);
			#$query = substr($query, 0, length($query) - $gaps_right_side);		
			die unless(length($ref) == length($query));
			
			if($alignment->query->name eq 'CHM13.gi|953910992|gb|LDOC03004332.1|')
			{
				#die Dumper($alignment_start_pos, $ref, $query);
			}			
		
			my $ref_pos = $alignment_start_pos - 1;
			my $running_gaps = 0;
			for(my $i = 0; $i < length($ref); $i++)
			{
				my $c_ref = substr($ref, $i, 1);
				
				if(!(($c_ref eq '-') or ($c_ref eq '*')))
				{
					$ref_pos++;
				}
			}
			
			my $alignment_last_pos = $ref_pos - 1;
			my $alignment_info_aref = [$ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos];
			push(@{$alignments_starting_at{$alignment_start_pos}}, $alignment_info_aref);
			
			print D join("\t", $ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos), "\n";
		}
	}
	else
	{
		%alignments_starting_at = %alignments_starting_at_test;
	}
	
	# reference gaps *before* the i-th reference character
	my $examine_gaps_n_alignment = 0;
	foreach my $alignment_start_pos (keys %alignments_starting_at)
	{
		foreach my $alignment (@{$alignments_starting_at{$alignment_start_pos}})
		{
			my $ref = $alignment->[0];
			
			my $start_pos = $alignment_start_pos - 1;
			my $ref_pos = $start_pos;
			my $running_gaps = 0;
			for(my $i = 0; $i < length($ref); $i++)
			{
				my $c_ref = substr($ref, $i, 1);
				
				if(($c_ref eq '-') or ($c_ref eq '*'))
				{
					$running_gaps++;
				}
				else
				{
					if(($ref_pos >= 5309528) and ($ref_pos <= 5309532))
					{
						# print "Gaps $ref_pos ". $alignment->[2] . ": " . $running_gaps . "\n";
					}	
					if($ref_pos != $start_pos)
					{
						if(not defined $gap_structure[$ref_pos])
						{
							$gap_structure[$ref_pos] = $running_gaps;
						}
						else
						{
							die "Gap structure mismatch at position $ref_pos - this is alignment $examine_gaps_n_alignment, have existing value $gap_structure[$ref_pos], want to set $running_gaps" unless($gap_structure[$ref_pos] == $running_gaps);
						}
					}
					$ref_pos++;
					$running_gaps = 0;
					$examine_gaps_n_alignment++;
				}
			}	
		}
	}
			
	close(D);
	
	print "Have loaded $n_alignments alignments -- $fn_for_BAM2VCF.\n";

	my $cmd = qq($bin_BAM2VCF --input $fn_for_BAM2VCF --referenceSequenceID $referenceSequenceID);
	
	print "Now executing: $cmd\n";	
	
	my $output_file = $fn_for_BAM2VCF . '.VCF';
	
	# todo
	#if(system($cmd))
	#{
	#	die "Command $cmd failed";
	#}
	
	unless(-e $output_file)
	{
		die "File $output_file not existing";
	}
	
	open(CHROUT, '<', $output_file) or die "Cannot open $output_file";
	while(<CHROUT>)
	{
		print OUT $_;
	}
}

close(OUT);

print "\n\nOK\n\n";


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
			$currentSequence =~ s/\s.+//;
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
}

sub printHaplotypesAroundPosition
{
	my $posI = shift;
	my $alignments_starting_at_href = shift;
	
	my @aroundPositions = grep {$_ >= 0} (($posI - 2) .. ($posI + 2));
	foreach my $startPos (keys %{$alignments_starting_at_href})
	{
		foreach my $alignment_aref (@{$alignments_starting_at_href->{$startPos}})
		{
			my $apparent_startPos = $alignment_aref->[3];
			if($apparent_startPos != $startPos)
			{
				warn Dumper("Apparent / startpos mismatch", $startPos, $apparent_startPos, $alignment_aref->[2]);
			}
			my $stopPos = $alignment_aref->[4];
			die unless($stopPos >= $startPos);
			my $interesting = 0;
			foreach my $interestingPos (@aroundPositions)
			{
				if(($interestingPos >= $startPos) and ($interesting <= $stopPos))
				{
					$interesting = 1;
				}
			}
			
			my $ref = $alignment_aref->[0];
			my $query = $alignment_aref->[1];
			
			if($interesting)
			{
				my %gt_per_position;
				my $ref_pos = $startPos - 1;
				my $running_allele = '';
				for(my $i = 0; $i < length($ref); $i++)
				{
					my $c_ref = substr($ref, $i, 1);
					my $c_query = substr($query, $i, 1);
					
					if(($c_ref eq '-') or ($c_ref eq '*'))
					{
						$running_allele .= $c_query;
					}
					else
					{
						if($running_allele)
						{
							$gt_per_position{$ref_pos} = $running_allele;
						}
						
						$running_allele = $c_query;
						$ref_pos++;						
					}
				}
				if($running_allele)
				{
					$gt_per_position{$ref_pos} = $running_allele;
				}
						
				print "Positions $alignment_aref->[2] \n";
				foreach my $interestingPos (@aroundPositions)
				{
					if(exists $gt_per_position{$interestingPos})
					{
						print "\t", $interestingPos, "\t", $gt_per_position{$interestingPos}, "\n";
					}
				}
			}
		}
	}
}

