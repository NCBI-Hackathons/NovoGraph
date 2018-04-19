#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max min all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

# example command: perl CRAM2MSA.pl --CRAM /data/projects/phillippy/projects/hackathon/intermediate_files/combined.cram --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa

my $CRAM;
my $referenceFasta;

my %targetPos_printAlignments;
while(<DATA>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @f = split(/\s+/, $line);
	die "Weird line in DATA (1): $line" unless(defined $f[1]);
	die "Weird line in DATA (2): $line" unless(defined $f[2]);
	die "Weird line in DATA (3): $line" unless(defined $f[3]);
	push(@{$targetPos_printAlignments{$f[1]}}, [@f[2 .. 3]]);
}	

GetOptions (
	'CRAM:s' => \$CRAM, 
	'referenceFasta:s' => \$referenceFasta, 
);

die "Please specify --CRAM" unless($CRAM);
die "Please specify --referenceFasta" unless($referenceFasta);

die "--CRAM $CRAM not existing" unless(-e $CRAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $CRAM);

my @sequence_ids = $sam->seq_ids();

print "Reading $referenceFasta\n";
my $reference_href = readFASTA($referenceFasta);
print "\t...done.\n";

my @referenceSequenceIDs = @sequence_ids;

my %alignments_per_referenceSequenceID;
my $total_alignments = 0;
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	next unless(exists $targetPos_printAlignments{$referenceSequenceID});
	
	my $fn_gaps = 'temp/gaps_' . $referenceSequenceID . '.txt';
	open(GAPSTRUCTURE, '>', $fn_gaps) or die "Cannot open $fn_gaps";
	
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
	
	my $n_alignments = 0;
	my %alignments_starting_at;

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
			my $c_query = substr($query, $i, 1);
			
			if(!(($c_ref eq '-') or ($c_ref eq '*')))
			{
				$ref_pos++;
			}
			
			if(($c_ref ne '-') and ($c_ref ne '*') and ($c_query ne '-') and ($c_query ne '*'))
			{
				
			}
		}
		
		my $alignment_last_pos = $ref_pos - 1;
		my $alignment_info_aref = [$ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos];
		push(@{$alignments_starting_at{$alignment_start_pos}}, $alignment_info_aref);
			
		$alignments_per_referenceSequenceID{$referenceSequenceID}[0]++;
		(my $query_nonGap = $query) =~ s/[\-_\*]//g;
		$alignments_per_referenceSequenceID{$referenceSequenceID}[1] += length($query_nonGap);
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
				
	$total_alignments += $n_alignments;
	print "Have loaded $n_alignments alignments\n";

	for(my $refPos = 0; $refPos <= $#gap_structure; $refPos++)
	{
		my $n_gaps = $gap_structure[$refPos];
		if($n_gaps)
		{
			print GAPSTRUCTURE join("\t", $referenceSequenceID, $refPos, $n_gaps), "\n";
		}
	}
	
	close(GAPSTRUCTURE);
	
	if(exists $targetPos_printAlignments{$referenceSequenceID})
	{
		foreach my $targetPosData (@{$targetPos_printAlignments{$referenceSequenceID}})
		{
			my $targetPos = $targetPosData->[0];
			my $outputFn = 'temp/' . $referenceSequenceID . '_around_' . $targetPos;
			outputMSAInto($targetPos, \%alignments_starting_at, $outputFn);
			# printHaplotypesAroundPosition($targetPos, \%alignments_starting_at);	
		}
	}
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

sub outputMSAInto
{
	my $posI = shift;
	my $alignments_starting_at_href = shift;
	my $outputFn = shift;
	
	open(OUT, '>', $outputFn) or die "Cannot open $outputFn";
	
	my %extractedSequences;
	my @aroundPositions = grep {$_ >= 0} (($posI - 2) .. ($posI + 2));
	my $minPos = min(@aroundPositions);
	my $maxPos = max(@aroundPositions);
	
	foreach my $startPos (keys %{$alignments_starting_at_href})
	{
		my $alignmentI = -1;
		my $n_alignments = scalar(@{$alignments_starting_at_href->{$startPos}});
		foreach my $alignment_aref (@{$alignments_starting_at_href->{$startPos}})
		{
			$alignmentI++;
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
				if(($interestingPos >= $startPos) and ($interestingPos <= $stopPos))
				{
					$interesting = 1;
				}
			}
			
			my $ref = $alignment_aref->[0];
			my $query = $alignment_aref->[1];
			
			if($interesting)
			{
				my %gt_per_position;
				my %gt_ref_per_position;
				my $ref_pos = $startPos - 1;
				my $running_allele = '';
				my $running_ref_allele = '';
				for(my $i = 0; $i < length($ref); $i++)
				{
					my $c_ref = substr($ref, $i, 1);
					my $c_query = substr($query, $i, 1);
					
					if(($c_ref eq '-') or ($c_ref eq '*'))
					{
						$running_allele .= $c_query;
						$running_ref_allele .= $c_ref;
					}
					else
					{
						if($running_allele)
						{
							$gt_per_position{$ref_pos} = $running_allele;
							$gt_ref_per_position{$ref_pos} = $running_ref_allele;
						}
						
						$running_allele = $c_query;
						$running_ref_allele = $c_ref;
						$ref_pos++;						
					}
				}
				if($running_allele)
				{
					$gt_per_position{$ref_pos} = $running_allele;
					$gt_ref_per_position{$ref_pos} = $running_ref_allele;
				}
						
				# print "Positions $alignment_aref->[2] \n";
				foreach my $interestingPos (@aroundPositions)
				{
					if(exists $gt_per_position{$interestingPos})
					{
						my $sequence_title = $alignment_aref->[2] . '_' . $minPos . ':' . $maxPos . '_' . $alignmentI . '_of_' . $n_alignments;
						$extractedSequences{$sequence_title}{$interestingPos} = $gt_per_position{$interestingPos};
						
						my $sequence_title_ref = 'ref' . '_' . $minPos . ':' . $maxPos;
						unless(exists $extractedSequences{$sequence_title_ref}{$interestingPos})
						{
							$extractedSequences{$sequence_title_ref}{$interestingPos} = $gt_ref_per_position{$interestingPos};						
						}
					}
				}
			}
		}
	}
	
	my %length_per_position;
	foreach my $sequence_title (keys %extractedSequences)
	{
		foreach my $interestingPos (@aroundPositions)
		{	
			if(exists $extractedSequences{$sequence_title}{$interestingPos})
			{
				my $l_allele = length($extractedSequences{$sequence_title}{$interestingPos});
				if(defined $length_per_position{$interestingPos})
				{
					unless($l_allele == $length_per_position{$interestingPos})
					{
						warn "Length mismatch $interestingPos $l_allele vs $length_per_position{$interestingPos}";
						if($l_allele > $length_per_position{$interestingPos})
						{
							$length_per_position{$interestingPos} = $l_allele;
						}
					}
				}
				else				
				{
					$length_per_position{$interestingPos} = $l_allele;
				}
				
			}
		}
	}
	
	my $seq_length;
	foreach my $sequence_title (keys %extractedSequences)
	{
		my @sequence_parts;
		foreach my $interestingPos (@aroundPositions)
		{	
			if(exists $extractedSequences{$sequence_title}{$interestingPos})
			{
				my $allele = $extractedSequences{$sequence_title}{$interestingPos};
				die unless(length($allele) <= $length_per_position{$interestingPos});
				$allele =~ s/\*/-/g;
				my $n_missing = $length_per_position{$interestingPos} - length($allele);
				my $star_sequence = ('*' x $n_missing);
				die unless(length($star_sequence) == $n_missing);
				$allele .= $star_sequence;
				push(@sequence_parts, $allele);								
			}
			else
			{
				my $n_sequence = ('N' x $length_per_position{$interestingPos});
				die unless(length($n_sequence) == $length_per_position{$interestingPos});
				push(@sequence_parts, $n_sequence);
			}
		} 
		my $sequence = join('', @sequence_parts);
		
		print OUT '>', $sequence_title, "\n", $sequence, "\n";
		$seq_length = length($sequence) unless(defined $seq_length);
		unless($seq_length == length($sequence))
		{
			warn "Length problem with $sequence_title - $seq_length vs " . length($sequence);
		}
	}
	
	
	
	
	close(OUT);
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
				if(($interestingPos >= $startPos) and ($interestingPos <= $stopPos))
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

__DATA__
34558 chr21 28891726   4729
51120 chr21 46060140  31837
50742 chr21 45802714 123647
242   chr21  5289128   1247
39875 chr21 34702674    946
39778 chr21 34700642    628
13832 chr21 10424743   2555
39509 chr21 34695843    820
39636 chr21 34697908    520
39788 chr21 34700829    822

