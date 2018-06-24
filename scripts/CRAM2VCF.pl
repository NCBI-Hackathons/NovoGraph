#!/usr/bin/env perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Aarti Jajoo (Baylor), Nancy Hansen (NIH), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes_CSHL/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

# Example command:
# perl CRAM2VCF.pl --CRAM /intermediate_files/combined_2.cram
#     --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa
#     --output VCF/graph_v2.vcf
#     --contigLengths /intermediate_files/postGlobalAlignment_readLengths_2


my $CRAM;
my $referenceFasta;
my $output;
my $bin_BAM2VCF = '../BAM2VCF/BAM2VCF';
my $contigLengths;
unless(-e $bin_BAM2VCF)
{
	die "BAM2VCF binary $bin_BAM2VCF not present - run 'make all' in the directory.";
}
 
 
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

my %expectedLengths;
if($contigLengths)
{
	open(L, '<', $contigLengths) or die "Cannot open --contigLengths $contigLengths";
	while(<L>)
	{
		my $l = $_; 
		chomp($l);
		my @f = split(/\t/, $l);
		$expectedLengths{$f[0]} = $f[1];
	}
	close(L);
}

my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $CRAM);

my @sequence_ids = $sam->seq_ids();

my $reference_href;

print "Reading $referenceFasta\n";
$reference_href = readFASTA($referenceFasta);
print "\t...done.\n";

my %targetPos_printAlignments;
open(OUT, ">", $output) or die "Cannot open $output";
print OUT qq(##fileformat=VCFv4.2
##fileDate=20161026
##source=BAM2VCF.pl
##reference=file://$referenceFasta), "\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";

my @referenceSequenceIDs = @sequence_ids;

my %alignments_starting_at_test;

my %alignments_per_referenceSequenceID;



my $total_alignments = 0;
mkdir('forVCF');
die unless(-e 'forVCF');

my $fn_cmds = $output . '_CRAM2VCF_commands.txt';
my $fn_cmds_cat = $fn_cmds . '.cat';
my $fn_gaps = $output . '_CRAM2VCF_gaps.txt';
open(CMDS, '>', $fn_cmds) or die "Cannot open $fn_cmds";
print CMDS 'ulimit -u 10000', "\n";
open(CMDSCAT, '>', $fn_cmds_cat) or die "Cannot open $fn_cmds_cat";
open(GAPSTRUCTURE, '>', $fn_gaps) or die "Cannot open $fn_gaps";
print CMDSCAT 'cat ';
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	print "Processing $referenceSequenceID .. \n";
	die "Sequence ID $referenceSequenceID not defined in $referenceFasta" unless(exists $reference_href->{$referenceSequenceID});

	my $l_ref_sequence = length($reference_href->{$referenceSequenceID});
	my @gap_structure;
	$#gap_structure = ($l_ref_sequence - 1);

	die unless(scalar(@gap_structure) == $l_ref_sequence);
	
	my $fn_for_BAM2VCF = $output . '.part_'. $referenceSequenceID;
	my $fn_for_BAM2VCF_SNPs = $output . '.part_'. $referenceSequenceID . '.SNPs';
	
	open(D, '>', $fn_for_BAM2VCF) or die "Cannot open $fn_for_BAM2VCF";
	open(D2, '>', $fn_for_BAM2VCF_SNPs) or die "Cannot open $fn_for_BAM2VCF_SNPs";
	
	print D $reference_href->{$referenceSequenceID}, "\n";
	my $n_alignments = 0;
	my %alignments_starting_at;
	if(not $testing)
	{
		print "\t\tSet sequence ID to $referenceSequenceID\n";
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
				my $expectedLength = $expectedLengths{$alignment->query->name};
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
				
			}
			die unless(defined $firstMatch);
			
			$ref = substr($ref, $firstMatch, $lastMatch - $firstMatch + 1);
			$query = substr($query, $firstMatch, $lastMatch - $firstMatch + 1);		
			
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
			die unless(length($ref) == length($query));
            
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
			
			print D join("\t", $ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos), "\n";
			
			$alignments_per_referenceSequenceID{$referenceSequenceID}[0]++;
			(my $query_nonGap = $query) =~ s/[\-_\*]//g;
			$alignments_per_referenceSequenceID{$referenceSequenceID}[1] += length($query_nonGap);
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
	
	$total_alignments += $n_alignments;
	print "Have loaded $n_alignments alignments -- $fn_for_BAM2VCF.\n";

	for(my $refPos = 0; $refPos <= $#gap_structure; $refPos++)
	{
		my $n_gaps = $gap_structure[$refPos];
		if($n_gaps)
		{
			print GAPSTRUCTURE join("\t", $referenceSequenceID, $refPos, $n_gaps), "\n";
		}
	}	
	my $cmd = qq($bin_BAM2VCF --input $fn_for_BAM2VCF --referenceSequenceID $referenceSequenceID &> VCF/output_${referenceSequenceID}.txt &);
	
	my $output_file = $fn_for_BAM2VCF . '.VCF';
	
	print CMDS $cmd, "\n";
	
	print CMDSCAT $output_file . ' ';

	print "Now we could be executing: $cmd (manual)\n";	

	if(exists $targetPos_printAlignments{$referenceSequenceID})
	{
		foreach my $targetPosData (@{$targetPos_printAlignments{$referenceSequenceID}})
		{
			my $targetPos = $targetPosData->[0];
			my $outputFn = 'temp/' . $referenceSequenceID . '_around_' . $targetPos;
			outputMSAInto($targetPos, \%alignments_starting_at, $outputFn);
		}
	}
}

close(OUT);

print CMDSCAT ' >> ' . $output;

my $fn_alignment_details = $output . '_alignmentsPerRefID';
open(ALIGNMENTDETAILS, '>', $fn_alignment_details) or die "Cannot open $fn_alignment_details";
print ALIGNMENTDETAILS join("\t", qw/referenceID alignments alignedBases/), "\n";
foreach my $referenceSequenceID (keys %alignments_per_referenceSequenceID)
{
	print ALIGNMENTDETAILS join("\t", $referenceSequenceID, @{$alignments_per_referenceSequenceID{$referenceSequenceID}}), "\n";
}
close(ALIGNMENTDETAILS);

print "\nTotal alignments: $total_alignments\n";
print "\nDetailed info on alignments and bases per reference ID: see file $fn_alignment_details\n";
print "\nGap structure is in $fn_gaps\n";
print "\nOK\n\n";

close(CMDS);
close(GAPSTRUCTURE);


print"\nNow launch launch_CRAM2VCF_C++.pl with the right parameters.\n\n";

sub readFASTA
{
	my $file = shift;	
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
