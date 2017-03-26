#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::Sam;

$| = 1;

# Example command:
# To check correctness of INPUT for mafft:
# 	./BAM2VCF.pl --BAM /data/projects/phillippy/projects/hackathon/Graph_Genomes_CSHL/scripts/uber_sorted.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output uber_vcf.vcf
# ./BAM2VCF.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.fa --output uber_vcf.vcf
my $BAM;
my $referenceFasta;
my $output;
my $bin_BAM2VCF = '../BAM2VCF/BAM2VCF';
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
	'BAM:s' => \$BAM, 
	'referenceFasta:s' => \$referenceFasta, 
	'output:s' => \$output,
);
	
die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --output" unless($output);

die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

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
if($BAM =~ /TRY2/)
{
	@referenceSequenceIDs = ("chr1");
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
#@referenceSequenceIDs = qw/chr21/;
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
		
	print "Have loaded $n_alignments alignments.\n";
	
	close(D);
	
	my $cmd = qq($bin_BAM2VCF --input $fn_for_BAM2VCF --referenceSequenceID $referenceSequenceID);
	
	print "Now executing: $cmd\n";
	
	if(system($cmd))
	{
		die "Command $cmd failed";
	}
	
	next;
	
	# my $last_all_equal = 0;
	my @open_haplotypes = (['', 0, -1]);
	my $start_open_haplotypes = 0;
	my $opened_alignments = 0;
	POSI: for(my $posI = 0; $posI < length($reference_href->{$referenceSequenceID}); $posI++)
	{
		if((($posI % 10000) == 0) or (scalar(@open_haplotypes) > 100))
		{
			print $posI, ", open haplotypes: ", scalar(@open_haplotypes), "\n";
		}
		
		foreach my $haplotype (@open_haplotypes)
		{
			# consume all gaps "before" the current reference position
			if($haplotype->[1] != 0)
			{
				if($haplotype->[2] == (length($haplotype->[1][0]) - 1))
				{
					my $n_gaps = $gap_structure[$posI-1];
					$n_gaps = 0 if(not defined $n_gaps);
					my $gaps = '-' x $n_gaps;
					die unless(length($gaps) == $n_gaps);
					$haplotype->[0] .= $gaps;
				}
				else
				{
					if((substr($haplotype->[1][0], $haplotype->[2], 1) eq '-') or (substr($haplotype->[1][0], $haplotype->[2], 1) eq '*'))
					{
						warn "Position $haplotype->[2] is gap in one of our haplotypes!";
					}
					
					my $nextPos = $haplotype->[2]+1;
					my $additionalExtension = '';
					while(($nextPos < length($haplotype->[1][0])) and ((substr($haplotype->[1][0], $nextPos, 1) eq '-') or (substr($haplotype->[1][0], $nextPos, 1) eq '*')))
					{
						$additionalExtension .= substr($haplotype->[1][1], $nextPos, 1);
						$nextPos++;
					}
					my $consumedUntil = $nextPos - 1;
					
					$haplotype->[0] .= $additionalExtension;
					$haplotype->[2] = $consumedUntil;
				}
			}	
			else
			{	
				if($posI > 0)
				{
					my $n_gaps = $gap_structure[$posI-1];
					$n_gaps = 0 if(not defined $n_gaps);
					my $gaps = '-' x $n_gaps;
					die unless(length($gaps) == $n_gaps);
					$haplotype->[0] .= $gaps;
				}
			}
		}
		
		# check that all extensions - i.e. up to the current reference position - have the same length
		my $assembled_h_length;
		foreach my $openHaplotype (@open_haplotypes)
		{
			unless(defined $assembled_h_length)
			{  
				$assembled_h_length = length($openHaplotype->[0]);
			}
			
			unless($assembled_h_length == length($openHaplotype->[0]))
			{
				print Dumper("Initial II length mismatch", $posI, $assembled_h_length, [@gap_structure[($posI-3) .. ($posI+1)]]), "\n";
				foreach my $oH2 (@open_haplotypes)
				{
					print "\t", length($oH2->[0]), "\t", " consumed until: ", $oH2->[2], ", of length ",  (($oH2->[1] == 0) ? "REF" : "nonRef $oH2->[1][2] / length " . length($oH2->[1][0])), , "\n";
				}	
				printHaplotypesAroundPosition($posI, \%alignments_starting_at);
				die;
			}
			# die Dumper("Initial II length mismatch", $posI, $assembled_h_length, length($openHaplotype->[0]), \@open_haplotypes, "Initial II length mismatch") unless($assembled_h_length == length($openHaplotype->[0]));
		}
		
		# for(my $existingHaploI = 0; $existingHaploI <= $#open_haplotypes; $existingHaploI++)
		# {
			# print "\t", $existingHaploI, "\t", $open_haplotypes[$existingHaploI][0], " ", $open_haplotypes[$existingHaploI][1], " ", $open_haplotypes[$existingHaploI][2], "\n";
		# }

		# copy in new haplotypes
		my $refC = substr($reference_href->{$referenceSequenceID}, $posI, 1);
		my @new_haplotypes = (exists $alignments_starting_at{$posI}) ? @{$alignments_starting_at{$posI}} : ();
		my $open_haplotypes_maxI = $#open_haplotypes;
		foreach my $new_haplotype (@new_haplotypes)
		{
			# print "Enter new haplotype ", $new_haplotype->[2], "\n";
			if($open_haplotypes_maxI > -1)
			{				
				$opened_alignments++;
				
				for(my $existingHaploI = 0; $existingHaploI <= $open_haplotypes_maxI; $existingHaploI++)
				{
					my $new_haplotype_copy_this = [$open_haplotypes[$existingHaploI][0], $new_haplotype, -1];
					push(@open_haplotypes, $new_haplotype_copy_this);					
				}
				
				my $open_span = $posI - $start_open_haplotypes;
				my $start_reference_extraction = $start_open_haplotypes;
				my $stop_reference_extraction = $posI - 1;
				my $referenceExtraction;
				die Dumper("Weird", $start_reference_extraction, $stop_reference_extraction) unless($stop_reference_extraction >= $start_reference_extraction);
				for(my $refI = $start_reference_extraction; $refI <= $stop_reference_extraction; $refI++)
				{
					$referenceExtraction .= substr($reference_href->{$referenceSequenceID}, $refI, 1);
					my $n_gaps = $gap_structure[$refI];
					$n_gaps = 0 if (not defined $n_gaps);
					my $gaps = '-' x $n_gaps;
					die unless(length($gaps) == $n_gaps);
					$referenceExtraction .= $gaps; 
					
				}
				#my $new_haplotype_referenceSequence = [substr($reference_href->{$referenceSequenceID}, $start_open_haplotypes, $open_span), $new_haplotype, -1];
				my $new_haplotype_referenceSequence = [$referenceExtraction, $new_haplotype, -1];
				#my $missing = $assembled_h_length - $open_span;
				#die Dumper($posI, $start_open_haplotypes, $missing, $open_span, $assembled_h_length) unless($missing >= 0);
				#my $missingStr = '*' x $missing;
				#die unless(length($missingStr) == $missing);
				#$new_haplotype_referenceSequence->[0] .= $missingStr;
				push(@open_haplotypes, $new_haplotype_referenceSequence);					
			}
		}
		
		# we exit these haplotypes be making them ref / another haplotype
		$open_haplotypes_maxI = $#open_haplotypes;
		foreach my $haplotype (@open_haplotypes)
		{
			if($haplotype->[1])
			{
				if($haplotype->[2] == (length($haplotype->[1][0]) - 1))
				{
					# print "exit one\n";
					$haplotype->[1] = 0;
					$haplotype->[2] = -1;
					
					for(my $existingHaploI = 0; $existingHaploI <= $open_haplotypes_maxI; $existingHaploI++)
					{
						next if($open_haplotypes[$existingHaploI] == $haplotype);
						my $new_haplotype_copy_this = [$haplotype->[0], $open_haplotypes[$existingHaploI]->[1], $open_haplotypes[$existingHaploI]->[2]];
						push(@open_haplotypes, $new_haplotype_copy_this);					
					}
				
				}
			}
		}
		
		# print "\tLength ", $assembled_h_length, "\n";
				

		if(1 == 0)
		{
			print "Haplotype info:\n";
			for(my $existingHaploI = 0; $existingHaploI <= $#open_haplotypes; $existingHaploI++)
			{
				print "\t", $existingHaploI, "\n";
				print "\t\t", $open_haplotypes[$existingHaploI][0], "\n";
				print "\t\t", $open_haplotypes[$existingHaploI][2], "\n";
				if($open_haplotypes[$existingHaploI][1])
				{
					my $ref_str = $open_haplotypes[$existingHaploI][1][0];
					my $haplo_str = $open_haplotypes[$existingHaploI][1][1];
					print "\t\t", $open_haplotypes[$existingHaploI][1][2], "\n";
					my $printFrom = $open_haplotypes[$existingHaploI][2];
					$printFrom = 0 if($printFrom < 0);
					print "\t\t", substr($ref_str, $printFrom, 10), "\n";
					print "\t\t", substr($haplo_str, $printFrom, 10), "\n";
				}
				else
				{
					print "\t\tREF\n";
				}
			}
			print "\n";
		}
		
		# we consume characters up to (but not including) the next reference non-gap
		my %extensions_nonRef;
		my $extensions_nonRef_length;
		foreach my $haplotype (@open_haplotypes)
		{
			my $extension;
			my $consumed_ref_start;
			my $consumed_ref = 0;
			my $consumed_ref_sequence;
			
			if($haplotype->[1] == 0)
			{

			}
			else
			{
				my $addIndex = 0;
				$consumed_ref_start = $haplotype->[2]+$addIndex+1;
				do {
					my $nextPosToConsume = $haplotype->[2]+$addIndex+1;
					die unless($nextPosToConsume < length($haplotype->[1][0]));
					my $refC = substr($haplotype->[1][0], $nextPosToConsume, 1);
					$consumed_ref_sequence .= $refC;
					if(($refC ne '-') and ($refC ne '*'))
					{
						$consumed_ref++;
					}
					my $hapC = substr($haplotype->[1][1], $nextPosToConsume, 1);
					$extension .= $hapC;
					$addIndex++;
				} while($consumed_ref < 1);
			}
			if($extension)
			{
				push(@{$extensions_nonRef{$extension}}, [$consumed_ref_start, $consumed_ref, $consumed_ref_sequence]);
				unless(defined $extensions_nonRef_length)
				{
					$extensions_nonRef_length = length($extension);
				}
				die Dumper("Length mismatch", $extension, \%extensions_nonRef, $posI, "Length mismatch") unless(length($extension) == $extensions_nonRef_length);
			}
		}
		
		# now carry out the actual extension
		my %extensions;
		foreach my $haplotype (@open_haplotypes)
		{
			my $extension;
			if($haplotype->[1] == 0)
			{
				my $refExt = $refC;
				if(defined $extensions_nonRef_length)
				{
					my $missing = $extensions_nonRef_length - length($refExt);
					die unless($missing >= 0);
					my $missingStr = '*' x $missing;
					die unless(length($missingStr) == $missing);
					$refExt .= $missingStr;
				}
				$extension = $refExt;
			}
			else
			{
				my $consumed_ref = 0;
				do {
					my $nextPosToConsume = $haplotype->[2]+1;
					die unless($nextPosToConsume < length($haplotype->[1][0]));
					my $refC = substr($haplotype->[1][0], $nextPosToConsume, 1);
					if(($refC ne '-') and ($refC ne '*'))
					{
						$consumed_ref++;
					}
					my $hapC = substr($haplotype->[1][1], $nextPosToConsume, 1);
					$extension .= $hapC;
					$haplotype->[2]++;
				} while($consumed_ref < 1);
			}
			die unless($extension);
			$haplotype->[0] .= $extension;
			$extensions{$extension}++;
		}
		
		die unless(scalar(keys %extensions));
		# print "Extensions:\n", join("\n", map {"\t'".$_."'"} keys %extensions), "\n\n";
		
		#my $this_all_equal = ( (scalar(keys %extensions) == 0) or ((scalar(keys %extensions) == 1) and (exists $extensions{$refC})) );
		my $this_all_equal = ((scalar(keys %extensions) == 1) and (exists $extensions{$refC}));
		if($posI == 0)
		{
			die unless($this_all_equal);
		}
		
		if($this_all_equal and $posI)
		{
			# close
			my $ref_span = $posI - $start_open_haplotypes;
			die unless($ref_span);
			my $reference_sequence = substr($reference_href->{$referenceSequenceID}, $start_open_haplotypes, $ref_span);
			my %alternativeSequences;
			my %uniqueRemainers;
			my @new_open_haplotypes;
			my $open_haplotypes_before = scalar(@open_haplotypes);			
			foreach my $haplotype (@open_haplotypes)
			{
				die Dumper("Length mismatch II", $ref_span+1, length($haplotype->[0]), "Length mismatch II") unless(length($haplotype->[0]) >= ($ref_span + 1));
				my $haplotype_coveredSequence = substr($haplotype->[0], 0, length($haplotype->[0])-1);
				$haplotype_coveredSequence =~ s/[\-\*]//g;
				if($haplotype_coveredSequence ne $reference_sequence)
				{
					$alternativeSequences{$haplotype_coveredSequence}++;
				}
				substr($haplotype->[0], 0, length($haplotype->[0])-1) = '';
				die unless(length($haplotype->[0]) == 1);
				
				my $k = join('--', @$haplotype);
				if(not $uniqueRemainers{$k})
				{
					push(@new_open_haplotypes, $haplotype);
					$uniqueRemainers{$k}++;
				}				
			}
			
			@open_haplotypes = @new_open_haplotypes;
			my $open_haplotypes_after = scalar(@open_haplotypes);
								
			if(keys %alternativeSequences)
			{
				my @alternativeAlleles = keys %alternativeSequences;
				
				# @alternativeAlleles = map {$_ =~ s/[\-\*]//g; $_} @alternativeAlleles;
				
				# print "Starting at position $start_open_haplotypes, have REF $reference_sequence and alternative sequences " . join(' / ', @alternativeAlleles) . "\n";
			
				if((length($reference_sequence) == 2) and (all {length($_) == 2} @alternativeAlleles))
				{
					my $firstRefChar = substr($reference_sequence, 0, 1);
					foreach my $alt (@alternativeAlleles)
					{
						die Dumper("Some problem with supposed SNP", $posI, $reference_sequence, \@alternativeAlleles, "Some problem with supposed SNP") unless(substr($alt, 0, 1) eq $firstRefChar);						
					}
					print OUT join("\t",
							$referenceSequenceID,
							$start_open_haplotypes+2,
							".",
							substr($reference_sequence, 1, 1),
							join(',', map {substr($_, 1, 1)} @alternativeAlleles),
							'.',
							'PASS',
							'.'
						), "\n";					
				}
				else
				{
					print OUT join("\t",
							$referenceSequenceID,
							$start_open_haplotypes+1,
							".",
							$reference_sequence,
							join(',', @alternativeAlleles),
							'.',
							'PASS',
							'.'
						), "\n";
				}
			}
			$start_open_haplotypes = $posI;
			
			if(($n_alignments == $opened_alignments) and ($open_haplotypes_after == 1))
			{
				last POSI;
			}
			# print "Went from $open_haplotypes_before to $open_haplotypes_after \n";
		}
		
		# $last_all_equal = $this_all_equal;
	}
	
}

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

