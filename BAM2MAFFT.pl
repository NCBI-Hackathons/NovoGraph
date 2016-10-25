#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;

unless(-e 'forMAFFT')
{
	mkdir('forMAFFT') or die "Cannot mkdir forMAFFT directory";
}

use Set::IntervalTree;

my $referenceFasta;
my $BAM;

$referenceFasta = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\reference.fa';
$BAM = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\ribosome.bam';
GetOptions (
 'referenceFasta:s' => \$referenceFasta, 
 'BAM:s' => \$BAM, 
);

{
	my $max_pos = 100;
	my @window_positions = (10, 20, 30);	
	my %window_switch_positions;	
}

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $reference_href = readFASTA($referenceFasta);
my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my $windows_info_fn = 'forMAFFT/_windowsInfo';
open(WINDOWS, '>', $windows_info_fn) or die "Cannot open $windows_info_fn";
print WINDOWS join("\t", "referenceContigID", "windowI", "firstPos_relative_to_ref", "lastPos_relative_to_ref"), "\n";

my $alignments_info_fn = 'forMAFFT/_alignments';
open(ALIGNMENTS, '>', $alignments_info_fn) or die "Cannot open $alignments_info_fn";
print ALIGNMENTS join("\t", "alignedSequenceID", "referenceContigID", "firstPositions_reference", "lastPosition_reference"), "\n";

my $alignments_gaps_info_fn = 'forMAFFT/_alignments_inWindow_onlyGaps';
open(ALIGNMENTSONLYGAPS, '>', $alignments_gaps_info_fn) or die "Cannot open $alignments_gaps_info_fn";
print ALIGNMENTSONLYGAPS join("\t", "alignedSequenceID", "referenceContigID", "windowI"), "\n";

my %saw_read_IDs;
my @sequence_ids = $sam->seq_ids();
foreach my $referenceSequenceID (@sequence_ids)
{	
	next unless(defined $reference_href->{$referenceSequenceID});

	my @alignments = $sam->get_features_by_location(-seq_id => $referenceSequenceID);
	
	print "Region $referenceSequenceID, have ", scalar(@alignments), " alignments.\n";

	my %sequences_per_window;
	my @window_positions;
	my @contig_coverage = ((0) x length($reference_href->{$referenceSequenceID}));
	my @contig_coverage_nonGap = ((0) x scalar(@contig_coverage));
	my $n_alignment = 0;
	foreach my $alignment (@alignments)
	{
		$n_alignment++;
		next unless(rand() <= 0.0001);
		my $alignment_start_pos = $alignment->start;
		my ($ref,$matches,$query) = $alignment->padded_alignment;
		# print $ref, "\n", $query, "\n\n";
		die unless(length($ref) == length($query));
		my $rel_idx = -1;	
		for(my $i = 0; $i < length($ref); $i++)
		{
			my $C_ref = substr($ref, $i, 1);
			my $C_query = substr($query, $i, 1);
			if($C_ref ne '-')
			{
				$rel_idx++;
			}
			
			if($rel_idx >= 0)
			{
				$contig_coverage[$alignment_start_pos + $rel_idx]++;
				if($C_query ne '-')
				{
					$contig_coverage_nonGap[$alignment_start_pos + $rel_idx]++;
				}				
			}
		}
	}
	
	my $scanTarget = 100;
	my $targetWindowSize = 10000;
	for(my $potentialWindowPos = 0; $potentialWindowPos <= $#contig_coverage; $potentialWindowPos++)
	{
		# print "Window: ", $potentialWindowPos, "\n";
		my $middleWindowPos = $potentialWindowPos + $targetWindowSize;
		my $minWindowsPos = $middleWindowPos - $scanTarget;
		my $maxWindowsPos = $middleWindowPos + $scanTarget;
		$minWindowsPos = 0 if($minWindowsPos < 0);
		$maxWindowsPos = $#contig_coverage if ($maxWindowsPos > $#contig_coverage);
		
		my @actualWindowPos_options;
		my %actualWindowPos_options_missing;
		for(my $actualWindowPos = $minWindowsPos; $actualWindowPos <= $maxWindowsPos; $actualWindowPos++)
		{
			my $rate_missing = 0;
			if($contig_coverage[$actualWindowPos] > 0)
			{
				die unless($contig_coverage_nonGap[$actualWindowPos] <= $contig_coverage[$actualWindowPos]);
				my $rate_missing = $contig_coverage_nonGap[$actualWindowPos] / $contig_coverage[$actualWindowPos];
				push(@actualWindowPos_options, $actualWindowPos);
				$actualWindowPos_options_missing{$rate_missing} = $rate_missing;
			}
		}
		
		if(scalar(@actualWindowPos_options) == 0)
		{
			last;
		}
		@actualWindowPos_options = sort {$actualWindowPos_options_missing{$a} <=> $actualWindowPos_options_missing{$b}} @actualWindowPos_options;
		if(scalar(@actualWindowPos_options))
		{
			die unless($actualWindowPos_options_missing{$actualWindowPos_options[0]} <= $actualWindowPos_options_missing{$actualWindowPos_options[0]});
		}
		my $selectedWindowPos = $actualWindowPos_options[0];
		push(@window_positions, $selectedWindowPos);
		$potentialWindowPos = $selectedWindowPos;
	}
	
	my $intervalTree_windows = Set::IntervalTree->new;	
	my %window_switch_positions;
	my $max_pos = $#contig_coverage;
	$window_switch_positions{$window_positions[0]} = 1;
	$intervalTree_windows->insert('w'.0,0,$window_positions[0]);
	$window_switch_positions{$window_positions[0]} = 'w1';
	for(my $windowI = 1; $windowI <= $#window_positions; $windowI++)
	{
		$intervalTree_windows->insert('w'.$windowI,$window_positions[$windowI-1],$window_positions[$windowI]);
		$window_switch_positions{$window_positions[$windowI]} = 'w'.($windowI+1);
		# print "Insert ", $window_positions[$windowI-1], " ", $window_positions[$windowI], " as $windowI\n";
	}
	$intervalTree_windows->insert('w'.($#window_positions+1), $window_positions[$#window_positions], $max_pos+1);
	# print "Insert ", $window_positions[$#window_positions], " ", $max_pos+1, " as ", $#window_positions+1, "\n";
	
	print WINDOWS join("\t", $referenceSequenceID, 0, 0, $window_positions[0] - 1), "\n";
	for(my $windowI = 0; $windowI <= $#window_positions; $windowI++)
	{
		if($windowI > 0)
		{
			die unless($window_positions[$windowI] > $window_positions[$windowI-1]);
			die unless(($window_positions[$windowI] - $window_positions[$windowI-1]) > 1);
		}
		
		{			
			my $windowPos = $window_positions[$windowI];
			my $retrieval_aref = $intervalTree_windows->fetch($windowPos,$windowPos+1);	
			die unless(scalar(@$retrieval_aref) == 1);
			
			die Dumper("Problem with windows I",  $windowI, [$windowPos,$windowPos+1], $retrieval_aref, \@window_positions) unless($retrieval_aref->[0] eq ('w'.($windowI+1)));
			die unless($window_switch_positions{$windowPos} eq ('w'.($windowI+1)));
		}
		
		{
			my $previousWindowLastPos = $window_positions[$windowI] - 1;
			my $retrieval_aref = $intervalTree_windows->fetch($previousWindowLastPos,$previousWindowLastPos+1);	
			die unless(scalar(@$retrieval_aref) == 1);
			die unless($retrieval_aref->[0] eq ('w'.$windowI));
		}
		
		print WINDOWS join("\t", $referenceSequenceID, $windowI + 1, $window_positions[$windowI], (($windowI != $#window_positions) ? $window_positions[$windowI+1] - 1 : $max_pos )), "\n";
	}
	
	print Dumper(\@window_positions, $#contig_coverage), "\n";
	
	die unless($#window_positions >= 0);
	
	foreach my $alignment (@alignments)
	{
		next unless(rand() <= 0.001);
		
		my $readID = $alignment->query->name;
		if(exists $saw_read_IDs{$readID})
		{
			warn "Observed more than one alignment for read ID $readID - this should not happen!";
			next;
		}
		$saw_read_IDs{$readID}++;
		
		my $alignment_start_pos = $alignment->start;
		my ($ref,$matches,$query) = $alignment->padded_alignment;
		die unless(length($ref) == length($query));
		my $rel_idx = -1;	
		my $firstPosition_reference;
		my $lastPosition_reference;
		my $currentWindow;
		my $unaccountedSequence;
		for(my $i = 0; $i < length($ref); $i++)
		{
			my $C_ref = substr($ref, $i, 1);
			my $C_query = substr($query, $i, 1);
			if($C_ref ne '-')
			{
				$rel_idx++;
			}
			
			if($rel_idx >= 0)
			{
				my $position_along_reference = $alignment_start_pos + $rel_idx;
				$firstPosition_reference = $position_along_reference unless(defined $firstPosition_reference);
				$lastPosition_reference = $position_along_reference;
				
				if(not defined $currentWindow)
				{
					my $result_windows = $intervalTree_windows->fetch($position_along_reference,$position_along_reference+1);
					die unless(scalar(@$result_windows) == 1);
					$currentWindow = $result_windows->[0];
				}
				else
				{
					if(exists $window_switch_positions{$position_along_reference})
					{
						$currentWindow = $window_switch_positions{$position_along_reference};
					}
				}
					
				if($unaccountedSequence)
				{
					$sequences_per_window{$currentWindow}{$readID} .= $unaccountedSequence;
					$unaccountedSequence = '';
				}
				
				$sequences_per_window{$currentWindow}{$readID} .= $C_query;				
			}
			else
			{
				$unaccountedSequence .= $C_query;
			}
		}
		
		die unless((defined $firstPosition_reference) and (defined $lastPosition_reference));
		print ALIGNMENTS join("\t", $readID, $referenceSequenceID, $firstPosition_reference, $lastPosition_reference), "\n";
	}

	foreach my $windowID (0 .. scalar(@window_positions))
	{
		my $window_idx_key = 'w' . $windowID;
		my $output_fn = 'forMAFFT/' . $referenceSequenceID . '_' . $windowID . '.fa';
		open(MAFFTOUT, '>', $output_fn) or die "Cannot open $output_fn";
		foreach my $sequenceID (keys %{$sequences_per_window{$window_idx_key}})
		{
			my $sequence_for_emission = $sequences_per_window{$window_idx_key}{$sequenceID};
			$sequence_for_emission =~ s/[\-_]//g;
			if(length($sequence_for_emission))
			{
				print MAFFTOUT '>', $sequenceID, "\n";
				print MAFFTOUT $sequences_per_window{$window_idx_key}{$sequenceID}, "\n";
			}
			else
			{
				print ALIGNMENTSONLYGAPS join("\t", $sequenceID, $referenceSequenceID, $windowID), "\n";			
			}
		}
		close(MAFFTOUT);		
	}
}

close(WINDOWS);
close(ALIGNMENTSONLYGAPS);
												 
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
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}

