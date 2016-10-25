#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;
$| = 1;

my $referenceFasta;
my $BAM;
my $outputDirectory = 'forMafft2';
GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputDirectory:s' => \$outputDirectory,	
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

unless((-e $outputDirectory) and (-d $outputDirectory))
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory directory";
}

my $reference_href = readFASTA($referenceFasta);
my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my $windows_info_fn =  $outputDirectory . '/_windowsInfo';
open(WINDOWS, '>', $windows_info_fn) or die "Cannot open $windows_info_fn";
print WINDOWS join("\t", "referenceContigID", "windowI", "firstPos_relative_to_ref", "lastPos_relative_to_ref", "lastPos_BAMcoverage", "lastPos_BAMcoverage_nonGap"), "\n";

my $alignments_info_fn = $outputDirectory . '/_alignments';
open(ALIGNMENTS, '>', $alignments_info_fn) or die "Cannot open $alignments_info_fn";
print ALIGNMENTS join("\t", "alignedSequenceID", "referenceContigID", "firstPositions_reference", "lastPosition_reference"), "\n";

my $alignments_gaps_info_fn = $outputDirectory . '/_alignments_inWindow_onlyGaps';
open(ALIGNMENTSONLYGAPS, '>', $alignments_gaps_info_fn) or die "Cannot open $alignments_gaps_info_fn";
print ALIGNMENTSONLYGAPS join("\t", "alignedSequenceID", "referenceContigID", "windowI"), "\n";

my %saw_read_IDs;
my @sequence_ids = $sam->seq_ids();
foreach my $referenceSequenceID (@sequence_ids)
{	
	next unless($referenceSequenceID =~ /chr/);
	unless(exists $reference_href->{$referenceSequenceID})
	{
		warn "No reference sequence for $referenceSequenceID";
		next;
	}

	print "Processing $referenceSequenceID", ", length ", length($reference_href->{$referenceSequenceID}), "\n";
	die "Length discrepancy between supplied FASTA reference and BAM index: " . $sam->length($referenceSequenceID) . " vs " . length($reference_href->{$referenceSequenceID}) unless($sam->length($referenceSequenceID) == length($reference_href->{$referenceSequenceID}));
	next unless(defined $reference_href->{$referenceSequenceID});

	# print "\tRegion $referenceSequenceID, have ", scalar(@alignments), " alignments.\n";
	
	my %sequences_per_window;
	my @window_positions;
	my @contig_coverage = ((0) x length($reference_href->{$referenceSequenceID}));
	my @contig_coverage_nonGap = ((0) x scalar(@contig_coverage));
		
	my @alignmentExtractionWindows = get_windows_for_refLength(length($reference_href->{$referenceSequenceID}));
	for(my $alignmentExtractionWindowI = 0; $alignmentExtractionWindowI <= $#alignmentExtractionWindows; $alignmentExtractionWindowI++)
	{
		my @alignmentExtractionWindow = @{$alignmentExtractionWindows[$alignmentExtractionWindowI]};
		
		my @alignments = $sam->get_features_by_location(-seq_id => $referenceSequenceID, -start => $alignmentExtractionWindow[0], -end => $alignmentExtractionWindow[1]);
	
		my $n_alignment = 0;
		foreach my $alignment (@alignments)
		{
			$n_alignment++;		
			print "\r\t\tProcessing ", $n_alignment, "/", scalar(@alignments), " alignments in window ", ($alignmentExtractionWindowI+1), "/", scalar(@alignmentExtractionWindows), "...   ";
			next unless(($alignment->start >= $alignmentExtractionWindow[0]) and ($alignment->start <= $alignmentExtractionWindow[1]));
			
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
	}
	print "\n";
	
	my $scanTarget = 100;
	my $targetWindowSize = 10000;
	for(my $potentialWindowPos = 0; $potentialWindowPos <= $#contig_coverage; $potentialWindowPos++)
	{
		# print "Window: ", $potentialWindowPos, "\n";
		my $middleWindowPos = $potentialWindowPos + $targetWindowSize;
		my $minWindowsPos = $middleWindowPos - $scanTarget;
		my $maxWindowsPos = $middleWindowPos + $scanTarget;
		last if($minWindowsPos >= $#contig_coverage);
		
		$minWindowsPos = 0 if($minWindowsPos < 0);
		$maxWindowsPos = $#contig_coverage if ($maxWindowsPos > $#contig_coverage);
		
		my @actualWindowPos_options;
		my %actualWindowPos_options_missing;
		die unless($minWindowsPos <= $maxWindowsPos);
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
			push(@actualWindowPos_options, $middleWindowPos);
			$actualWindowPos_options_missing{$middleWindowPos} = 0;
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
	
	print "\tRegion $referenceSequenceID, have ", scalar(@window_positions), " windows.\n";
	# print "Window starting positions (first window starting at 0 is implicit):\n";
	# print join("\n", map {' - ' . $_} @window_positions), "\n";
	die unless($#window_positions >= 0);
	
	my $intervalTree_windows = Set::IntervalTree->new;	
	my %window_switch_positions;
	my $max_pos = $#contig_coverage;
	$window_switch_positions{$window_positions[0]} = 1;
	$intervalTree_windows->insert('w0', 0, $window_positions[0]);
	# print "Insert ", 0, " ", $window_positions[0], " as w0\n";	
	$window_switch_positions{$window_positions[0]} = 'w1';
	for(my $windowI = 1; $windowI <= $#window_positions; $windowI++)
	{
		my $windowID_for_tree = 'w'.$windowI;
		$intervalTree_windows->insert($windowID_for_tree, $window_positions[$windowI-1], $window_positions[$windowI]);
		# print "Insert ", $window_positions[$windowI-1], " ", $window_positions[$windowI], " as $windowID_for_tree\n";
		$window_switch_positions{$window_positions[$windowI]} = 'w'.($windowI+1);		
	}
	
	my $last_windowID_for_tree = 'w'.($#window_positions+1);
	$intervalTree_windows->insert($last_windowID_for_tree, $window_positions[$#window_positions], $max_pos+1);
	# print "Insert ", $window_positions[$#window_positions], " ", $max_pos+1, " as $last_windowID_for_tree\n";
	
	my $firstWindow_lastPos = $window_positions[0] - 1;
	print WINDOWS join("\t", $referenceSequenceID, 0, 0, $firstWindow_lastPos, $contig_coverage[$firstWindow_lastPos], $contig_coverage_nonGap[$firstWindow_lastPos]), "\n";
	
	my $first_window_idx_key = 'w0';
	my $first_window_referenceSequence = substr($reference_href->{$referenceSequenceID}, 0, $window_positions[0]);
	$sequences_per_window{$first_window_idx_key}{ref} = $first_window_referenceSequence;
	
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
			
			die Dumper("Problem with windows I",  $windowI, [$windowPos,$windowPos+1], $retrieval_aref) unless($retrieval_aref->[0] eq ('w'.($windowI+1)));
			die unless($window_switch_positions{$windowPos} eq ('w'.($windowI+1)));
		}
		
		{
			my $previousWindowLastPos = $window_positions[$windowI] - 1;
			my $retrieval_aref = $intervalTree_windows->fetch($previousWindowLastPos,$previousWindowLastPos+1);	
			die unless(scalar(@$retrieval_aref) == 1);
			die unless($retrieval_aref->[0] eq ('w'.$windowI));
		}
		
		my $window_lastPos = (($windowI != $#window_positions) ? $window_positions[$windowI+1] - 1 : $max_pos );
		
		print WINDOWS join("\t",
			$referenceSequenceID,
			$windowI + 1,
			$window_positions[$windowI],
			$window_lastPos,
			$contig_coverage[$window_lastPos],
			$contig_coverage_nonGap[$window_lastPos],
		), "\n";
		
		my $window_idx_key = 'w' . ($windowI + 1);
		my $referenceSequence = substr($reference_href->{$referenceSequenceID}, $window_positions[$windowI], $window_lastPos - $window_positions[$windowI] + 1);
		$sequences_per_window{$window_idx_key}{ref} = $referenceSequence;
	}
	
	
	for(my $alignmentExtractionWindowI = 0; $alignmentExtractionWindowI <= $#alignmentExtractionWindows; $alignmentExtractionWindowI++)
	{
		my @alignmentExtractionWindow = @{$alignmentExtractionWindows[$alignmentExtractionWindowI]};
		
		my @alignments = $sam->get_features_by_location(-seq_id => $referenceSequenceID, -start => $alignmentExtractionWindow[0], -end => $alignmentExtractionWindow[1]);
	
		my $n_alignment = 0;
		foreach my $alignment (@alignments)
		{
			$n_alignment++;
			print "\r\t\tProcessing ", $n_alignment, "/", scalar(@alignments), " alignments in window ", ($alignmentExtractionWindowI+1), "/", scalar(@alignmentExtractionWindows), "...   ";
			next unless(($alignment->start >= $alignmentExtractionWindow[0]) and ($alignment->start <= $alignmentExtractionWindow[1]));

			my $readID = $alignment->query->name;
			if(exists $saw_read_IDs{$readID})
			{
				die "Observed more than one alignment for read ID $readID - this should not happen!";
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
	}
	print "\n";

	foreach my $windowID (0 .. scalar(@window_positions))
	{
		my $window_idx_key = 'w' . $windowID;
		my $output_fn = $outputDirectory . '/' . $referenceSequenceID . '_' . $windowID . '.fa';
		my $is_open = 0;
		foreach my $sequenceID (keys %{$sequences_per_window{$window_idx_key}})
		{
			my $sequence_for_emission = $sequences_per_window{$window_idx_key}{$sequenceID};
			$sequence_for_emission =~ s/[\-_]//g;
			if(length($sequence_for_emission))
			{
				unless($is_open)
				{
					open(MAFFTOUT, '>', $output_fn) or die "Cannot open $output_fn";
					$is_open = 1;
				}
				print MAFFTOUT '>', $sequenceID, "\n";
				print MAFFTOUT $sequence_for_emission, "\n";
			}
			else
			{
				print ALIGNMENTSONLYGAPS join("\t", $sequenceID, $referenceSequenceID, $windowID), "\n";			
			}
		}
		
		if($is_open)
		{
			close(MAFFTOUT);		
		}
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

sub get_windows_for_refLength
{
	my $refLength = shift;
	die unless(defined $refLength);
	my $desiredWindowLength = 100000;
	
	if($refLength <= $desiredWindowLength)
	{
		return ([1, $refLength]);
	}
	
	my @forReturn;
	my $currentEnd = $desiredWindowLength;
	while($currentEnd < $refLength)
	{
		push(@forReturn, [$currentEnd - $desiredWindowLength + 1, $currentEnd]);
		$currentEnd += $desiredWindowLength;
	}
	if($currentEnd > $refLength)
	{
		$currentEnd = $refLength;
	}
	die if($forReturn[$#forReturn][1] == $currentEnd);
	
	my $lastWindowStart = $forReturn[$#forReturn][1] + 1;
	die unless($lastWindowStart <= $currentEnd);
	push(@forReturn, [$lastWindowStart, $currentEnd]);
	
	return @forReturn;
}