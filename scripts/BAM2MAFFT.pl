#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use Bio::DB::HTS;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;
use File::Path;

$| = 1;

## perl BAM2MAFFT.pl --BAM <path to BAM, output of FIND_GLOBAL_ALIGNMENTS.pl>
##                  --referenceFasta <path to reference FASTA>
##                  --readsFasta <path to contigs FASTA>
##                  --outputDirectory <path to output directory for MAFFT, e.g. '/forMAFFT'>
##                  --inputTruncatedReads <path to output of FIND_GLOBAL_ALIGNMENTS.pl, 'outputTruncatedReads'>
##
## Example command
## ./BAM2MAFFT.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /data/projects/phillippy/projects/hackathon/shared/contigs/AllContigs.fa --outputDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT --inputTruncatedReads /data/projects/phillippy/projects/hackathon/intermediate_files/truncatedReads
 

my $referenceFasta;
my $BAM;
my $outputDirectory;
my $inputTruncatedReads;
my $readsFasta;

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputDirectory:s' => \$outputDirectory,	
	'readsFasta:s' => \$readsFasta,	
	'inputTruncatedReads:s' => \$inputTruncatedReads,	
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --inputTruncatedReads" unless($inputTruncatedReads);
die "Please specify --outputDirectory" unless($outputDirectory);
die "Please specify --readsFasta" unless($readsFasta);
die "--readsFasta $readsFasta not existing" unless(-e $readsFasta);

die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

die "--inputTruncatedReads $inputTruncatedReads not existing" unless(-e $inputTruncatedReads);

die "Security check - $outputDirectory will be deleted" unless($outputDirectory =~ /mafft/i);
rmtree($outputDirectory);
unless((-e $outputDirectory) and (-d $outputDirectory))
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory directory";
}

my $reads_href = readFASTA($readsFasta, 0);

my $reference_href = readFASTA($referenceFasta);
my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $BAM);

my %truncatedReads;
open(TRUNCATED, '<', $inputTruncatedReads) or die "Cannot open $inputTruncatedReads";
while(<TRUNCATED>)
{
	my $line = $_;
	chomp($line);
	$truncatedReads{$line}++;
}

close(TRUNCATED);
my $windows_info_fn =  $outputDirectory . '/_windowsInfo';
open(WINDOWS, '>', $windows_info_fn) or die "Cannot open $windows_info_fn";
print WINDOWS join("\t", "referenceContigID", "chrDir", "windowI", "firstPos_relative_to_ref", "lastPos_relative_to_ref", "lastPos_BAMcoverage", "lastPos_BAMcoverage_nonGap"), "\n";

my $alignments_info_fn = $outputDirectory . '/_alignments';
open(ALIGNMENTS, '>', $alignments_info_fn) or die "Cannot open $alignments_info_fn";
print ALIGNMENTS join("\t", "alignedSequenceID", "referenceContigID", "chrDir", "firstPositions_reference", "lastPosition_reference"), "\n";

my $alignments_gaps_info_fn = $outputDirectory . '/_alignments_inWindow_onlyGaps';
open(ALIGNMENTSONLYGAPS, '>', $alignments_gaps_info_fn) or die "Cannot open $alignments_gaps_info_fn";
print ALIGNMENTSONLYGAPS join("\t", "alignedSequenceID", "referenceContigID", "chrDir", "windowI"), "\n";

my %saw_ref_IDs;
my %saw_read_IDs;
my @sequence_ids = $sam->seq_ids();
foreach my $referenceSequenceID (@sequence_ids)
{	
	next unless($referenceSequenceID =~ /chr[XY\d]+/);
	unless(exists $reference_href->{$referenceSequenceID})
	{
		warn "No reference sequence for $referenceSequenceID";
		next;
	}
	next unless(length($reference_href->{$referenceSequenceID}) > 20000);
	
	my $chrDir = $referenceSequenceID;
	$chrDir =~ s/\W//g;
	die "Duplicate directory $chrDir?" if(exists $saw_ref_IDs{$chrDir});
	mkdir( $outputDirectory . '/' . $chrDir) or die "Cannot mkdir $chrDir";
	
	print "Processing $referenceSequenceID", ", length ", length($reference_href->{$referenceSequenceID}), "\n";
	die "Length discrepancy between supplied FASTA reference and BAM index: " . $sam->length($referenceSequenceID) . " vs " . length($reference_href->{$referenceSequenceID}) unless($sam->length($referenceSequenceID) == length($reference_href->{$referenceSequenceID}));
	next unless(defined $reference_href->{$referenceSequenceID});
	
	my %sequences_per_window;
	my @window_positions;
	my @contig_coverage = ((0) x length($reference_href->{$referenceSequenceID}));
	my @contig_coverage_nonGap = ((0) x scalar(@contig_coverage));
	
	my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);
	my $n_alignment = 0;
	while(my $alignment = $alignment_iterator->next_seq)
	{
		$n_alignment++;		
		print "\r\t\tProcessing ", $n_alignment, "...   ";
		
		my $cigar  = $alignment->cigar_str;
		while($cigar =~ /^(\d+)([MIDHS])/)
		{
			my $number = $1;
			my $action = $2;
			$cigar =~ s/^(\d+)([MIDHS])//;		
			if(($action eq 'H') or ($action eq 'S'))
			{
				die "BAM $BAM contains H or S in CIGAR string - illegal, we want global, non-clipped alignments.";
			}
		}

	
		my $alignment_start_pos = $alignment->start;
		my ($ref,$matches,$query) = $alignment->padded_alignment;
        
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
	print "\n";
	
	my $scanTarget = 100;
	my $targetWindowSize = 10000;
	for(my $potentialWindowPos = 0; $potentialWindowPos <= $#contig_coverage; $potentialWindowPos++)
	{
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
	die unless($#window_positions >= 0);
	
	my $intervalTree_windows = Set::IntervalTree->new;	
	my %window_switch_positions;
	my $max_pos = $#contig_coverage;
	$window_switch_positions{$window_positions[0]} = 1;
	$intervalTree_windows->insert('w0', 0, $window_positions[0]);
	$window_switch_positions{$window_positions[0]} = 'w1';
	for(my $windowI = 1; $windowI <= $#window_positions; $windowI++)
	{
		my $windowID_for_tree = 'w'.$windowI;
		$intervalTree_windows->insert($windowID_for_tree, $window_positions[$windowI-1], $window_positions[$windowI]);
		$window_switch_positions{$window_positions[$windowI]} = 'w'.($windowI+1);		
	}
	
	my $last_windowID_for_tree = 'w'.($#window_positions+1);
	$intervalTree_windows->insert($last_windowID_for_tree, $window_positions[$#window_positions], $max_pos+1);
	
	my $firstWindow_lastPos = $window_positions[0] - 1;
	print WINDOWS join("\t", $referenceSequenceID, $chrDir, 0, 0, $firstWindow_lastPos, $contig_coverage[$firstWindow_lastPos], $contig_coverage_nonGap[$firstWindow_lastPos]), "\n";
	
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
			$chrDir, 
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
	
	my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);
	my $n_alignment = 0;
	while(my $alignment = $alignment_iterator->next_seq)
	{
		$n_alignment++;		
		print "\r\t\tProcessing ", $n_alignment, "...   ";

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
		print ALIGNMENTS join("\t", $readID, $referenceSequenceID, $chrDir, $firstPosition_reference, $lastPosition_reference), "\n";
	}	
	print "\n";

	my %saw_sequence_already;
	my %runningSequencesForReconstruction;
	foreach my $windowID (0 .. scalar(@window_positions))
	{
		my %sequenceID_not_seen = map {$_ => 1} keys %runningSequencesForReconstruction;
		
		my $window_idx_key = 'w' . $windowID;
		my $output_fn = $outputDirectory . '/' . $chrDir . '/' . $referenceSequenceID . '_' . $windowID . '.fa';
		my $is_open = 0;
		foreach my $sequenceID (keys %{$sequences_per_window{$window_idx_key}})
		{
			my $sequence_for_emission = $sequences_per_window{$window_idx_key}{$sequenceID};
			#$sequence_for_emission =~ s/[\-_]//g;
			$sequenceID_not_seen{$sequenceID} = 0;
			
			my $sequence_for_emission_noGaps = $sequence_for_emission;
			$sequence_for_emission_noGaps =~ s/[\-_]//g;
			$runningSequencesForReconstruction{$sequenceID} .= $sequence_for_emission_noGaps;
			
			my $sequenceID_for_print = $sequenceID;
			if(not $saw_sequence_already{$sequenceID})
			{
				# this is the first time we see that sequence
				die unless(length($sequence_for_emission_noGaps));
				$saw_sequence_already{$sequenceID} = 1;
				$sequenceID_for_print .= "_FIRST";
			}
			if(length($sequence_for_emission))
			{
				unless($is_open)
				{
					open(MAFFTOUT, '>', $output_fn) or die "Cannot open $output_fn";
					$is_open = 1;
				}
				print MAFFTOUT '>', $sequenceID_for_print, "\n";
				print MAFFTOUT $sequence_for_emission, "\n";
			}
			else
			{
				print ALIGNMENTSONLYGAPS join("\t", $sequenceID, $referenceSequenceID, $chrDir, $windowID), "\n";			
			}
		}
		
		if($is_open)
		{
			close(MAFFTOUT);		
		}
		
		foreach my $sequenceID (grep {$sequenceID_not_seen{$_}} keys %sequenceID_not_seen)
		{
			unless(exists $reads_href->{$sequenceID})
			{
				warn "No truth sequence for $sequenceID";
				delete($runningSequencesForReconstruction{$sequenceID});				
				next;
			}
			
			my $trueSequence = $reads_href->{$sequenceID};
			my $trueSequence_revCmp = reverseComplement($trueSequence);
			
			my $supposedSequence = $runningSequencesForReconstruction{$sequenceID};
			if(($supposedSequence eq $trueSequence) or ($supposedSequence eq $trueSequence_revCmp))
			{
				# print "OK";
			}
			else
			{
				if(not exists $truncatedReads{$sequenceID})
				{
					print "Disagreement for $sequenceID\n";
					print "\t", "length(\$trueSequence)", ": ", length($trueSequence), "\n";
					print "\t", "length(\$supposedSequence)", ": ", length($supposedSequence), ", ", substr($supposedSequence, 0, 10), "\n";
				}
			}
			delete($runningSequencesForReconstruction{$sequenceID});
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

