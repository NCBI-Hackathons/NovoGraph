#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
# use Bio::DB::HTS;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;
use File::Path;

$| = 1;

## BAM2MAFFT.pl --SAM <path to SAM, output of FIND_GLOBAL_ALIGNMENTS.pl>
##                  --referenceFasta <path to reference FASTA>
##                  --readsFasta <path to contigs FASTA>
##                  --outputDirectory <path to output directory for MAFFT, e.g. '/forMAFFT'>
##                  --inputTruncatedReads <path to output of FIND_GLOBAL_ALIGNMENTS.pl, 'outputTruncatedReads'>
##
## Example command
## ./BAM2MAFFT.pl --SAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT.sam 
##                --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa 
##                --readsFasta /data/projects/phillippy/projects/hackathon/shared/contigs/AllContigs.fa 
##                --outputDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT 
##                --inputTruncatedReads /data/projects/phillippy/projects/hackathon/intermediate_files/truncatedReads
 
my $referenceFasta;
my $SAM;
my $outputDirectory;
my $inputTruncatedReads;
my $readsFasta;
my $processNonChrReferenceContigs = 0;
my $bin_sam2alignment;
my $samtools_path;
GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'SAM:s' => \$SAM, 
	'outputDirectory:s' => \$outputDirectory,	
	'readsFasta:s' => \$readsFasta,	
	'inputTruncatedReads:s' => \$inputTruncatedReads,	
	'processNonChrReferenceContigs:s' => \$processNonChrReferenceContigs,	
	'sam2alignment_executable:s' => \$bin_sam2alignment,	
	'samtools_path:s' => \$samtools_path,		
);

die "Please specify --SAM" unless($SAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --inputTruncatedReads" unless($inputTruncatedReads);
die "Please specify --outputDirectory" unless($outputDirectory);
die "Please specify --readsFasta" unless($readsFasta);
die "--readsFasta $readsFasta not existing" unless(-e $readsFasta);

die "--SAM $SAM not existing" unless(-e $SAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

die "--inputTruncatedReads $inputTruncatedReads not existing" unless(-e $inputTruncatedReads);
die "--sam2alignment_executable $bin_sam2alignment not present; Please run 'make' in the directory /src and specify the path of the executable via --sam2alignment_executable" unless(-e $bin_sam2alignment);

unless($samtools_path)
{
	$samtools_path = `which samtools`;
	$samtools_path =~ s/[\n\r]//g;
	unless($samtools_path and -x $samtools_path)
	{
		die "Can't determine path to samtools - please specify --samtools_path";
	}
}	
die unless(-x $samtools_path);


die "Security check - $outputDirectory will be deleted" unless($outputDirectory =~ /mafft/i);
rmtree($outputDirectory);
unless((-e $outputDirectory) and (-d $outputDirectory))
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory directory";
}

indexFastaIfNecessary($readsFasta);

my $fn_SAM_alignments = $SAM . '.alignments';
my $cmd_doConvert = qq($bin_sam2alignment $SAM $referenceFasta > $fn_SAM_alignments);
system($cmd_doConvert) and die "Cannot execute command: $cmd_doConvert";

my $reference_href = readFASTA($referenceFasta);

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



my %processedReferenceIDs;
my %saw_ref_IDs;
my @runningAlignments;
my $runningReferenceID;
my $process_collected_read_data = sub {
	if(scalar(@runningAlignments) == 0)
	{
		return;
	}
	
	die unless(defined $runningReferenceID);
	if($processedReferenceIDs{$runningReferenceID})
	{
		die "Reference ID $runningReferenceID has already been processed - is the input sorted?";
	}
	$processedReferenceIDs{$runningReferenceID}++;
	
	my $referenceSequenceID = $runningReferenceID;
	
	unless($processNonChrReferenceContigs)
	{
		next unless($referenceSequenceID =~ /chr[XY\d]+/);
	}
	
	unless(exists $reference_href->{$referenceSequenceID})
	{
		die "No reference sequence for $referenceSequenceID";
		next;
	}
	next unless(length($reference_href->{$referenceSequenceID}) > 20000);
	
	my $chrDir = $referenceSequenceID;
	$chrDir =~ s/\W//g;
	die "Duplicate directory $chrDir?" if(exists $saw_ref_IDs{$chrDir});
	mkdir( $outputDirectory . '/' . $chrDir) or die "Cannot mkdir $chrDir";
	
	print "Processing $referenceSequenceID", ", length ", length($reference_href->{$referenceSequenceID}), "\n";
	
	my %sequences_per_window;
	my @window_positions;
	my @contig_coverage = ((0) x length($reference_href->{$referenceSequenceID}));
	my @contig_coverage_nonGap = ((0) x scalar(@contig_coverage));
	
	my $n_alignment = 0;
	foreach my $alignment (@runningAlignments)
	{
		$n_alignment++;		
		print "\r\t\tProcessing ", $n_alignment, "...   ";
		my $readID = $alignment->{readID};
	
		my $read_sequence = getSequenceFromIndexedFasta($readsFasta, $readID);
		unless(($alignment->{firstPos_read} == 0) and ($alignment->{lastPos_read} == ($alignment->{readLength}-1)))
		{
			die Dumper("Illegal - we want global, non-clipped alignments");
		}


		my $alignment_start_pos = $alignment->{firstPos_reference};
		my $ref = $alignment->{alignment_reference};
		my $query = $alignment->{alignment_read};
        
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
	
	
	
	$n_alignment = 0;
	foreach my $alignment (@runningAlignments)
	{
		$n_alignment++;		
		print "\r\t\tProcessing ", $n_alignment, "...   ";
		my $readID = $alignment->{readID};
			
		my $alignment_start_pos = $alignment->{firstPos_reference};
		my $ref = $alignment->{alignment_reference};
		my $query = $alignment->{alignment_read};
        
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
			my $read_sequence = getSequenceFromIndexedFasta($readsFasta, $sequenceID);
		
			# unless(exists $reads_href->{$sequenceID})
			# {
				# die "No truth sequence for $sequenceID";
				# delete($runningSequencesForReconstruction{$sequenceID});				
				# next;
			# }
			
			my $trueSequence = uc($read_sequence);
			my $trueSequence_revCmp = uc(reverseComplement($trueSequence));
			
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
	
	@runningAlignments = ();
};

my %processedReadIDs;
open(ALIGNMENTS, '<', $fn_SAM_alignments) or die "Cannot open $fn_SAM_alignments";
while(<ALIGNMENTS>)
{
	my $header = $_; chomp($header);
	my $aligned_reference = <ALIGNMENTS>; chomp($aligned_reference);
	my $aligned_read = <ALIGNMENTS>; chomp($aligned_read);
	my $aligned_qualities = <ALIGNMENTS>; chomp($aligned_qualities);
	die unless($aligned_reference and $aligned_read and $aligned_qualities);
	die unless(length($aligned_reference) == length($aligned_read));
	print length($aligned_reference), "\n";
	
	my @header_line_fields = split(/ /, $header);
	my $readID = $header_line_fields[0];
	my $refPositions = $header_line_fields[1];
	die unless($refPositions =~ /^(\w+):(\d+)-(\d+)$/);
	my $referenceSequenceID = $1;
	my $alignmentStart_1based = $2;
	my $alignmentStop_1based = $3;
	die unless($alignmentStart_1based <= $alignmentStop_1based); # otherwise RC
	
	my $readPositions = $header_line_fields[2];
	die unless($readPositions =~ /^read:(\d+)-(\d+)$/);
	my $alignmentStartInRead_1based = $1;
	my $alignmentStopInRead_1based = $2;
	
	my $readLength = $header_line_fields[8];
	die unless($readLength =~ /readLength=(\d+)/);
	$readLength = $1;
	
	if($processedReadIDs{$readID})
	{
		die "Read ID $readID appears more than once - abort!";
	}
	$processedReadIDs{$readID}++;
	
	if((defined $runningReferenceID) and ($runningReferenceID ne $referenceSequenceID))
	{
		$process_collected_read_data->();
	}
	$runningReferenceID = $referenceSequenceID;
	
	my $read_href = convertAlignmentToHash($readID, [$referenceSequenceID, $alignmentStart_1based, $alignmentStop_1based], [$alignmentStartInRead_1based, $alignmentStopInRead_1based], [$aligned_reference, $aligned_read], $readLength);	
	push(@runningAlignments, $read_href);
}
close(ALIGNMENTS);
$process_collected_read_data->();


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
	$kMer =~ tr/ACGTacgt/TGCAtgca/;
	return scalar(reverse($kMer));
}

my $warnings = 0;
sub convertAlignmentToHash
{
	my $readID = shift;
	my $referenceCoordinates_aref = shift;
	my $readCoordinates_aref = shift;
	my $alignment = shift;
	my $readLength = shift;
	
	die unless(length($alignment->[0]) eq length($alignment->[1]));
	die unless(defined $readLength);
	
	# $readID, [$referenceSequenceID, $alignmentStartInRead_1based, $alignmentStop_1based], [$alignmentStartInRead_1based, $alignmentStopInRead_1based], [$aligned_reference, $aligned_read]
	
	# positions
	
	my $chromosome = $referenceCoordinates_aref->[0];
	my $firstPos_reference_0based = $referenceCoordinates_aref->[1] - 1;
	die unless($firstPos_reference_0based >= 0);
	
	my $lastPos_reference_0based = $referenceCoordinates_aref->[2] - 1;
	
	die unless($firstPos_reference_0based <= $lastPos_reference_0based);

	my $firstPos_read_0based = $readCoordinates_aref->[0] - 1;
	my $lastPos_read_0based = $readCoordinates_aref->[1] - 1;
	
	my $isReverseComplement = 0;
	if($firstPos_read_0based > $lastPos_read_0based)
	{
		$isReverseComplement = 1;
		my $third = $firstPos_read_0based;
		$firstPos_read_0based = $lastPos_read_0based;
		$lastPos_read_0based = $third;
	}
	die Dumper("Read positions", $referenceCoordinates_aref, $readCoordinates_aref) unless($firstPos_read_0based <= $lastPos_read_0based);

	my $strand = ($isReverseComplement) ? -1 : 1;
	die unless(($strand == 1) or ($strand == -1));
	$strand = ($strand == 1) ? '+' : '-';
	
	# alignments
	
	my $alignment_reference = uc($alignment->[0]);
	my $alignment_read = uc($alignment->[1]);
	
	# check that all alignment sequences are OK
	
	die "Missing reference sequence for $chromosome" unless(exists $reference_href->{$chromosome});
	my $supposed_reference_sequence = uc(substr($reference_href->{$chromosome}, $firstPos_reference_0based, $lastPos_reference_0based - $firstPos_reference_0based + 1));
	
	my $alignment_reference_noGaps = $alignment_reference;
	$alignment_reference_noGaps =~ s/\-//g;
		
	unless($alignment_reference_noGaps eq $supposed_reference_sequence)
	{
		print "Reference mismatch for read $readID\n";
		#print "\t", "CIGAR: ", $inputAlignment->cigar_str, "\n";
		# print "\t", "Softclip remove: ", join("\t", $remove_softclipping_front, $remove_softclipping_back), "\n";			
		#print "\t", $alignment_reference, "\n";
		#print "\t", $alignment_read, "\n";
		print "\t", "firstPos_reference_0based: ", $firstPos_reference_0based, "\n";
		print "\t", "REF: ", substr($alignment_reference_noGaps, 0, 20), " .. ", substr($alignment_reference_noGaps, length($alignment_reference_noGaps)-10, 10), " ", length($alignment_reference_noGaps),  "\n";
		print "\t", "REF2: ", substr($reference_href->{$chromosome}, $firstPos_reference_0based, 20),  " ", length($supposed_reference_sequence), "\n";
		#print "\t", "ALG: ", $supposed_reference_sequence, "\n";
		print "\t", "strand: ", $strand, "\n";			
		print "\n";
		$warnings++;
		die if($warnings > 0);
		return undef;
	}


	my $raw_read_sequence = uc(getSequenceFromIndexedFasta($readsFasta, $readID));

	if(not $truncatedReads{$readID})
	{
		my $supposed_read_sequence = substr($raw_read_sequence, $firstPos_read_0based, $lastPos_read_0based - $firstPos_read_0based + 1);
		
		if($strand eq '-')
		{
			$supposed_read_sequence = reverseComplement($supposed_read_sequence);
		}
		
		my $alignment_read_noGaps = $alignment_read;
		$alignment_read_noGaps =~ s/\-//g;	
		
		unless($alignment_read_noGaps eq $supposed_read_sequence)
		{
			print "Sequence mismatch for read $readID\n";

			# print "\t", "CIGAR: ", $inputAlignment->cigar_str, "\n";
			# print "\t", "Softclip remove: ", join("\t", $remove_softclipping_front, $remove_softclipping_back), "\n";
			print "\t", "firstPos_read_0based: ", $firstPos_read_0based, "\n";
			print "\t", "lastPos_read_0based: ", $lastPos_read_0based, "\n";
			print "\t", "alignment_read_noGaps : ", length($alignment_read_noGaps), "\n";
			print "\t", "supposed_read_sequence: ", length($supposed_read_sequence), "\n\n";
			#print "\t", "alignment_ref : ", $alignment_reference, "\n";
			#print "\t", "alignment_read: ", $alignment_read, "\n";
			print "\t", "alignment_read_noGaps: ", substr($alignment_read_noGaps, 0, 10), "\n";
			print "\t", "supposed_read_sequence: ", substr($supposed_read_sequence, 0, 10), "\n";
			print "\n";
			$warnings++;
			die if($warnings > 0);
			return undef;
		}

		die "Mismatch query\n$alignment_read_noGaps\n$supposed_read_sequence" unless($alignment_read_noGaps eq $supposed_read_sequence);
	}
	
	my $n_matches = 0;
	my $n_mismatches = 0;
	my $n_insertions = 0;
	my $n_deletions = 0;
	my $n_gaps;

	my $alignment_length = length($alignment_reference);
	die Dumper("Lenght mismatch: " . length($alignment_read) . " v/s $alignment_length", $readID, $chromosome, $firstPos_reference_0based) unless(length($alignment_read) == $alignment_length);
	
	my $runningPos_query = 0;
	my $runningPos_reference = $firstPos_reference_0based - 1;
	my $runningPos_read = $firstPos_read_0based - 1;
	for(my $alignmentPosI = 0; $alignmentPosI < $alignment_length; $alignmentPosI++)
	{
		my $c_ref = substr($alignment_reference, $alignmentPosI, 1);
		my $c_query = substr($alignment_read, $alignmentPosI, 1);
		die if(($c_ref eq '-') and ($c_query eq '-'));
		if($c_ref eq '-')
		{
			$n_insertions++;			
			$runningPos_read++;
		}
		elsif($c_query eq '-')
		{
			$n_deletions++;
			$runningPos_reference++;
			
		}
		else
		{
			$runningPos_reference++;
			$runningPos_read++;
			if($c_ref eq $c_query)
			{
				$n_matches++;
			}
			else
			{
				$n_mismatches++;
			}
		}
	}
	
	die unless(defined $runningPos_read);
	die unless($runningPos_read == $lastPos_read_0based);
	die unless($runningPos_reference == $lastPos_reference_0based);
	
	$n_gaps = $n_insertions + $n_deletions;

	return {
		readID => $readID,
		chromosome => $chromosome,
		firstPos_reference => $firstPos_reference_0based,
		lastPos_reference => $lastPos_reference_0based,
		firstPos_read => $firstPos_read_0based,
		lastPos_read => $lastPos_read_0based,
		strand => $strand,
		n_matches => $n_matches,
		n_mismatches => $n_mismatches,
		n_gaps => $n_gaps,
		alignment_reference => $alignment_reference,
		alignment_read => $alignment_read,
		readLength => $readLength,
	};
}

sub indexFastaIfNecessary
{
	my $fasta = shift;
	die unless(defined $fasta);
	unless((-e $fasta . '.fai') and ((stat($fasta . '.fai'))[9] > (stat($fasta))[9]))
	{
		my $index_cmd = qq($samtools_path faidx $fasta);
		system($index_cmd) and die "Could not execute command: $index_cmd";
	}
}

sub getSequenceFromIndexedFasta
{
	my $fasta = shift;
	my $seqID = shift;
	die unless(defined $fasta);
	die unless(defined $seqID);

	die "FASTA $fasta not indexed" unless(-e $fasta . '.fai');
	my $extract_cmd = qq($samtools_path faidx $fasta "$seqID");

	open(SAMTOOLS, "$extract_cmd |") or die "Can't open samtools command $extract_cmd";

	my %R;
	my $currentSequence;
	while(<SAMTOOLS>)
	{
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
	
	close(SAMTOOLS);	
	
	unless(exists $R{$seqID})
	{
		die "Sequence $seqID could not be found in FASTA output";
	}
	
	return $R{$seqID};
}
 
 