#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
# use Bio::DB::HTS;
use Getopt::Long;   
$| = 1;

## Usage:
## checkBAM_SVs_and_INDELs.pl --BAM <path to sorted by name, contigs BAM>
##                            --referenceFasta <path to reference FASTA>
##                            --readsFasta <path to FASTA of all contigs>
##
## Example command: 
## 	./checkBAM_SVs_and_INDELs_linear.pl --BAM /data/projects/phillippy/projects/hackathon/shared/alignments/SevenGenomesPlusGRCh38Alts.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /data/projects/phillippy/projects/hackathon/shared/contigs/AllContigs.fa

my $referenceFasta;
my $BAM;
my $outputFile = 'fromBAM_lengthStatistics.txt';
my $readsFasta;
my $lenientOrder = 1;
my $bin_sam2alignment;
my $samtools_path;
my $printDetailedAlignmentData;

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputFile:s' => \$outputFile,	
	'readsFasta:s' => \$readsFasta,	
	'lenientOrder:s' => \$lenientOrder,
	'sam2alignment_executable:s' => \$bin_sam2alignment,
	'samtools_path:s' => \$samtools_path,	
    'printDetailedAlignmentData:s' => \$printDetailedAlignmentData,
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --readsFasta" unless($readsFasta);
die "Please specify --output" unless($outputFile); 
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);
die "--readsFasta $readsFasta not existing" unless(-e $readsFasta);
die "--sam2alignment_executable $bin_sam2alignment not present; Please run 'make' in the directory /src." unless(-e $bin_sam2alignment);

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

my $outputFile2 = 'fromBAM_lengthStatistics.txt.individualAlignments';

my $reference_href; 
my $reads_href;

print "Read $referenceFasta\n";
$reference_href = readFASTA($referenceFasta, 0);
print "\tdone.\n";

print "Read $readsFasta\n";
$reads_href = readFASTA($readsFasta, 0);
print "\tdone.\n";

foreach my $readID (keys %$reads_href)
{
	$reads_href->{$readID} = uc($reads_href->{$readID});
}

my $fn_BAM_SAM = $BAM . '.extract.SAM';
my $cmd_extract_SAM = qq($samtools_path view $BAM > $fn_BAM_SAM);
system($cmd_extract_SAM) and die "Cannot execute command: $cmd_extract_SAM"; 


my $fn_BAM_alignments = $BAM . '.extract.alignments';
my $cmd_doConvert = qq($bin_sam2alignment $fn_BAM_SAM $referenceFasta > $fn_BAM_alignments);
system($cmd_doConvert) and die "Cannot execute command: $cmd_doConvert";

print $fn_BAM_SAM, "\n", $fn_BAM_alignments, "\n";

indexFastaIfNecessary($readsFasta);
foreach my $readID (keys %$reads_href)
{
	my $seqFromFasta = uc(getSequenceFromIndexedFasta($readsFasta, $readID));
	die Dumper("Mismatch for read $readID") unless($reads_href->{$readID} eq $seqFromFasta); 
}

#debugTest/_combined.cram.extract.alignments system($cmd_extract_SAM) and die "Cannot execute command: $cmd_extract_SAM";


my $runningReadID;
my %read_reference_positions;
my %read_got_primary_alignment;
my %processed_readID;

my %minusOne_perRead;
my $read_got_primary = 0;
my $read_no_primary = 0;
my %histograms;

if ($printDetailedAlignmentData)
{
    open(OUT2, '>', $outputFile2) or die "Cannot open $outputFile2";
}


my $process_collected_read_data = sub {
	die "Weird number of collected reads: ". scalar(keys %read_reference_positions) unless(scalar(keys %read_reference_positions) == 1);
	my $readID = (keys %read_reference_positions)[0];
	print "Analyzing alignments for $readID ... - have " . scalar(@{$read_reference_positions{$readID}}) . " read positions and " . scalar(keys %read_reference_positions) . " entries.\n";

	die unless(scalar(@{$read_reference_positions{$readID}}) == length($reads_href->{$readID}));
	for(my $posI = 0; $posI <= $#{$read_reference_positions{$readID}}; $posI++)
	{
		$read_reference_positions{$readID}[$posI] = -1 unless(defined $read_reference_positions{$readID}[$posI]);
	}
	
	if($read_got_primary_alignment{$readID})
	{
		$read_got_primary++;
		
		print OUT2 "Read ", $readID, "\n" if ($printDetailedAlignmentData);
		my $first_definedPosition;
		my $last_definedPosition;
		for(my $readPos = 0; $readPos <= $#{$read_reference_positions{$readID}}; $readPos++)
		{
			if($read_reference_positions{$readID}[$readPos] != -1)
			{
				$first_definedPosition = $readPos unless(defined $first_definedPosition);
				$last_definedPosition = $readPos;
			}
			print OUT2 $readPos, "\t", $read_reference_positions{$readID}[$readPos], "\n" if ($printDetailedAlignmentData);
		}
		
		my $padding_front;
		my $padding_end;
		if(defined $first_definedPosition)
		{
			my %local_histograms;
			my $allConsistent = 1;
			
			$padding_front = $first_definedPosition;
			$padding_end = $#{$read_reference_positions{$readID}} - $last_definedPosition;
			
			my $runningGapLength = 0;
			my $lastRefPos;
			for(my $readPos = $padding_front; $readPos <= $last_definedPosition; $readPos++)
			{
				if($read_reference_positions{$readID}[$readPos] != -1)
				{
					my $diff_refPos = (defined $lastRefPos) ? ($read_reference_positions{$readID}[$readPos] - $lastRefPos) : 1;
					my $deletion_relative_to_reference = $diff_refPos - 1;
					
					$allConsistent = 0 if($deletion_relative_to_reference < 0);
					$runningGapLength -= $deletion_relative_to_reference;
					
					if($runningGapLength)
					{	
						die unless($lastRefPos >= 0);
						$local_histograms{$runningGapLength}++;
					}
					
					$runningGapLength = 0;
					$lastRefPos = $read_reference_positions{$readID}[$readPos];					
				}
				else
				{
					$runningGapLength++;
				}
				
			}		
			
			die unless($runningGapLength == 0);
			
			$histograms{consistent}{$allConsistent}++;
			
			foreach my $k (keys %local_histograms)
			{
				if($allConsistent)
				{
					$histograms{INDELs}{$k} += $local_histograms{$k};
					if($k == -1)
					{
						$minusOne_perRead{$readID} += $local_histograms{$k};				
					}
				}
			}
		
		}
		else
		{
			$padding_front = $#{$read_reference_positions{$readID}}+1;
			$padding_end = 0;
		}
		
		$histograms{frontPadding}{$padding_front}++;
		$histograms{endPadding}{$padding_end}++;	
		my $combinedPadding = $padding_front + $padding_end;
		my $combinedPadding_perc = int(($combinedPadding / ($#{$read_reference_positions{$readID}}+1)) * 100);
		$histograms{combinedPaddingPerc}{$combinedPadding_perc}++;		
		
	}
	else
	{
		$read_no_primary++;
	}
};

open(ALIGNMENTS, '<', $fn_BAM_alignments) or die "Cannot open $fn_BAM_alignments";
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
	
	die unless($header =~ /secondary=([10])/);
	my $isPrimary  = (! $1);
	print $isPrimary, "\n";
	
	if((defined $runningReadID) and ($runningReadID ne $readID))
	{
		$process_collected_read_data->();
		$processed_readID{$runningReadID}++;
		%read_reference_positions = ();
	}
	$runningReadID = $readID;
	if(exists $processed_readID{$readID})
	{
		die "Warning - $readID has already been processed - is this input BAM sorted by read ID?";
	}
	print "Collecting alignment $readID ... - running entries " . scalar(keys %read_reference_positions) . " - processed " . scalar(keys %processed_readID) . "\n";
		
	my $read_href = convertAlignmentToHash($readID, [$referenceSequenceID, $alignmentStart_1based, $alignmentStop_1based], [$alignmentStartInRead_1based, $alignmentStopInRead_1based], [$aligned_reference, $aligned_read]);	

	(my $read_noGaps = $aligned_read) =~ s/-//g;
	
	
	if($read_href)
	{
		if($isPrimary)
		{
			$read_got_primary_alignment{$readID} = 1;
		}
		
		if(not defined $read_reference_positions{$readID})
		{
			$read_reference_positions{$readID} = [];
			$#{$read_reference_positions{$readID}} = (length($reads_href->{$readID})-1);
		}
			
		die unless($#{$read_reference_positions{$readID}} == (length($reads_href->{$readID})-1));	
		
		my $runningRefPos = $read_href->{firstPos_reference};
		my $runningReadPos = $read_href->{firstPos_read};
		for(my $alignmentPosI = 0; $alignmentPosI < length($read_href->{alignment_reference}); $alignmentPosI++)
		{
			my $c_ref = substr($read_href->{alignment_reference}, $alignmentPosI, 1);
			my $c_read = substr($read_href->{alignment_read}, $alignmentPosI, 1);
			
			my $referencePos;
			my $readPos;
			if($c_ref eq '-')
			{
				$referencePos = -1;
			}
			else
			{
				$referencePos = $runningRefPos;
			}
			
			if($c_read eq '-')
			{
				$readPos = -1;
			}
			else
			{
				$readPos = $runningReadPos; 
			}		
			
			if($readPos != -1)
			{
				die unless($readPos >= 0); 
				die unless($readPos <= $#{$read_reference_positions{$readID}});
				if(defined $read_reference_positions{$readID}[$readPos]) 
				{
					#warn "Read position $readPos in $readID covered by multiple alignments - existing reference value $read_reference_positions{$readID}[$readPos], want to set to $referencePos";
				}
				$read_reference_positions{$readID}[$readPos] = $referencePos;
			}
			
			if($c_ref ne '-')
			{
				$runningRefPos++;
			}	
			if($c_read ne '-')
			{
				$runningReadPos++;
			}	
		}
	}
}
close(ALIGNMENTS);
unlink($fn_BAM_alignments) if not $printDetailedAlignmentData;
unlink($fn_BAM_SAM) if not $printDetailedAlignmentData;

$process_collected_read_data->();

$histograms{primary}{1} = $read_got_primary;
$histograms{primary}{0} = $read_no_primary;


open(OUT, '>', $outputFile) or die "Cannot open $outputFile";

foreach my $category (keys %histograms)
{
	foreach my $l (keys %{$histograms{$category}})
	{
		print OUT join("\t", $category, $l, $histograms{$category}{$l}), "\n";
	}
}
close(OUT);

print "\n\nProduced output file $outputFile (only from reads that have primary information)\n";
print "\t", "read_got_primary", ": ", $read_got_primary, "\n";
print "\t", "read_no_primary", ": ", $read_no_primary, "\n";

my $outputFile_3 = $outputFile . ".minus1PerRead";
open(O, '>', $outputFile_3) or die "Cannot open $outputFile_3";
foreach my $rID (keys %minusOne_perRead)
{
	print O join("\t", $rID, $minusOne_perRead{$rID}), "\n";
}
close(O);

print "\n\nProduced output file $outputFile_3\n";

my $warnings = 0;
sub convertAlignmentToHash
{
	my $readID = shift;
	my $referenceCoordinates_aref = shift;
	my $readCoordinates_aref = shift;
	my $alignment = shift;
	
	die unless(length($alignment->[0]) eq length($alignment->[1]));
	
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


	unless(exists $reads_href->{$readID})
	{
		die "Missing read sequence for $readID";
		return undef;
	}
	

	my $raw_read_sequence = $reads_href->{$readID};
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
		print "\t", "strand: ", $strand, "\n";
		print "\n";
		$warnings++;
		die if($warnings > 0);
		return undef;
	}

	die "Mismatch query\n$alignment_read_noGaps\n$supposed_read_sequence" unless($alignment_read_noGaps eq $supposed_read_sequence);

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
 
sub readFASTA
{
	my $file = shift;	
	my $keepCompleteIdentier = shift;
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if($_ =~ /\r/)
		{
			die "It seems that $file contains Windows-style line endings. It is highly recommended to remove these. I will abort, comment out this line if you want to continue regardless.";
		}
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
			unless($keepCompleteIdentier)
			{
				$currentSequence =~ s/\s.+//;
			}
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
