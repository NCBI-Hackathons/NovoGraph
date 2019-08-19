#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use Data::Dumper;
use Bio::DB::HTS;
use Getopt::Long;   
$| = 1;

## Usage:
## checkBAM_SVs_and_INDELs.pl --BAM <path to sorted, contigs BAM>
##                            --referenceFasta <path to reference FASTA>
##                            --readsFasta <path to FASTA of all contigs>
##
## Example command: 
## 	./checkBAM_SVs_and_INDELs.pl --BAM /data/projects/phillippy/projects/hackathon/shared/alignments/SevenGenomesPlusGRCh38Alts.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /data/projects/phillippy/projects/hackathon/shared/contigs/AllContigs.fa

my $referenceFasta;
my $BAM;
my $outputFile = 'fromBAM_lengthStatistics.txt';
my $readsFasta;
my $lenientOrder = 1;

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputFile:s' => \$outputFile,	
	'readsFasta:s' => \$readsFasta,	
	'lenientOrder:s' => \$lenientOrder
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --readsFasta" unless($readsFasta);
die "Please specify --output" unless($outputFile); 
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);
die "--readsFasta $readsFasta not existing" unless(-e $readsFasta);

my $outputFile2 = 'fromBAM_lengthStatistics.txt.individualAlignments';

my $reference_href;
my $reads_href;

my %read_reference_positions;
my %read_got_primary_alignment;


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
my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $BAM, -expand_flags => 1);
my $iterator = $sam->features(-iterator=>1);

my $warnedAboutOrder;

while(my $alignment = $iterator->next_seq)
{
	my $isPrimary  = (! $alignment->get_tag_values('NOT_PRIMARY'));
	# die if(not $isPrimary);
	# next unless($alignment->seq_id eq 'chr20');
		
	my $readID = $alignment->query->name;
		
	print "Collecting alignment $readID ... \n";
		
	my $read_href = convertAlignmentToHash($alignment);	
		
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
					warn "Read position $readPos in $readID covered by multiple alignments - existing reference value $read_reference_positions{$readID}[$readPos], want to set to $referencePos";
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



my %minusOne_perRead;
my $read_got_primary = 0;
my $read_no_primary = 0;
my %histograms;

open(OUT2, '>', $outputFile2) or die "Cannot open $outputFile2";

foreach my $readID (keys %read_reference_positions)
{
	print "Analyzing alignments for $readID ... \n";

	die unless(scalar(@{$read_reference_positions{$readID}}) == length($reads_href->{$readID}));
	for(my $posI = 0; $posI <= $#{$read_reference_positions{$readID}}; $posI++)
	{
		$read_reference_positions{$readID}[$posI] = -1 unless(defined $read_reference_positions{$readID}[$posI]);
	}
	
	if($read_got_primary_alignment{$readID})
	{
		$read_got_primary++;
		
		print OUT2 "Read ", $readID, "\n";
		my $first_definedPosition;
		my $last_definedPosition;
		for(my $readPos = 0; $readPos <= $#{$read_reference_positions{$readID}}; $readPos++)
		{
			if($read_reference_positions{$readID}[$readPos] != -1)
			{
				$first_definedPosition = $readPos unless(defined $first_definedPosition);
				$last_definedPosition = $readPos;
			}
			print OUT2 $readPos, "\t", $read_reference_positions{$readID}[$readPos], "\n";
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
}

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

my $outputFile_2 = $outputFile . ".minus1PerRead";
open(O, '>', $outputFile_2) or die "Cannot open $outputFile_2";
foreach my $rID (keys %minusOne_perRead)
{
	print O join("\t", $rID, $minusOne_perRead{$rID}), "\n";
}
close(O);

print "\n\nProduced output file $outputFile_2\n";

my $warnings = 0;
sub convertAlignmentToHash
{
	my $inputAlignment = shift;
	
	return undef if($inputAlignment->unmapped);
	
	my $readID = $inputAlignment->query->name;
	my $chromosome = $inputAlignment->seq_id;
	my $firstPos_reference = $inputAlignment->start - 1;
	die unless($firstPos_reference >= 0);
	
	my $lastPos_reference;
	my $firstPos_read = 0;
	my $lastPos_read;
	my $strand = $inputAlignment->strand;
	die unless(($strand == 1) or ($strand == -1));
	$strand = ($strand == 1) ? '+' : '-';
	
	my $n_matches = 0;
	my $n_mismatches = 0;
	my $n_insertions = 0;
	my $n_deletions = 0;
	my $n_gaps;
	my $alignment_reference;
	my $alignment_read;
	
	my $cigar  = $inputAlignment->cigar_str;
	my $remove_softclipping_front = 0;
	my $remove_softclipping_back = 0;	
	while($cigar =~ /^(\d+)([MIDHS])/)
	{
		my $number = $1;
		my $action = $2;
		$cigar =~ s/^(\d+)([MIDHS])//;		
		if(($action eq 'H') or ($action eq 'S'))
		{
			$firstPos_read += $number;
			if($action eq 'S')
			{
				$remove_softclipping_front += $number;
			}
		}
		else
		{
			last;
		}
	}
	
	{
		my $cigar  = $inputAlignment->cigar_str;
		while($cigar =~ /(\d+)([MIDHS])$/)
		{
			my $number = $1;
			my $action = $2;
			$cigar =~ s/(\d+)([MIDHS])$//;					
			if(($action eq 'H') or ($action eq 'S'))
			{
				if($action eq 'S')
				{
					$remove_softclipping_back += $number; 
				}
			}
			else
			{
				last;
			}
		}	
	}
	
	my ($ref, $matches, $query) = $inputAlignment->padded_alignment;
	unless(defined $ref)
	{
		warn "No reference sequence for $chromosome?";
		return undef;
	}
	
	if($remove_softclipping_front)
	{
		my $softclip_remove_ref = substr($ref, 0, $remove_softclipping_front);
		die "Weird softclipping for read" unless($softclip_remove_ref =~ /^\-+$/);
		substr($ref, 0, $remove_softclipping_front) = '';
		substr($query, 0, $remove_softclipping_front) = '';
	}

	if($remove_softclipping_back)
	{
		my $softclip_remove_ref = substr($ref, length($ref) - $remove_softclipping_back);
		die unless(length($softclip_remove_ref) == $remove_softclipping_back);
		
		die "Weird softclipping for read" unless($softclip_remove_ref =~ /^\-+$/);
		substr($ref, length($ref) - $remove_softclipping_back) = '';
		substr($query, length($query) - $remove_softclipping_back) = '';
	}
	
	$alignment_reference = $ref;
	$alignment_read = $query;
	
	my $alignment_length = length($ref);
	die Dumper("Lenght mismatch: " . length($query) . " v/s $alignment_length", $readID, $chromosome, $firstPos_reference) unless(length($query) == $alignment_length);
	
	my $runningPos_query = 0;
	my $runningPos_reference = $firstPos_reference - 1;
	my $runningPos_read = $firstPos_read;
	for(my $alignmentPosI = 0; $alignmentPosI < $alignment_length; $alignmentPosI++)
	{
		my $c_ref = substr($ref, $alignmentPosI, 1);
		my $c_query = substr($query, $alignmentPosI, 1);
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
		
		$lastPos_reference = $runningPos_reference;
		$lastPos_read = $runningPos_read;
	}
	
	die unless(defined $lastPos_reference);
	die unless(defined $runningPos_read);
	$lastPos_read--;
	die unless($lastPos_read >= $firstPos_read);
	
	$n_gaps = $n_insertions + $n_deletions;
	
	die "Missing reference sequence for $chromosome" unless(exists $reference_href->{$chromosome});
	my $supposed_reference_sequence = substr($reference_href->{$chromosome}, $firstPos_reference, $lastPos_reference - $firstPos_reference + 1);
	
	unless(exists $reads_href->{$readID})
	{
		warn "Missing read sequence for $readID";
		return undef;
	}
	
	my $read_sequence = $reads_href->{$readID};
	if($strand eq '-')
	{
		$read_sequence = reverseComplement($read_sequence);
	}
	my $supposed_read_sequence = substr($read_sequence, $firstPos_read, $lastPos_read - $firstPos_read + 1);
	
	my $alignment_reference_noGaps = $alignment_reference;
	$alignment_reference_noGaps =~ s/\-//g;
	
	my $alignment_read_noGaps = $alignment_read;
	$alignment_read_noGaps =~ s/\-//g;
	
	if(index($supposed_reference_sequence, "N") == -1)
	{
		unless($alignment_reference_noGaps eq $supposed_reference_sequence)
		{
			print "Reference mismatch for read $readID\n";
			print "\t", "CIGAR: ", $inputAlignment->cigar_str, "\n";
			print "\t", "Softclip remove: ", join("\t", $remove_softclipping_front, $remove_softclipping_back), "\n";			
			print "\t", $ref, "\n";
			print "\t", $query, "\n";
			print "\t", "\$inputAlignment->start: ", $inputAlignment->start, "\n";
			print "\t", "REF: ", $alignment_reference_noGaps, "\n";
			print "\t", "RE2: ", substr($reference_href->{$chromosome}, $firstPos_reference - 10, 20), "\n";
			print "\t", "ALG: ", $supposed_reference_sequence, "\n";
			print "\t", "strand: ", $strand, "\n";			
			print "\n";
			$warnings++;
			die if($warnings > 10);
			return undef;
		}

		unless($alignment_read_noGaps eq $supposed_read_sequence)
		{
			print "Sequence mismatch for read $readID\n";

			print "\t", "CIGAR: ", $inputAlignment->cigar_str, "\n";
			print "\t", "Softclip remove: ", join("\t", $remove_softclipping_front, $remove_softclipping_back), "\n";
			print "\t", "lastPos_read: ", $lastPos_read, "\n";
			print "\t", "firstPos_read: ", $firstPos_read, "\n";
			print "\t", "alignment_read_noGaps : ", $alignment_read_noGaps, "\n";
			print "\t", "supposed_read_sequence: ", $supposed_read_sequence, "\n\n";
			print "\t", "alignment_ref : ", $ref, "\n";
			print "\t", "alignment_query: ", $query, "\n";
			print "\t", "strand: ", $strand, "\n";
			print "\n";
			$warnings;;
			die if($warnings > 10);
			return undef;
		}
	}
	
	die "Mismatch query" unless($alignment_read_noGaps eq $supposed_read_sequence);
	
	return {
		readID => $readID,
		chromosome => $chromosome,
		firstPos_reference => $firstPos_reference,
		lastPos_reference => $lastPos_reference,
		firstPos_read => $firstPos_read,
		lastPos_read => $lastPos_read,
		strand => $strand,
		n_matches => $n_matches,
		n_mismatches => $n_mismatches,
		n_gaps => $n_gaps,
		alignment_reference => $alignment_reference,
		alignment_read => $alignment_read,
	};
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
