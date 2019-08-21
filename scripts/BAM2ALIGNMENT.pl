#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use Data::Dumper;
use Getopt::Long;   
$| = 1;

## Usage:
## BAM2ALIGNMENT.pl --BAM <path to sorted contigs BAM> 
##                  --referenceFasta <path to reference FASTA> 
##                  --readsFasta <path to contigs FASTA> 
##                  --outputFile <path to text file>
##
## Example command:
## ./BAM2ALIGNMENT.pl --BAM /home/data/alignments/SevenGenomesPlusGRCh38Alts.bam 
##                    --referenceFasta /home/data/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa 
##                    --readsFasta /home/data/contigs/AllContigs.fa 
##                    --outputFile /intermediate_files/AlignmentInput.txt


my $referenceFasta;
my $BAM;
my $outputFile;
my $readsFasta;
my $bin_sam2alignment;
my $samtools_path;

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputFile:s' => \$outputFile,	
	'readsFasta:s' => \$readsFasta,	
	'sam2alignment_executable:s' => \$bin_sam2alignment,
	'samtools_path:s' => \$samtools_path,		
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

my $fn_BAM_SAM = $BAM . '.extract.SAM';
my $cmd_extract_SAM = qq($samtools_path view $BAM > $fn_BAM_SAM);
system($cmd_extract_SAM) and die "Cannot execute command: $cmd_extract_SAM"; 


my $fn_BAM_alignments = $BAM . '.extract.alignments';
my $cmd_doConvert = qq($bin_sam2alignment $fn_BAM_SAM $referenceFasta > $fn_BAM_alignments);
system($cmd_doConvert) and die "Cannot execute command: $cmd_doConvert";

print $fn_BAM_SAM, "\n", $fn_BAM_alignments, "\n";

print "Read $referenceFasta\n";
my $reference_href = readFASTA($referenceFasta, 0);
print "\tdone.\n";

indexFastaIfNecessary($readsFasta);

my @out_headerfields = qw/readID chromosome firstPos_reference lastPos_reference firstPos_read lastPos_read strand n_matches n_mismatches n_gaps alignment_reference alignment_read completeReadSequence_plus completeReadSequence_minus/;

my $headerFn = $outputFile . '.header';
open(HEADER, '>', $headerFn) or die "Cannot open $headerFn";
print HEADER join("\t", @out_headerfields), "\n";
close(HEADER);

open(OUT, '>', $outputFile) or die "Cannot open $outputFile";

my %processedReadIDs;
my $runningReadID;
my %printed_complete_sequence;
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
	
	if($runningReadID ne $readID)
	{
		if($runningReadID)
		{		
			die "Ordering constrains violated - we want a BAM ordered by read ID - violating read ID $readID" if($processedReadIDs{$readID});
			$processedReadIDs{$runningReadID}++;
		}
		$runningReadID = $readID;
	}
	
	my $read_href = convertAlignmentToHash($readID, [$referenceSequenceID, $alignmentStart_1based, $alignmentStop_1based], [$alignmentStartInRead_1based, $alignmentStopInRead_1based], [$aligned_reference, $aligned_read]);	
	
	if($read_href)
	{
		if($printed_complete_sequence{$readID})
		{
			$read_href->{completeReadSequence_plus} = '';
			$read_href->{completeReadSequence_minus} = '';
		}
		else
		{
			my $sequence_read = uc(getSequenceFromIndexedFasta($readsFasta, $readID));
			$read_href->{completeReadSequence_plus} = $sequence_read;
			$read_href->{completeReadSequence_minus} = reverseComplement($sequence_read);
			$printed_complete_sequence{$readID}++;
		}
		
		print OUT join("\t", map {die "Field $_ undefined" unless(defined $read_href->{$_}); $read_href->{$_}} @out_headerfields), "\n";
	}	
}

# It would seem that the sorting steps are not really necessary

my $sorted_outputFile = $outputFile.'.sorted';

my $sort_cmd = qq(sort $outputFile > $sorted_outputFile);
if(system($sort_cmd))
{
	die "Failed during command: $sort_cmd";
}
unless(-e $sorted_outputFile)
{
	die "Expected output file $sorted_outputFile not existing";
}

my $combined_outputFile = $outputFile . '.sortedWithHeader';
my $combine_cmd = qq(cat $headerFn $sorted_outputFile > $combined_outputFile);
if(system($combine_cmd))
{
	die "Failed during command: $combine_cmd";
}

print "\n\nProduced output file $combined_outputFile\n";

## remove extraneous files, e.g.
## '../intermediate_files/AlignmentInput.txt'
## '../intermediate_files/AlignmentInput.txt.header'
## '../intermediate_files/AlignmentInput.txt.sorted'
## only 'combined_outputFile' should remain

unlink($headerFn);
unlink($outputFile);
unlink($sorted_outputFile);


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

	my $raw_read_sequence = uc(getSequenceFromIndexedFasta($readsFasta, $readID));
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

	my $firstPos_read_0based_forReturn = $firstPos_read_0based;
	my $lastPos_read_0based_forReturn = $lastPos_read_0based;
	if($strand eq '-')
	{ 
		$firstPos_read_0based_forReturn = length($raw_read_sequence) - $lastPos_read_0based - 1;
		$lastPos_read_0based_forReturn = length($raw_read_sequence) - $firstPos_read_0based - 1;
	} 
	else
	{
		die unless($strand eq '+');
	}
	return {
		readID => $readID,
		chromosome => $chromosome,
		firstPos_reference => $firstPos_reference_0based,
		lastPos_reference => $lastPos_reference_0based,
		firstPos_read => $firstPos_read_0based_forReturn,
		lastPos_read => $lastPos_read_0based_forReturn,
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
	$kMer =~ tr/ACGTacgt/TGCAtgca/;
	return scalar(reverse($kMer));
}

sub indexFastaIfNecessary
{
	my $fasta = shift;
	die unless(defined $fasta);

	unless(-e $fasta . '.fai')
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
 
 