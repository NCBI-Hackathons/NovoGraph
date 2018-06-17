#!/usr/bin/perl

use strict;
use Data::Dumper;
use Bio::DB::Sam;
use Getopt::Long;   
$| = 1;

# Example command:
# To check correctness of INPUT for mafft:
# 	./BAM2AARTI.pl --BAM /home/data/alignments/SevenGenomesPlusGRCh38Alts.bam --referenceFasta /home/data/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /home/data/contigs/AllContigs.fa

my $referenceFasta;
my $BAM;
my $outputFile = 'forAarti.txt';
my $readsFasta;
my $lenientOrder = 1;

# $referenceFasta = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\reference.fa';
# $BAM = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\ribosome.bam';

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

print "Read $referenceFasta\n";
my $reference_href = readFASTA($referenceFasta, 0);
print "\tdone.\n";

print "Read $readsFasta\n";
my $reads_href = readFASTA($readsFasta, 0);
print "\tdone.\n";

my @out_headerfields = qw/readID chromosome firstPos_reference lastPos_reference firstPos_read lastPos_read strand n_matches n_mismatches n_gaps alignment_reference alignment_read completeReadSequence_plus completeReadSequence_minus/;

my $headerFn = $outputFile . '.header';
open(HEADER, '>', $headerFn) or die "Cannot open $headerFn";
print HEADER join("\t", @out_headerfields), "\n";
close(HEADER);

open(OUT, '>', $outputFile) or die "Cannot open $outputFile";

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);
my $iterator = $sam->features(-iterator=>1);

my $warnedAboutOrder;
my %processedReadIDs;
my $runningReadID;
my %printed_complete_sequence;
while(my $alignment = $iterator->next_seq)
{
	if($alignment->seq_id ne 'chr1')
	{
		#warn "For testing purposes, stop after chr1"; # todo
		#last;
	}
	
	my $readID = $alignment->query->name;
	if($runningReadID ne $readID)
	{
		if($runningReadID)
		{		
			if($lenientOrder)
			{
				if(not $warnedAboutOrder)
				{
					die "Ordering constrains violated - we want a BAM ordered by read ID - violating read ID $readID" if($processedReadIDs{$readID});
					$warnedAboutOrder++;
				}
			}
			else
			{
				die "Ordering constrains violated - we want a BAM ordered by read ID - violating read ID $readID" if($processedReadIDs{$readID});
			}
			
			$processedReadIDs{$runningReadID}++;
		}
		$runningReadID = $readID;
	}
	my $read_href = convertAlignmentToHash($alignment);	
	
	if($read_href)
	{
		if($printed_complete_sequence{$readID})
		{
			$read_href->{completeReadSequence_plus} = '';
			$read_href->{completeReadSequence_minus} = '';
		}
		else
		{
			die unless($reads_href->{$readID});
			$read_href->{completeReadSequence_plus} = $reads_href->{$readID};
			$read_href->{completeReadSequence_minus} = reverseComplement($reads_href->{$readID});
			$printed_complete_sequence{$readID}++;
		}
		
		print OUT join("\t", map {die "Field $_ undefined" unless(defined $read_href->{$_}); $read_href->{$_}} @out_headerfields), "\n";
	}
}

my $sorted_outputFile = $outputFile.'.sorted';
#die "Implement this properly - first line issue!";
# sed '1d' forAarti > AartiInput.forSort

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
	die unless(length($query) == $alignment_length);
	
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

# my $bam          = Bio::DB::Bam->open($BAM);
# my $header       = $bam->header;
# my $target_count = $header->n_targets;
# my $target_names = $header->target_name;
# while (my $align = $bam->read1) {
	# my $seqid     = $target_names->[$align->tid];
	# my $start     = $align->pos+1;
	# my $end       = $align->calend;
	# my $cigar     = $align->cigar_str;
	# my $readID = $align->query->name;
	
	# die ref($align);
	# die $align->padded_alignment;
	# die $readID;
# }
 
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
