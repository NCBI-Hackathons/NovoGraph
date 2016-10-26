#!/usr/bin/perl

use strict;
use Data::Dumper;
use Bio::DB::Sam;
use Getopt::Long;   


my $referenceFasta;
my $BAM;
my $outputFile = 'forAarti.txt';
my $paranoid = 1;
my $readsFasta;
my $lenientOrder = 1;

$referenceFasta = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\reference.fa';
$BAM = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\ribosome.bam';

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputFile:s' => \$outputFile,	
	'paranoid:s' => \$paranoid,	
	'readsFasta:s' => \$readsFasta,	
	'lenientOrder:s' => \$lenientOrder,	
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $reference_href;
my $reads_href;
if($paranoid)
{
	$reference_href = readFASTA($referenceFasta, 0);
	# $reads_href = readFASTA($readsFasta, 1);
}

my @out_headerfields = qw/readID chromosome firstPos_reference lastPos_reference firstPos_read lastPos_read strand n_matches n_mismatches n_gaps alignment_reference alignment_read/;
		
open(OUT, '>', $outputFile) or die "Cannot open $outputFile";
print OUT join("\t", @out_headerfields), "\n";

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);
my $iterator = $sam->features(-iterator=>1);

my $warnedAboutOrder;
my %processedReadIDs;
my $runningReadID;
while(my $alignment = $iterator->next_seq)
{
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
		print OUT join("\t", map {die "Field $_ undefined" unless(defined $read_href->{$_}); $read_href->{$_}} @out_headerfields), "\n";
	}
}


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
	while($cigar =~ /^(\d+)([MIDHS])/)
	{
		my $number = $1;
		my $action = $2;
		if(($action eq 'H') or ($action eq 'S'))
		{
			$firstPos_read += $number;
			$cigar =~ s/^(\d+)([MIDHS])//;
		}
		else
		{
			last;
		}
	}
		
	my ($ref, $matches, $query) = $inputAlignment->padded_alignment;
	unless(defined $ref)
	{
		warn "No reference sequence for $chromosome?";
		return undef;
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
	$n_gaps = $n_insertions + $n_deletions;
	
	if($paranoid)
	{
		die "Missing reference sequence for $chromosome" unless(exists $reference_href->{$chromosome});
		my $supposed_reference_sequence = substr($reference_href->{$chromosome}, $firstPos_reference, $lastPos_reference - $firstPos_reference + 1);
		
		# die "Missing read sequence for $readID" unless(exists $reads_href->{$readID});
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
				print $ref, "\n";
				print $query, "\n";
				print "\$inputAlignment->start: ", $inputAlignment->start, "\n";
				print "REF: ", $alignment_reference_noGaps, "\n";
				print "RE2: ", substr($reference_href->{$chromosome}, $firstPos_reference - 10, 20), "\n";
				print "ALG: ", $supposed_reference_sequence, "\n";
				print "\n\n\n";
				exit;
				die "Mismatch reference";
			}
			print "Ref ok!\n";
		}
		#print "Ref ok!\n";
		#die "Mismatch query" unless($alignment_read_noGaps eq $supposed_read_sequence);
	}
	
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