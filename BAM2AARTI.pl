#!/usr/bin/perl

use strict;
use Data::Dumper;
use Bio::DB::Sam;
use Getopt::Long;   


my $referenceFasta;
my $BAM;
my $outputFile = 'forAarti.txt';
my $paranoid = 0;
my $readsFasta;

$referenceFasta = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\reference.fa';
$BAM = 'C:\\Users\\AlexanderDilthey\\Desktop\\Temp\\Ribosomes\\ribosome.bam';

GetOptions (
	'referenceFasta:s' => \$referenceFasta, 
	'BAM:s' => \$BAM, 
	'outputFile:s' => \$outputFile,	
	'paranoid:s' => \$paranoid,	
	'readsFasta:s' => \$readsFasta,	
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
	$reads_href = readFASTA($readsFasta, 1);
}

my @out_headerfields = qw/readID firstPos_reference lastPos_reference firstPos_read lastPos_read strand n_matches n_mismatches n_gaps alignment_reference $alignment_reference alignment_read alignment_read/;
		
open(OUT, '>', $outputFile) or die "Cannot open $outputFile";
print OUT join("\t", @out_headerfields), "\n";

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);
my $iterator = $sam->features(-iterator=>1);

while(my $alignment = $iterator->next_seq)
{
	my $read_href = convertAlignmentToHash($alignment);	
	print OUT join("\t", map {die unless(defined $read_href->{$_}); $read_href->{$_}} @out_headerfields), "\n";
}


sub convertAlignmentToHash
{
	my $inputAlignment = shift;
    
	
	my $readID = $inputAlignment->query->name;
	my $firstPos_reference = $inputAlignment->start;
	my $lastPos_reference;
	my $firstPos_read;
	my $lastPos_read;
	my $strand = $inputAlignment->strand;
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
		
		
		
	}
	
	my ($ref, $matches, $query) = $inputAlignment->padded_alignment;
	$alignment_reference = $ref;
	$alignment_read = $query;
	
	my $alignment_length = length($ref);
	die unless(length($query) == $alignment_length);
	
	my $runningPos_query = 0;
	my $runningPos_reference = $firstPos_reference;
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
	
	return {
		readID => $readID,
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
