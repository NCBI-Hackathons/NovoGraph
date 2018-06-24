#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max/;
use List::MoreUtils qw/mesh/;
use Bio::DB::Sam;

$| = 1;

# Example command:
# 	./validate_BAM_MSA.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/combined_sorted.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa


my $BAM;
my $referenceFasta;

GetOptions (
	'BAM:s' => \$BAM, 
	'referenceFasta:s' => \$referenceFasta, 
);
	
die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);

die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my @sequence_ids = $sam->seq_ids();

#print "Reading $referenceFasta\n";
my $reference_href = readFASTA($referenceFasta);
#print "\t...done.\n";

foreach my $referenceSequenceID (@sequence_ids) # todo reinstate
{
	print "Processing $referenceSequenceID .. \n";
	die "Sequence ID $referenceSequenceID not defined in $referenceFasta" unless(exists $reference_href->{$referenceSequenceID});
	
	my $l_ref_sequence = length($reference_href->{$referenceSequenceID});
	my @gap_structure;
	$#gap_structure = ($l_ref_sequence - 1);
	#print "Set...";
	for(my $i = 0; $i < $l_ref_sequence; $i++)
	{
		$gap_structure[$i] = -1;
	}
	#print " .. done.\n";
	die unless(scalar(@gap_structure) == $l_ref_sequence);
	
	my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);

	my $n_alignments = 0;
	my %alignments_starting_at;
	while(my $alignment = $alignment_iterator->next_seq)
	{
		$n_alignments++;		

		my $alignment_start_pos = $alignment->start;
		my ($ref,$matches,$query) = $alignment->padded_alignment;
		
		my $gaps_left_side = 0;
		my $gaps_right_side = 0;
		my $firstMatch = -1;
		for(my $i = 0; $i < length($ref); $i++)
		{
			if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*') and (substr($query, $i, 1) ne '-') and (substr($query, $i, 1) ne '*'))
			{
				$firstMatch = $i;
				last; # first match
			}
			
			my $isDeletion = (((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*')) and ((substr($query, $i, 1) eq '-') or (substr($query, $i, 1) eq '*')));
			my $isInsertion = (((substr($query, $i, 1) ne '-') and (substr($query, $i, 1) ne '*')) and ((substr($ref, $i, 1) eq '-') or (substr($ref, $i, 1) eq '*')));
			die if($isDeletion and $isInsertion);
			
			if($isDeletion)
			{
				# warn Dumper("Deletions before start?", $ref, $query, $alignment->query->name, $alignment->cigar_str);
			}
			if($isInsertion)
			{
				# warn Dumper("One of the insertions!", $ref, $query, $alignment->query->name, $alignment->cigar_str);
			}
			
		}
		die if($firstMatch == -1);
		
		$ref = substr($ref, $firstMatch);
		$query = substr($query, $firstMatch);
		
		
		for(my $i = 0; $i < length($ref); $i++)
		{
			if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*'))
			{
				last;
			}
			else
			{
				$gaps_left_side++;
			}
		}
		
		for(my $i = length($ref) - 1; $i >= 0; $i--)
		{
			if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*'))
			{
				last;
			}
			else
			{
				$gaps_right_side++;
			}
		}
		
		my $starstar_left_side = 0;
		my $starstar_right_side = 0;
		for(my $i = 0; $i < length($ref); $i++)
		{
			if((substr($ref, $i, 1) eq '*') and (substr($query, $i, 1) eq '*'))
			{
				$starstar_left_side++;
			}
			else
			{
				last;
			}
		}
		
		for(my $i = length($ref) - 1; $i >= 0; $i--)
		{
			if((substr($ref, $i, 1) eq '*') and (substr($query, $i, 1) eq '*'))
			{
				$starstar_right_side++;
			}
			else
			{
				last;
			}
		}
		
		
		print "Sequence ", $alignment->query->name, " remove $gaps_left_side left and $gaps_right_side right \n";
		$ref = substr($ref, $gaps_left_side, length($ref) - $gaps_left_side - $gaps_right_side);
		$query = substr($query, $gaps_left_side, length($query) - $gaps_left_side - $gaps_right_side);		
		#$ref = substr($ref, 0, length($ref) - $gaps_right_side);
		#$query = substr($query, 0, length($query) - $gaps_right_side);		
		die unless(length($ref) == length($query));
		
		push(@{$alignments_starting_at{$alignment_start_pos}}, [$ref, $query, $alignment->query->name]);
		
		$n_alignments++;
		
		my $start_pos = $alignment_start_pos-2;
		
		my $ref_pos = $start_pos;
		my $running_gaps = 0;
		for(my $i = 0; $i < length($ref); $i++)
		{
			my $c_ref = substr($ref, $i, 1);
			
			if(($c_ref eq '-') or ($c_ref eq '*'))
			{
				$running_gaps++;
			}
			else
			{
				if(not defined $gap_structure[$ref_pos])
				{
					die "Position $ref_pos not defined - length $l_ref_sequence";
				}
				if(($ref_pos >= 0) and ($ref_pos != $start_pos))
				{				
					if($gap_structure[$ref_pos] == -1)
					{
						$gap_structure[$ref_pos] = $running_gaps;	
					}
					else
					{
						die Dumper("Gap structure mismatch at position $referenceSequenceID $ref_pos - this is alignment $n_alignments, have existing value $gap_structure[$ref_pos], want to set $running_gaps", $alignment_start_pos-1, $alignment->query->name, $alignment->cigar_str, $ref, $query) unless($gap_structure[$ref_pos] == $running_gaps);
					}
				}
				$ref_pos++;
				$running_gaps = 0;
				
			}
		}
	}
}	

print "\n\nOK\n\n";

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


