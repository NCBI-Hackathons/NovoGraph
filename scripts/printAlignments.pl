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
# To check correctness of INPUT for mafft:
# 	./printAlignments.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.fa
# 	./printAlignments.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_4.bam.sorted.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_4.fa



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
#my $reference_href = readFASTA($referenceFasta);
#print "\t...done.\n";

foreach my $referenceSequenceID (@sequence_ids)
{
	#print "Processing $referenceSequenceID .. \n";
	#die "Sequence ID $referenceSequenceID not defined in $referenceFasta" unless(exists $reference_href->{$referenceSequenceID});
	
	
	my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);

	my $n_alignments = 0;
	my %alignments_starting_at;
	while(my $alignment = $alignment_iterator->next_seq)
	{
		$n_alignments++;		

		my $alignment_start_pos = $alignment->start;
		my ($ref,$matches,$query) = $alignment->padded_alignment;
		
		my $query_start = $alignment->query->start;     
		my $query_end   = $alignment->query->end;
	
		my $query_end_2 = calc_final_refpos($alignment_start_pos, $alignment->cigar_str);
		print $alignment->query->name, " from ", $alignment_start_pos, " to ", $query_end_2, "\n";
		#print $alignment->cigar_str, "\n";
		print $ref, "\n";
		print $query, "\n";
	}
}	



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
			# last if(($currentSequence ne 'chr1') or ($currentSequence ne 'ref')); # todo remove
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

sub calc_final_refpos {
    my $pos = shift;
    my $orig_cigar = shift;
    my $cigar = $orig_cigar;

	my $sawMatch = 0;
    my $thispos = $pos - 1;
    while ($cigar) {
        my ($nbases, $nextop) = ($cigar =~ s/^(\d+)([MIDSHP])//) ? ($1, $2) : (undef, undef);
        if (!$nbases) {
            die "Unparsable cigar string $orig_cigar!\n";
        }

        if ($nextop eq 'M') {
            $thispos += $nbases;
			$sawMatch = 1;
        }
        if ($nextop eq 'D') {
			if($sawMatch)
			{
				$thispos += $nbases;
			}
        }
    }

    return $thispos;
}


