#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max/;
use List::MoreUtils qw/mesh/;

$| = 1;

## Usage:
## FIND_GLOBAL_ALIGNMENTS.pl --alignmentsFile <path to output of BAM2ALIGNMENT.pl, with extension *.sortedWithHeader>
##                           --referenceFasta <path to reference FASTA>
##                           --outputFile <name of output BAM>
##                           --outputTruncatedReads <name of text outfile, e.g. 'truncatedReads'> 
##                           --outputReadLengths <name of text outfile, e.g. 'postGlobalAlignment_readLengths'>
##                           --CIGARscript_path <path to script dealWithTooManyCIGAROperations.pl>
##
## Example command:
## ./FIND_GLOBAL_ALIGNMENTS.pl --alignmentsFile /intermediate_files/AartiInput.sortedWithHeader 
##                             --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa 
##                             --outputFile /intermediate_files/forMAFFT.bam 
##                             --outputTruncatedReads intermediate_files/truncatedReads 
##                             --outputReadLengths /intermediate_files/postGlobalAlignment_readLengths
##                             --CIGARscript_path dealWithTooManyCIGAROperations.pl
##

my $alignmentsFile;
my $referenceFasta;
my $outputFile = 'output.sam';
my $lenientOrder = 1;
my $outputTruncatedReads;
my $outputReadLengths;
my $endsFree_reference = 1;
# my $CIGARscript_path; 
my $samtools_path;

my $S_match = 1;
my $S_mismatch = -1;
my $S_gap = -1;

GetOptions (
	'alignmentsFile:s' => \$alignmentsFile, 
	'referenceFasta:s' => \$referenceFasta, 
	'outputFile:s' => \$outputFile,	
	'outputTruncatedReads:s' => \$outputTruncatedReads,
	'outputReadLengths:s' => \$outputReadLengths,
#	'CIGARscript_path:s' => \$CIGARscript_path,
	'samtools_path:s' => \$samtools_path,	
);

die "Please specify --alignmentsFile" unless($alignmentsFile);
die "--alignmentsFile $alignmentsFile not existing" unless(-e $alignmentsFile);

die "Please specify --referenceFasta" unless($referenceFasta);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

die "Please specify --outputFile" unless($outputFile);

die "Please specify --outputTruncatedReads" unless($outputTruncatedReads);
die "Please specify --outputReadLengths" unless($outputReadLengths);

# die "Please specify path to script dealWithTooManyCIGAROperations.pl --CIGARscript_path" unless($CIGARscript_path);

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


print "Read $referenceFasta\n";
my $reference_href = readFASTA($referenceFasta, 0);
print "\tdone.\n";

# my $outputFile_sam_unsorted = $outputFile . '.sam.unfiltered';
my $outputFile_sam_unsorted = $outputFile . '.unsorted';

open(SAMOUTPUT, '>', $outputFile_sam_unsorted) or die "Cannot open $outputFile_sam_unsorted";
print SAMOUTPUT "\@HD\tVN:1.5", "\n";
foreach my $refChromosome (keys %$reference_href)
{
	print SAMOUTPUT "\@SQ\tSN:${refChromosome}\tLN:", length($reference_href->{$refChromosome}),"\n";
}

open(INPUT, '<', $alignmentsFile) or die "$alignmentsFile not existing";
my $alignments_headerLine = <INPUT>;
chomp($alignments_headerLine);
my @alignments_headerFields = split(/\t/, $alignments_headerLine);
die unless($alignments_headerFields[0] eq 'readID');

my $n_output_alignments = 0;
my $n_chains_sum = 0;
my $n_alignments_leftGapsRemoved = 0;
my $n_alignments_rightGapsRemoved = 0;
my $n_alignments_leftAndRightGapsRemoved = 0;

open(TRUNCATED, ">", $outputTruncatedReads) or die "Cannot open $outputTruncatedReads";
open(LENGTHS, ">", $outputReadLengths) or die "Cannot open $outputReadLengths";

my %processed_readIDs;
my @lines_current_read;
my $currentReadID = '';
my $processReadLines = sub {
	
	die unless(@lines_current_read);
	
	my $completeReadSequence_plus;
	my $completeReadSequence_minus;
	my $readID;
	my %alignments_perChr_perStrand;

	for(my $lineI = 0; $lineI <= $#lines_current_read; $lineI++)
	{
		my $line = $lines_current_read[$lineI];
		next unless($line);
		
		my @line_fields = split(/\t/, $line, -1);
		die unless($#line_fields == $#alignments_headerFields);
		my %line = (mesh @alignments_headerFields, @line_fields);
		
		if($lineI == 0)
		{
			$readID = $line{readID};
			print "Processing $readID\n";			
		}
		else
		{
			die unless($readID eq $line{readID});
		}
				
		if($line{completeReadSequence_plus})
		{
			die unless($line{completeReadSequence_minus});
			die if(defined $completeReadSequence_plus);
			$completeReadSequence_plus = $line{completeReadSequence_plus};
			$completeReadSequence_minus = $line{completeReadSequence_minus};
		}
		
		unless((defined $line{chromosome}) and (defined $line{strand}))
		{
			die Dumper("Lines with undefined fields", \@alignments_headerFields, \@line_fields, \%line);
		}
			
		push(@{$alignments_perChr_perStrand{$line{chromosome}}{$line{strand}}}, \%line);
	}
	
	if($processed_readIDs{$readID})
	{
		die "Read ID $readID has been seen already - are you using a sorted input file?";
	}
	$processed_readIDs{$readID}++;
	
	die Dumper("No complete sequences", \@lines_current_read) unless((defined $completeReadSequence_plus) and (defined $completeReadSequence_minus));
	die Dumper("No defined complete read sequence?", \@lines_current_read) unless($completeReadSequence_plus and $completeReadSequence_minus);
	
	my @finalScores;
	my @finalScores_origins;
	my @finalScores_chromosome;
	my @finalScores_strand;
	
	my %chainAddr_2_i;
	foreach my $chromosome (keys %alignments_perChr_perStrand)
	{
		die unless(exists $reference_href->{$chromosome});
		
		foreach my $strand (keys %{$alignments_perChr_perStrand{$chromosome}})
		{
			my @chains = enrichChains($alignments_perChr_perStrand{$chromosome}{$strand}, $completeReadSequence_plus, $completeReadSequence_minus, $reference_href, $chromosome);
			
			die unless(scalar(@chains));
			 
			foreach my $chain (@chains)
			{
				die unless((defined $chain->{n_matches}) and (defined $chain->{n_mismatches}) and (defined $chain->{n_gaps}));
				my $chain_score = $S_match * $chain->{n_matches} + $S_mismatch * $chain->{n_mismatches} + $S_gap * $chain->{n_gaps};
				$chain->{chainTraversalScore} = $chain_score;
			}
			
			@chains = sort {
				die unless(defined $a->{firstPos_reference});
				die unless(defined $b->{firstPos_reference});
				die unless(defined $a->{firstPos_read});
				die unless(defined $b->{firstPos_read});				
				
				if($a->{firstPos_reference} == $b->{firstPos_reference})
				{
					$a->{firstPos_read} <=> $b->{firstPos_read}
				}
				else
				{
					$a->{firstPos_reference} <=> $b->{firstPos_reference}				
				}
			} @chains;

			%chainAddr_2_i = ();
			if($strand eq '+')
			{
				#print "Collected chains:\n";
				for(my $chainI = 0; $chainI <= $#chains; $chainI++)
				{
					my $chain = $chains[$chainI];
					$chain->{enrich} = 0 unless(defined $chain->{enrich});
					#print "Chain ${chainI}; reference ", $chain->{firstPos_reference}, " - ", $chain->{lastPos_reference}, ", read ", $chain->{firstPos_read}, " - ", $chain->{lastPos_read}, "; score ", $chain->{chainTraversalScore}, "; enrich $chain->{enrich}\n";
					$chainAddr_2_i{$chain} = $chainI;
				}
				#print "";
				print "Sequence length: ", length($completeReadSequence_plus), ", reference length: ", length($reference_href->{$chromosome}), "\n\n";
				
			}
			for(my $chainI = 0; $chainI <= $#chains; $chainI++)
			{
				my $chain = $chains[$chainI];
				
				my @inputScores;
				my @inputScoreOrigins;
				my @inputScores_preGaps;
				
				if($endsFree_reference)
				{
					push(@inputScores, $S_gap * $chain->{firstPos_read});
				}
				else
				{
					push(@inputScores, $S_gap * ($chain->{firstPos_read} + $chain->{firstPos_reference}));
				}
				
				push(@inputScoreOrigins, 0);
				push(@inputScores_preGaps, 0);
				
				for(my $chainII = 0; $chainII < $#chains; $chainII++)
				{
					if($chainII < $chainI)
					{
						my $chain2 = $chains[$chainII];
						if(($chain2->{lastPos_reference} < $chain->{firstPos_reference}) and ($chain2->{lastPos_read} < $chain->{firstPos_read}))
						{
							my $delta_reference = $chain->{firstPos_reference} - $chain2->{lastPos_reference} - 1;
							my $delta_read = $chain->{firstPos_read} - $chain2->{lastPos_read} - 1;
							die unless(defined $chain2->{chainOutputScore});
							push(@inputScores, $chain2->{chainOutputScore} + $S_gap * ($delta_read + $delta_reference));
							push(@inputScores_preGaps, $chain2->{chainOutputScore});
							push(@inputScoreOrigins, $chain2);
						}
					}
					else
					{
						my $chain2 = $chains[$chainII];
						die if(($chain2->{lastPos_reference} < $chain->{firstPos_reference}) and ($chain2->{lastPos_read} < $chain->{firstPos_read}));			
					}
				}
				die unless($#inputScores == $#inputScoreOrigins);
				
				if(0 and $chainI == 40)
				{
					# print "Chain $chainI incoming scores:\n";	
					for(my $fromChainI = 0; $fromChainI <= $#inputScores; $fromChainI++)
					{
						my $origin = $inputScoreOrigins[$fromChainI];
						$origin = ($origin) ? $chainAddr_2_i{$origin} : "ORIGIN";
						# print " - $origin $inputScores[$fromChainI] -- $inputScores_preGaps[$fromChainI] \n";
					}
				}
				my $max_input = which_max(\@inputScores);
				my $bestInputScore = $inputScores[$max_input];
				my $bestInputScore_origin = $inputScoreOrigins[$max_input];
				
				$chain->{chainOutputScore} = $bestInputScore + $chain->{chainTraversalScore};
				$chain->{chainOutputScoreOrigin} = $bestInputScore_origin;
			}
			
			my @finalOutputScores;
			my @finalOutputScoreOrigins;
			
			if($endsFree_reference)
			{
				push(@finalOutputScores, $S_gap * (length($completeReadSequence_plus)));
			}
			else
			{
				push(@finalOutputScores, $S_gap * (length($completeReadSequence_plus) + length($reference_href->{$chromosome})));
			}
			push(@finalOutputScoreOrigins, 0);
			
			for(my $chainI = 0; $chainI <= $#chains; $chainI++)
			{
				my $chain = $chains[$chainI];
				die unless(defined $chain->{chainOutputScore});
				die unless(defined $chain->{chainOutputScoreOrigin});
				my $final_delta_read = length($completeReadSequence_plus) - $chain->{lastPos_read} - 1;
				my $final_delta_ref = length($reference_href->{$chromosome}) - $chain->{lastPos_reference} - 1;
				if($endsFree_reference)
				{
					push(@finalOutputScores, $chain->{chainOutputScore} + $S_gap * $final_delta_read);
				}
				else
				{
					push(@finalOutputScores, $chain->{chainOutputScore} + $S_gap * ($final_delta_read + $final_delta_ref));
				}
				push(@finalOutputScoreOrigins, $chain);
			}			
			
			my $i_of_max_finalScore = which_max(\@finalOutputScores);
			my $maxFinalScore = $finalOutputScores[$i_of_max_finalScore];
			my $maxFinalScore_origin = $finalOutputScoreOrigins[$i_of_max_finalScore];

			push(@finalScores, $maxFinalScore);
			push(@finalScores_origins, $maxFinalScore_origin);
			push(@finalScores_chromosome, $chromosome);
			push(@finalScores_strand, $strand);
		}
	}
	
	die unless(scalar(@finalScores));
	my $i_of_max_finalfinalScore = which_max(\@finalScores);

	my $strand = $finalScores_strand[$i_of_max_finalfinalScore];
	my $chromosome = $finalScores_chromosome[$i_of_max_finalfinalScore];
	
	die unless(($strand eq '+') or ($strand eq '-'));
	my $useReadSequence = ($strand eq '+') ? $completeReadSequence_plus : $completeReadSequence_minus;
	
	my $bt_reference = '';
	my $bt_contig = '';
	
	my $next_bt_position = $finalScores_origins[$i_of_max_finalfinalScore];
	my $last_emitted_read_position = length($useReadSequence);
	my $last_emitted_reference_position = -1;
	my $min_emitted_reference_position;
	my $max_emitted_reference_position;
	
	my $leftOut_delta_reference = 0;
	my $n_thisBacktrace_chains = 0;
	while($last_emitted_read_position > 0)
	{
		die unless(defined $next_bt_position);
		if($next_bt_position == 0)
		{
			my $missing_contig_sequence = substr($useReadSequence, 0, $last_emitted_read_position);
			my $equivalent_reference_sequence_gap = ('-' x length($missing_contig_sequence));
			
			# print "Left end, add delta ", $last_emitted_reference_position, "\n";
			$leftOut_delta_reference += ($last_emitted_reference_position);	
			
			$bt_contig .= reverse($missing_contig_sequence);
			$bt_reference .= $equivalent_reference_sequence_gap;
			
			$last_emitted_read_position = 0;
			$last_emitted_reference_position = -1;
		}
		else
		{
			my $chain_for_consumption = $next_bt_position;
			$n_thisBacktrace_chains++;
			
			if(exists $chainAddr_2_i{$chain_for_consumption})
			{
				# print "Bt chain $chainAddr_2_i{$chain_for_consumption} \n";
			}
			my $delta_read = $last_emitted_read_position - $chain_for_consumption->{lastPos_read} - 1;
			if($n_thisBacktrace_chains == 1)
			{
				# print "Right end, add delta ", (length($reference_href->{$chromosome}) - $chain_for_consumption->{lastPos_reference} - 1), "\n";
				
				$leftOut_delta_reference += (length($reference_href->{$chromosome}) - $chain_for_consumption->{lastPos_reference} - 1);
			}	
			die unless($delta_read >= 0);
			
			my $betweenChains_read = '';
			my $betweenChains_reference = '';
			
			my $jumpedOverRead_read = substr($useReadSequence, $chain_for_consumption->{lastPos_read} + 1, $last_emitted_read_position  - 1 - ($chain_for_consumption->{lastPos_read} + 1) + 1);
			die unless(length($jumpedOverRead_read) == $delta_read);			
			my $jumpedOverRead_reference = ('-' x $delta_read);
			die unless(length($jumpedOverRead_reference) == $delta_read);			
			
			$betweenChains_read .= $jumpedOverRead_read;
			$betweenChains_reference .= $jumpedOverRead_reference;
			
			if($last_emitted_reference_position != -1)
			{
				my $delta_reference = $last_emitted_reference_position - $chain_for_consumption->{lastPos_reference} - 1;
				
				my $jumpedOverReference_reference =  substr($reference_href->{$chromosome}, $chain_for_consumption->{lastPos_reference}+1, $last_emitted_reference_position - $chain_for_consumption->{lastPos_reference} - 1);
				die unless(length($jumpedOverReference_reference) == $delta_reference);				
				my $jumpedOverReference_read = ('-' x $delta_reference);
				die unless(length($jumpedOverReference_read) == $delta_reference);				
						
				$betweenChains_read .= $jumpedOverReference_read;
				$betweenChains_reference .= $jumpedOverReference_reference;				
			}

			die unless(length($betweenChains_read) == length($betweenChains_reference));
			
			$bt_contig .= reverse($chain_for_consumption->{alignment_read} . $betweenChains_read);
			$bt_reference .= reverse($chain_for_consumption->{alignment_reference} . $betweenChains_reference);
			
			die Dumper("Error I", $chain_for_consumption->{firstPos_read}, $last_emitted_read_position) unless($chain_for_consumption->{firstPos_read} <= $last_emitted_read_position);
			$last_emitted_read_position = $chain_for_consumption->{firstPos_read};
			
			if($last_emitted_reference_position != -1)
			{
				die Dumper("Error 2", $chain_for_consumption->{firstPos_reference}, $last_emitted_reference_position)  unless($chain_for_consumption->{firstPos_reference} <= $last_emitted_reference_position);
			}
			$last_emitted_reference_position = $chain_for_consumption->{firstPos_reference};
			die Dumper("Error 3", $last_emitted_read_position)  unless($last_emitted_reference_position >= 0);
			
			$next_bt_position = $chain_for_consumption->{chainOutputScoreOrigin};
			
			if(not defined $max_emitted_reference_position)
			{
				$max_emitted_reference_position = $chain_for_consumption->{lastPos_reference};
			}
			
			if((not defined $min_emitted_reference_position) or ($chain_for_consumption->{firstPos_reference} < $min_emitted_reference_position))
			{
				$min_emitted_reference_position = $chain_for_consumption->{firstPos_reference};
			}
		}
	}
	
	if($last_emitted_reference_position != -1)
	{
		#print "Final step with last_emitted_read_position = $last_emitted_read_position and last_emitted_reference_position = $last_emitted_reference_position :\n";
		#print $last_emitted_reference_position, "\n";
		$leftOut_delta_reference += ($last_emitted_reference_position);	
	}
	
	$bt_contig = reverse($bt_contig);
	$bt_reference = reverse($bt_reference);
	
	my $bt_contig_noGaps = $bt_contig;
	$bt_contig_noGaps =~ s/\-//g;
	die Dumper("Error 4", $useReadSequence, $bt_contig_noGaps, length($bt_contig_noGaps), length($useReadSequence)) unless($bt_contig_noGaps eq $useReadSequence);
	
	my $bt_reference_noGaps = $bt_reference;
	$bt_reference_noGaps =~ s/\-//g;
	my $reference_extract;
	if($bt_reference_noGaps)
	{
		die unless(defined $min_emitted_reference_position);
		die unless(defined $max_emitted_reference_position);
		die unless($min_emitted_reference_position <= $max_emitted_reference_position);
		$reference_extract = substr($reference_href->{$chromosome}, $min_emitted_reference_position, $max_emitted_reference_position - $min_emitted_reference_position + 1);
		my $reference_extract_expected_length = $max_emitted_reference_position - $min_emitted_reference_position + 1;
		 die unless(length($reference_extract) == $reference_extract_expected_length);
		if(index($reference_extract, "N") == -1)
		{
			if(($bt_reference_noGaps ne $reference_extract) and (length($bt_reference_noGaps) == length($reference_extract)))
			{
				for(my $i = 0; $i < length($bt_reference_noGaps); $i++)
				{
					my $c1 = substr($bt_reference_noGaps, $i, 1);
					my $c2 = substr($reference_extract, $i, 1);
					if($c1 ne $c2)
					{
						print "Mismatch position $i - $c1 vs $c2\n";
					}
				}
			}	
			die Dumper("Reference sequence mismatch", length($bt_reference_noGaps), length($reference_extract)) unless($bt_reference_noGaps eq $reference_extract);
			#die Dumper("Reference sequence mismatch",$bt_reference_noGaps, $reference_extract, length($bt_reference_noGaps), length($reference_extract)) unless($bt_reference_noGaps eq $reference_extract);
		}
	}
	
	
	die unless(length($bt_contig) == length($bt_reference));
	my $score_reconstructed = 0;
	my $n_mismatches = 0;
	for(my $i = 0; $i < length($bt_contig); $i++)
	{
		my $c1 = substr($bt_contig, $i, 1);
		my $c2 = substr($bt_reference, $i, 1);
		
		if(($c1 eq '-') or ($c2 eq '-'))
		{
			$score_reconstructed += $S_gap;
			$n_mismatches++;
		}
		else
		{
			if($c1 eq $c2)
			{
				$score_reconstructed += $S_match;
			}
			else
			{
				$score_reconstructed += $S_mismatch;
				$n_mismatches++;
			}
		}
	
	}
	
	my $score_before_endsFreeChange = $score_reconstructed;
	
	if(not $endsFree_reference)
	{
		$score_reconstructed += ($leftOut_delta_reference * $S_gap);
	}	
	
	unless($score_reconstructed == $finalScores[$i_of_max_finalfinalScore])
	{
		die Dumper("Score mismatch", length($useReadSequence), $bt_reference, $bt_contig, $score_reconstructed, $finalScores[$i_of_max_finalfinalScore]);
	}
	
	if($bt_reference_noGaps)
	{
		my $remove_columns_front = 0;
		while(substr($bt_reference, $remove_columns_front, 1) eq '-')
		{
			$remove_columns_front++;
		}
		
		my $remove_columns_back = 0;
		while(substr($bt_reference, length($bt_reference) - $remove_columns_back - 1, 1) eq '-')
		{
			$remove_columns_back++;
		}
		
		my $alignment_reference_forSAM = substr($bt_reference, $remove_columns_front, length($bt_reference) - $remove_columns_front - $remove_columns_back);
		my $alignment_contig_forSAM = substr($bt_contig, $remove_columns_front, length($bt_contig) - $remove_columns_front - $remove_columns_back);
		die unless(length($alignment_reference_forSAM) == length($alignment_contig_forSAM));
		
		my $alignment_reference_forSAM_noGaps = $alignment_reference_forSAM;
		$alignment_reference_forSAM_noGaps =~ s/-//g;
		if(index($reference_extract, "N") == -1)
		{
			die unless($alignment_reference_forSAM_noGaps eq $reference_extract);
		}

		die if($alignment_reference_forSAM =~ /^-/);
		die if($alignment_reference_forSAM =~ /-$/);
		
		my $alignment_contig_forSAM_noGaps = $alignment_contig_forSAM;
		$alignment_contig_forSAM_noGaps =~ s/-//g;		
		my $tlen = $max_emitted_reference_position - $max_emitted_reference_position + 1;
		
		my @CIGAR;
		my $addCharacter = sub {
			my $character = shift;
			die unless(length($character) == 1);
			if((scalar(@CIGAR) == 0) or ($CIGAR[$#CIGAR][1] ne $character))
			{
				push(@CIGAR, [0, $character]);
			}
			$CIGAR[$#CIGAR][0]++;
		};
		
		for(my $i = 0; $i < length($alignment_reference_forSAM); $i++)
		{
			my $c_ref = substr($alignment_reference_forSAM, $i, 1);
			my $c_read = substr($alignment_contig_forSAM, $i, 1);
			
			die if(($c_ref eq '-') and ($c_read eq '-'));
			if($c_ref eq '-')
			{
				$addCharacter->('I');
			}
			elsif($c_read eq '-')
			{
				$addCharacter->('D');			
			}
			else
			{
				$addCharacter->('M');							
			}
		}
		
		my $CIGAR = join('', map {die unless(scalar(@{$_}) == 2); $_->[0] . $_->[1]} @CIGAR);
		die if($CIGAR =~ /H/i);
		
		my @fields_for_output = (
			$readID,
			(2 | 64),
			$chromosome,
			$min_emitted_reference_position+1,
			255,
			$CIGAR,
			'*',
			'0',
			$tlen,
			$alignment_contig_forSAM_noGaps,
			'*'
		);
		
		print SAMOUTPUT join("\t", @fields_for_output), "\n";
		
		$n_chains_sum += $n_thisBacktrace_chains;		
		$n_output_alignments++;
		
		if($remove_columns_front and $remove_columns_back)
		{
			$n_alignments_leftAndRightGapsRemoved++;
		}
		elsif($remove_columns_front)
		{ 
			$n_alignments_leftGapsRemoved++;
		}
		elsif($remove_columns_back)
		{
			$n_alignments_rightGapsRemoved++;
		}
		if($remove_columns_front or $remove_columns_back)
		{
			print TRUNCATED $readID, "\n";
		}
		
		print LENGTHS join("\t", $readID, length($alignment_contig_forSAM_noGaps)), "\n";
	}
};


while(<INPUT>)
{
	if(($. % 10000) == 0)
	{
		print "Read line $.\n";
	}
	
	my $line = $_;
	chomp($line);
	next unless($line);
	die unless($line =~ /^(.+?)\t/);
	my $readID = $1;
	if($readID ne $currentReadID)
	{
		if(@lines_current_read)
		{
			# print "Call $readID line $.\n";
			$processReadLines->();
		}
		@lines_current_read = ();
		$currentReadID = $readID;
	}
	
	push(@lines_current_read, $line);
}
close(INPUT);
if(@lines_current_read)
{
	# print "Call end\n";
	$processReadLines->();
}
 
close(SAMOUTPUT);

print "\nDone.\n\nStatistics:\n";
print "\tAlignments: ", $n_output_alignments, "\n";
print "\tAvg chains per alignment: ", ($n_chains_sum / $n_output_alignments), "\n";
print "\tAlignments with both ends trimmed: ", $n_alignments_leftAndRightGapsRemoved, "\n";
print "\tAlignments with left end trimmed: ", $n_alignments_leftGapsRemoved, "\n";
print "\tAlignments with right end trimmed: ", $n_alignments_rightGapsRemoved, "\n";

print "\n\nDone. Produced (unsorted) SAM file $outputFile_sam_unsorted - this file will be deleted.\n\n";

# my $outputFile_sam_unsorted_filtered = $outputFile_sam_unsorted . ".filtered";
# my $cmd_filter = qq(perl $CIGARscript_path --input ${outputFile_sam} --output ${outputFile_sam_filtered});
# print "Filtering SAM with command: $cmd_filter\n";
# die "Filtering failed" unless(system($cmd_filter) == 0);

# my $bam_unsorted = $outputFile_sam_unsorted_filtered . '.bam';
# my $cmd_bam_conversion = qq($samtools_path view -S -b -o${bam_unsorted} $outputFile_sam_unsorted_filtered);

# print "Converting to BAM with command:\n\t$cmd_bam_conversion\n\n";
# die "BAM conversion failed" unless(system($cmd_bam_conversion) == 0);

my $cmd_sam_sort = qq($samtools_path sort -o${outputFile} -O sam $outputFile_sam_unsorted);
die "SAM sorting failed" unless(system($cmd_sam_sort) == 0);

# print "Sorting and indexing BAM with command:\n\t$cmd_bam_sort\n\n";
# die "BAM sorting/indexing failed" unless(system($cmd_bam_sort) == 0);

print "Produced (sorted) SAM file $outputFile\n\n";
unlink($outputFile_sam_unsorted);

sub which_max
{
	my $input_aref = shift;
	die unless(@$input_aref);
	my $max = max(@$input_aref);
	for(my $i = 0; $i <= $#{$input_aref}; $i++)
	{
		return $i if($input_aref->[$i] == $max);
	}
	die;
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

sub enrichChains
{
	my $chainsIn_aref = shift;
	my $completeReadSequence_plus = shift;
	my $completeReadSequence_minus = shift;
	my $reference_href = shift;
	my $chromosome = shift;
	die if($completeReadSequence_plus eq $completeReadSequence_minus);
	die unless(length($completeReadSequence_plus) == length($completeReadSequence_minus));
	die unless(defined $completeReadSequence_minus);
	
	my @chainsIn = @$chainsIn_aref;
	
	my %beginCoordinates_reference;
	my %beginCoordinates_read;
	foreach my $chain (@chainsIn)
	{
		$beginCoordinates_reference{$chain->{firstPos_reference}}++;
		$beginCoordinates_read{$chain->{firstPos_read}}++;
	}
	
	
	my @chainsOut;	
	for(my $chainI = 0; $chainI <= $#chainsIn; $chainI++)
	{
		my %createdEndPointForReferenceCoordinate;
		my %createdEndPointForReadCoordinate;
		
		my $chain = $chainsIn[$chainI];
		
		my $useReadSequence = ($chain->{strand} eq '+') ? $completeReadSequence_plus : $completeReadSequence_minus;				
	
	
		{
			my $supposedReferenceSequence = substr($reference_href->{$chromosome}, $chain->{firstPos_reference}, $chain->{lastPos_reference} - $chain->{firstPos_reference} + 1);		
			if(index($supposedReferenceSequence, 'N') == -1)
			{
				my $chainSequence_reference_noGaps = $chain->{alignment_reference};
				$chainSequence_reference_noGaps =~ s/-//g;
				die Dumper("Error 5", length($supposedReferenceSequence), length($chainSequence_reference_noGaps), $supposedReferenceSequence, $chainSequence_reference_noGaps) unless($supposedReferenceSequence eq $chainSequence_reference_noGaps);		
			}
		}
		push(@chainsOut, $chain);
		
		
		my $n_matches = 0;
		my $n_mismatches = 0;
		my $n_gaps = 0;
		my $alignment_reference = '';
		my $alignment_read = '';
		my $lastPos_reference = undef;
		my $lastPos_read = undef;
		for(my $alignmentPos = 0; $alignmentPos < length($chain->{alignment_reference}); $alignmentPos++)
		{
			my $character_reference = substr($chain->{alignment_reference}, $alignmentPos, 1);
			my $character_read = substr($chain->{alignment_read}, $alignmentPos, 1);
			
			$alignment_reference .= $character_reference;
			$alignment_read .= $character_read;
			
			if($character_reference eq $character_read)
			{
				$n_matches++;
			}
			else
			{
				if(($character_reference eq '-') or ($character_read eq '-'))
				{
					$n_gaps++;
				}
				else
				{
					$n_mismatches++;
				}
			}
			
			if($character_reference ne '-')
			{
				if(defined $lastPos_reference)
				{
					$lastPos_reference++;
				}
				else
				{
					$lastPos_reference = $chain->{firstPos_reference};
				}
				
			}
			if($character_read ne '-')
			{
				if(defined $lastPos_read)
				{
					$lastPos_read++;
				}
				else
				{
					$lastPos_read = $chain->{firstPos_read};
				}			
			}
			
			my $includeCurrentPositionAsEndpoint = 0;
			if((defined $lastPos_reference) and (defined $lastPos_read))
			{
				if(exists $beginCoordinates_reference{$lastPos_reference + 1})
				{
					if(not $createdEndPointForReferenceCoordinate{$lastPos_reference + 1})
					{
						$includeCurrentPositionAsEndpoint = 1;
					}
				}
				
				if(exists $beginCoordinates_read{$lastPos_read + 1})
				{
					if(not $createdEndPointForReadCoordinate{$lastPos_read + 1})
					{
						$includeCurrentPositionAsEndpoint = 1;
					}
				}
			}		

			if($includeCurrentPositionAsEndpoint)
			{
				$createdEndPointForReferenceCoordinate{$lastPos_reference + 1} = 1;
				$createdEndPointForReadCoordinate{$lastPos_read + 1} = 1;
				
				my $supposedReferenceSequence = substr($reference_href->{$chromosome}, $chain->{firstPos_reference}, $lastPos_reference - $chain->{firstPos_reference} + 1);
				if($supposedReferenceSequence =~ /^[ACGT]+$/i)
				{								
					my $chainSequence_reference_noGaps = $alignment_reference;
					$chainSequence_reference_noGaps =~ s/-//g;
					die Dumper("Sequence mismatch", $supposedReferenceSequence, $chainSequence_reference_noGaps, $chain, $chromosome) unless($supposedReferenceSequence eq $chainSequence_reference_noGaps);
					# print "Checked reference sequence.\n";
				}
				
				my $supposedReadSequence = substr($useReadSequence, $chain->{firstPos_read}, $lastPos_read - $chain->{firstPos_read} + 1);
				my $chainSequence_read_noGaps = $alignment_read;
				$chainSequence_read_noGaps =~ s/-//g;
				
				die Dumper("Error 6", ['Strand', $chain->{strand}, (($completeReadSequence_minus eq $useReadSequence) ? '-' : '+')], $chain->{chromosome}, $character_reference, $character_read, $alignmentPos, $chain->{firstPos_read}, $chain->{lastPos_read}, $lastPos_read, length($supposedReadSequence), length($chainSequence_read_noGaps), ['Alignments', substr($alignment_reference, 0, 20), substr($alignment_read, 0, 20)], [substr($supposedReadSequence, 0, 20), substr($chainSequence_read_noGaps, 0, 20)], [substr($supposedReadSequence, length($supposedReadSequence) - 20, 20), substr($chainSequence_read_noGaps, length($chainSequence_read_noGaps) - 20, 20)]) unless($supposedReadSequence eq $chainSequence_read_noGaps);				

				push(@chainsOut, {
					readID => $chain->{readID},
					chromosome => $chain->{chromosome},
					firstPos_reference => $chain->{firstPos_reference},
					lastPos_reference => $lastPos_reference,
					firstPos_read => $chain->{firstPos_read},
					lastPos_read => $lastPos_read, 
					strand => $chain->{strand},
					n_matches => $n_matches,
					n_mismatches => $n_mismatches,
					n_gaps => $n_gaps,
					alignment_reference => $alignment_reference,
					alignment_read => $alignment_read,
					enrich => 1,
				});
			}			
		}
	}
	
	print "\tChain enrichment: from ", scalar(@chainsIn), " to ", scalar(@chainsOut), " chains.\n";
	return @chainsOut;
}

