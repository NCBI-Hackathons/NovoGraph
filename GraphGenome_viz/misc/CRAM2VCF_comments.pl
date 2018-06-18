#!/usr/bin/env perl

## Dilthey comments, 16 Jan 2018

## ## Data preparation:

## So we start at the position of Q4 - i.e. the line $query_noGaps =~ s/[\-_\*]//g; .
##
## The subsequent blocks have the job of cutting the contig alignments to areas that begin
## and end with matches (for we don't know how to treat dangling sequences that aren't anchored
## to the reference), and to make sure nothing weird has happened during these transformations.
## This block ends around line 318 (in your version).
##
## Then we determine the last reference position of the remaining (i.e. after removal of non-match
## positions at beginning/end) contig-to-reference alignment - this is done in the block starting with
##
## '''
## for(my $i = 0; $i < length($ref); $i++)
## '''
##
## .. wherein we just increment our 'reference position' variable whenever the alignment doesn't
##    contain a non-gap reference character
##
## (it's useful to think about these alignments the way they're displayed by e.g. BLAST, like
##
## '''
## ACGTAC-TT
## ACT-AACTT
## '''
##
## , and we step through these in a column-wise fashion)
##
## Finally, we store the alignments with some additional data like last position, indexed by
## starting position, like so:
##
## '''
## push(@{$alignments_starting_at{$alignment_start_pos}}, $alignment_info_aref); ## add elements to array
## '''
##
## ... and we also store these data in a file that's then processed by the BAM2VCF C++ component:
##
## '''
## print D join("\t", $ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos), "\n"; ### write out in VCF format
## '''
##
## (so this is just another step in preparation for the VCF generation, not the VCF generation
## itself - which is later called with
## my $cmd = qq($bin_BAM2VCF --input $fn_for_BAM2VCF --referenceSequenceID $referenceSequenceID);.
##
##
## ## Validating the gap structure of the alignments:
##
## This is the whole business from line 376 to line 433.
##
## This is to make sure that we're dealing with a valid multiple sequence alignment of the input
## sequences, after the modifications to the alignments (removal of non-match begin/end positions)
## that are carried out in this script.
##
## What I mean by this is that the data structure we're encoding is something like this:
##
## '''
## ref     ACGTAC-TT
## contig1 ACT-AACTT
## contig2 ATTTAACTT
## contig3 ACT-AA-TT
## '''
##
## ... but BAM/CRAM only store pairwise alignments, i.e.
##
## '''
## ref     ACGTAC-TT
## contig1 ACT-AACTT
##
## ref     ACGTAC-TT
## contig2 ATTTAACTT
##
## ref     ACGTAC-TT
## contig3 ACT-AA-TT
## '''
##
##  (note the double-gap column in the last alignment, which is unusual but can be encoded
##  via the CIGAR string in BAM/CRAM).
##
## ... so what we're trying to do is make sure that the individual pairwise alignments can be transformed
## into a valid MSA structure, which is the case iff
## - the reference sequences of the alignments are identical and
## - the reference alignments encode identical gap structures along the reference coordinate
##
## What complicates things a bit is that not all pairwise alignments start at the same position per se,
## so we have to do this gap-counting business to make sure everything lines up.
##
## ## Intermediate output
##
## So one question is whether the output from CRAM2VCF is complete (i.e. we're not losing contigs).
## I think I did check this some time ago, but it should be obvious from the output if you re-run the pipeline.
##
## Then the next step is to look at the actual VCF generation, which happens in here:
##
## https://github.com/NCBI-Hackathons/Graph_Genomes_CSHL/blob/master/BAM2VCF/BAM2VCF.cpp
##
## ## VCF generation
##
## So what happens in the C++ script?
## We load the output file from CRAM2VCF, which basically contains alignments by alignment starting position
## We filter these so that they don't contain gaps that are 'too large'
## - if(max_running_gap_length <= max_gap_length) --
## max_gap_length might be a crucial parameter for fixing our issue that we should play around with.
##
## We call the function produceVCF to produce the VCF.
##
## Function produceVCF:
## We do some preliminary checks (alignments are gap-structure compatible, basically identical to the
## last part of CRAM2VCF) and we calculate and output positional coverage. This would probably be very
## useful for tracking down missing variants, i.e. by checking whether the contigs we expect to be at
## a position really go into the VCF generation process.
##
## '''
## if((pI >= 15000000) && (pI <= 17000000)) std::cout << "\t" << "Window starting at pI = " << pI << " => avg. coverage " << avg_coverage << "\n"
## '''
##
##
## We then start building the graph reference position by reference position.
##
## 'Graph building' means that we compute possible haplotypes spanned by the input alignments
##  (allowing for recombination) and encode them as VCF.
##
## The key data structure here is
##
## '''
## using openHaplotype = std::tuple<std::string, const startingHaplotype*, int>;
## std::vector<openHaplotype> open_haplotypes;
## '''
##
## ... where open_haplotypes always keeps track of the graph sequences we're currently constructing.
## The first element (std::string) keeps the sequence being constructed; startingHaplotype* is a pointer
## to the alignment we're copying sequence from into the haplotype; the int element is a pointer that tells
## us how far into the alignment we've progressed.
##
## If the startingHaplotype element is 0, this means that we're copying from the reference sequence.
##
## I'll try to add more comments to the C++ source file directly!
## 



# Example command:
# To check correctness of INPUT for mafft:
# 	./BAM2VCF.pl --BAM /data/projects/phillippy/projects/hackathon/Graph_Genomes_CSHL/scripts/uber_sorted.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output uber_vcf.vcf
# ./BAM2VCF.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.fa --output uber_vcf.vcf


### standard libraries/modules, including Perl's HTSLIB

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;


$| = 1;

### perl voodo
### $| is the built-in variable for autoflush, cf. "perl pre-defined variables", perlvar
### $| = 1; forces a flush after every write or print, so the output apperas as soon as it's generated rather than being buffered
### $OUTPUT_AUTOFLUSH; If set to nonzero, forces an fflush(3) after every write or print on the currently selected output channel.

### VARIABLES: CRAM, hg38 fasta, output VCF (?), subdirectory of BAM2VCF, contigLengths (?) (what's this file/file format?)

## RE: what's this file/file format?
## contig ID and length are tab-separated -- no header
## e.g. 
## > CHM13.gi|953902711|gb|LDOC03012091.1|   16569
## > CHM13.gi|953902713|gb|LDOC03012090.1|   48337
## > CHM13.gi|953902714|gb|LDOC03012089.1|   2193864


### RE: $output ---- we use `my $fn_for_BAM2VCF = $output . '.part_'. $referenceSequenceID;` below


my $CRAM;
my $referenceFasta;
my $output;
my $bin_BAM2VCF = '../BAM2VCF/BAM2VCF';
my $contigLengths;                          ### optional parameter

unless(-e $bin_BAM2VCF)   ### if not a subdirectory, kill; see checks below
{
	die "BAM2VCF binary $bin_BAM2VCF not present - run 'make all' in the directory.";
}

#printHaplotypesAroundPosition(2, {
	# 0 => [
		# ['AAAC', 'TTTT', 'r0', 0, 3],
	# ],
	# 2 => [
		# ['ACGT', 'TTTT', 'r1', 2, 5],
		# ['AC-GT', 'TTATT', 'r2', 2, 5],	
		# ['ACGT', 'AA_A', 'r3', 2, 5],		
	# ],
	#2 => [
		# ['-CGT', 'TTTT', 'r4', 3, 5],	
	#	['ACG-', 'TTTT', 'r5', 2, 4],	
	#],
	#0 => [
	#	 ['AAAC', 'TTT-', 'r7', 0, 3],
	#],	
#});
 
### kill program if CRAM, referenceFasta, output, contigLengths doesn't exist

GetOptions (
	'CRAM:s' => \$CRAM, 
	'referenceFasta:s' => \$referenceFasta, 
	'output:s' => \$output,
	'contigLengths:s' => \$contigLengths, 	   ### optional parameter
);
	
die "Please specify --CRAM" unless($CRAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --output" unless($output);

die "--CRAM $CRAM not existing" unless(-e $CRAM);    ### '-e', checks whether exists at certain path
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);



my %expectedLengths;    ### hash defined, 'expectedLengths'
if($contigLengths)      ### read in tab-delimited file
{
	open(L, '<', $contigLengths) or die "Cannot open --contigLengths $contigLengths";
	while(<L>)
	{
		my $l = $_; 
		chomp($l);
		my @f = split(/\t/, $l);
		$expectedLengths{$f[0]} = $f[1];    ### get value of expectedLengths
	}
	close(L);
}

### -> operator like C, means 'deference'
### in this case, means 'Class -> method()'

my $sam = Bio::DB::HTS->new(-fasta => $referenceFasta, -bam => $CRAM);   ### standard Bio::DB::HTS constructor; 'prepare' FASTA, SAM for loading

my @sequence_ids = $sam->seq_ids();    ### define array 'sequence_ids', sequence IDs from the BAM

my $testing = 0;   ## if == 1, trips the $testing conditionals below; includes 'try # 2' and unit tests

my $reference_href;   ## variable for 'reference FASTA hg38'
if(not $testing)    ## not sure why this would be in 'testing' ....
{
	print "Reading $referenceFasta\n";
	$reference_href = readFASTA($referenceFasta);       ### load the reference FASTA; cf readFasta() function defined below
	print "\t...done.\n";
}

open(OUT, ">", $output) or die "Cannot open $output";   
### snazzy hack for printing out a VCF in the correct format
print OUT qq(##fileformat=VCFv4.2                       
##fileDate=20161026
##source=BAM2VCF.pl
##reference=file://$referenceFasta), "\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "\n";
#foreach my $referenceSequenceID (@sequence_ids)
my @referenceSequenceIDs = @sequence_ids;     ### sequence IDs for the FASTA
if($CRAM =~ /TRY2/)                           ### several tries for this script...from comment, looks like this should be removed
{
	# @referenceSequenceIDs = ("chr1");
}
my %alignments_starting_at_test;
if($testing)                                  ### quick unit test à la perl
{
	# ?CGT-A--CGT ref
	# -CTT-A--CGT r0
	#      AAACGT r1
	#   GTAT      r2
	#
		 
	@referenceSequenceIDs = qw/test1/;
	$reference_href = {
		test1 => "ACGTACGT",
	};
	
	%alignments_starting_at_test = (
	 1 => [
		 ['CGT-A--CGT', 'CTT-A--CGT', 'r0', 1, 7],
	 ],
	 4 => [
		 ['A--CGT', 'AAACGT', 'r1', 4, 7],
	 ],	 
	 2 => [
		 ['GT-A', 'GTAT', 'r2', 2, 4],
	 ],		 
	);
}


mkdir('forVCF');           ## mkdir 'forVCF'
die unless(-e 'forVCF');
# @referenceSequenceIDs = qw/chr21/;
foreach my $referenceSequenceID (@referenceSequenceIDs)
{
	print "Processing $referenceSequenceID .. \n";
	die "Sequence ID $referenceSequenceID not defined in $referenceFasta" unless(exists $reference_href->{$referenceSequenceID});  ### if come across a strange sequence ID, kill it
	

	## Dilthey comments, 14Jan18:
	## '$#gap_structure = ($l_ref_sequence - 1);' sets the size of the array 
	## i.e. after this command, the array will have ($l_ref_sequence - 1)+1 elements 
	## this is to speed up things a bit, because it preallocates the array. 
	## After having executed this command, the new array > elements will be undef. 
	## I think the -1 was another way to express 'not defined', but using the inbuilt undef seemed more practical to me in the end.

	my $l_ref_sequence = length($reference_href->{$referenceSequenceID});   ### length of sequence ID
	my @gap_structure;
	$#gap_structure = ($l_ref_sequence - 1); ## preallocates size of the array, which has '($l_ref_sequence - 1)+1' elements
	print "Set...";


	#for(my $i = 0; $i < $l_ref_sequence; $i++)
	#{
	#	#$gap_structure[$i] = -1;
	#}
	print " .. done.\n";
	die unless(scalar(@gap_structure) == $l_ref_sequence);  ## length of array 'gap_structure' should equal than of 'l_ref_sequence'; if not, kill the program
	
	my $fn_for_BAM2VCF = $output . '.part_'. $referenceSequenceID;
	open(D, '>', $fn_for_BAM2VCF) or die "Cannot open $fn_for_BAM2VCF";
	print D $reference_href->{$referenceSequenceID}, "\n";
	my $n_alignments = 0;          ### n_alignments=0; used for iteration for below; just like C/C++ iteration
	my %alignments_starting_at;
	if(not $testing)
	{
		my $alignment_iterator = $sam->features(-seq_id => $referenceSequenceID, -iterator => 1);	### iterate through all reference seq IDs
		while(my $alignment = $alignment_iterator->next_seq)
		{
			$n_alignments++;		## iteration n_alignments

			my $alignment_start_pos = $alignment->start - 1; 
			my ($ref,$matches,$query) = $alignment->padded_alignment;     
			unless(length($ref) == length($query))                        ### were there any warnings on the first run?
			{
				warn Dumper("Alignment length mismatch", length($ref), length($query), $alignment->query->name);
				next;
			}	

			## point of `if($contigLengths){}`
			## this is to make sure we're using correct contig sequences.
			
			if($contigLengths)
			## perl search-and-replace regex, '$string =~ s/regex/replacement/g;'
			{                                     
				my $query_noGaps = $query;
				$query_noGaps =~ s/[\-_\*]//g;   ## what this does: removes all -, _, * characters (i.e. substitute with zero-length string). 
				## Then we check whether the length of this string agrees with what we expect from the contig lengths file - a double-check against data corruption or loss.
				my $expectedLength = $expectedLengts{$alignment->query->name};
				die "No length for " . $alignment->query->name unless(defined $expectedLength);
				unless($expectedLength == length($query_noGaps))
				{
					die Dumper("Sequence length mismatch", $query_noGaps, length($query_noGaps), $expectedLength, length($alignment->query->dna), $alignment->query->name);
				}
			}

			
			
			my $ref_preAll = $ref;         ### cf. `my ($ref,$matches,$query) = $alignment->padded_alignment;` above  
			my $query_preAll = $query;
			
			my $firstMatch;
			my $lastMatch;
			### PSEUDOCODE
			###
			### for i in length(ref):
			###    if ( (ref[i:i+1] != '_') & (ref[i, i+1] != '*') & (query[i:i+1] != '_') & (query[i, i+1] != '*')):
			###        firstMatch = i   ### unless(defined $firstMatch);
			###        lastMatch = i
			###
			###    ### I guess the point of this is for the check here: 'die if($isDeletion and $isInsertion);''
			###    isDeletion = ( (ref[i:i+1] != '_') & (ref[i, i+1] != '*') & (query[i:i+1] != '_') | (query[i, i+1] != '*')) [Boolean]   ## NB: or
			###    isInsertion = ( (query[i:i+1] != '_') & (query[i, i+1] != '*') & (ref[i:i+1] != '_') | (ref[i, i+1] != '*')) [Boolean]  ## NB: or
			###
			for(my $i = 0; $i < length($ref); $i++)   ### why are we using 'ref' and 'query' instead of the variables above? 
			{
				if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*') and (substr($query, $i, 1) ne '-') and (substr($query, $i, 1) ne '*'))
				{
					$firstMatch = $i unless(defined $firstMatch);   ### why would 'firstMatch' be defined? 
					$lastMatch = $i;
				}
				
				my $isDeletion = (((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*')) and ((substr($query, $i, 1) eq '-') or (substr($query, $i, 1) eq '*')));
				my $isInsertion = (((substr($query, $i, 1) ne '-') and (substr($query, $i, 1) ne '*')) and ((substr($ref, $i, 1) eq '-') or (substr($ref, $i, 1) eq '*')));
				die if($isDeletion and $isInsertion);   ### kill script if DEL == INS [both Booleans]
				
				## if isDeletion == True....then what? if isInsertion == True...then what? 
				if($isDeletion)                ### these conditionals should be commented out
				{
					# warn Dumper("Deletions before start?", $ref, $query, $alignment->query->name, $alignment->cigar_str);
				}
				if($isInsertion)
				{
					# warn Dumper("One of the insertions!", $ref, $query, $alignment->query->name, $alignment->cigar_str);
				}	
			}
			die unless(defined $firstMatch);
			
			$ref = substr($ref, $firstMatch, $lastMatch - $firstMatch + 1);        ### I'm guessing here is where we skip regions?
			$query = substr($query, $firstMatch, $lastMatch - $firstMatch + 1);		
			
			if($alignment->query->name eq 'Korean.gi|1078261939|gb|LPVO02001249.1|')   ### checking for strange alignment names? Looks like should be removed for final script
			{
				# warn Dumper($alignment->query->name, $ref_preAll, $query_preAll, $ref, $query);
			}
			my $gaps_left_side = 0;    ### both at 0 to begin
			my $gaps_right_side = 0;
			
			for(my $i = 0; $i < length($ref); $i++)    
			{
				### PSEUDO CODE
				### for i in length(ref):
				###     if ((ref[i, i+1] != '-') & (ref[i, i+1] != '*')):
				###         break
				###     else:
				###         gaps_right_side += 1
				###
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
				if((substr($ref, $i, 1) ne '-') and (substr($ref, $i, 1) ne '*'))    ### if ((ref[i, i+1] != '-') & (ref[i, i+1] != '*')):
				{
					last;   ### 'last;' means 'break;'
				}
				else
				{
					$gaps_right_side++;
				}
			}
			
			my $starstar_left_side = 0;
			my $starstar_right_side = 0;
			for(my $i = 0; $i < length($ref); $i++)
			    ### PSEUDO CODE
				### for i in length(ref):
				###     if ((ref[i, i+1] != '*') & (query[i, i+1] != '*')):
				###        starstar_left_side += 1
				###     else:
				###         break
				###
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
			
			
			die unless($gaps_left_side == 0);    ### gaps_left_side == 0 or kill 
			die unless($gaps_right_side == 0);   ### gaps_right_side == 0 or kill 
			
			### PSEUDO CODE
			### ref = ref[int(gaps_left_side), int( len(ref) - int(gaps_left_side) - int(gaps_right_side) )]
			### query = query[int(gaps_left_side), int( len(query) - int(gaps_left_side) - int(gaps_right_side) )]
			###
			$ref = substr($ref, $gaps_left_side, length($ref) - $gaps_left_side - $gaps_right_side);
			$query = substr($query, $gaps_left_side, length($query) - $gaps_left_side - $gaps_right_side);		
			#$ref = substr($ref, 0, length($ref) - $gaps_right_side);
			#$query = substr($query, 0, length($query) - $gaps_right_side);		
			die unless(length($ref) == length($query));   ### ref and query should length(ref) == length(query)
			
			if($alignment->query->name eq 'CHM13.gi|953910992|gb|LDOC03004332.1|')   ### should be commented out
			{
				#die Dumper($alignment_start_pos, $ref, $query);
			}			
		
			my $ref_pos = $alignment_start_pos - 1;   ### cf. above, `my $alignment_start_pos = $alignment->start - 1; `
			my $running_gaps = 0;
			###
			### PSEUDOCODE
			### ref_pos = alignment_start_pos -1  
			### running_gaps = 0
			### for i in length(ref):
			###     c_ref = ref[i, i+1]
			###     if ((c_ref == '_') | (c_ref == '*'))
			###         ref_pos += 1   ### move to the next
			### 
			for(my $i = 0; $i < length($ref); $i++)     ### EVAN: I don't understand what this should do
			{
				my $c_ref = substr($ref, $i, 1);
				
				if(!(($c_ref eq '-') or ($c_ref eq '*')))
				{
					$ref_pos++;
				}
			}
			
			my $alignment_last_pos = $ref_pos - 1;
			my $alignment_info_aref = [$ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos];
			push(@{$alignments_starting_at{$alignment_start_pos}}, $alignment_info_aref);  ## add elements to array
			
			print D join("\t", $ref, $query, $alignment->query->name, $alignment_start_pos, $alignment_last_pos), "\n";   ### write out in VCF format
		}
	}
	else   ### testing
	{
		%alignments_starting_at = %alignments_starting_at_test;
	}
	
	# reference gaps *before* the i-th reference character
    ## PSEUDO CODE
    ### examine_gaps_n_alignment = 0
    ### for alignment_start_pos in alignments_starting_at:
    ###     for alignment 


	my $examine_gaps_n_alignment = 0;
	foreach my $alignment_start_pos (keys %alignments_starting_at)
	{
		foreach my $alignment (@{$alignments_starting_at{$alignment_start_pos}})
		{
			my $ref = $alignment->[0];
			
			my $start_pos = $alignment_start_pos - 1;
			my $ref_pos = $start_pos;
			my $running_gaps = 0;
			###
            ### running_gaps = 0
			### for i in length(ref):
			###    c_ref = ref[i, i+1]
			###    if ((c_ref == '-') | (c_ref == '*')):
			###        running_gaps += 1   ### move to the next    
			###    else:
			###        if (res_pos != start_pos):
			###            if gap_structure[ref_pos] == NULL: ##???
			###                gap_structure[ref_pos] = running_gaps
			###            else: 
			###                if gap_structure[$ref_pos] == running_gaps:   ### it already is assigned this? huh? I'm not sure how this would occur given the if statement
			###                    break:
			###                else: 
			###                    kill program, because "Gap structure mismatch at position $ref_pos - this is alignment etc. etc. "
			###
			for(my $i = 0; $i < length($ref); $i++)
			{
				my $c_ref = substr($ref, $i, 1);
				
				if(($c_ref eq '-') or ($c_ref eq '*'))
				{
					$running_gaps++;
				}
				else
				{    
					if(($ref_pos >= 5309528) and ($ref_pos <= 5309532))   ### EVAN: why these integers? 5309528 and 5309532?  Comment this out? 
					{
						# print "Gaps $ref_pos ". $alignment->[2] . ": " . $running_gaps . "\n";
					}	
					if($ref_pos != $start_pos)
					{
						if(not defined $gap_structure[$ref_pos])
						{
							$gap_structure[$ref_pos] = $running_gaps;
						}
						else
						{
							die "Gap structure mismatch at position $ref_pos - this is alignment $examine_gaps_n_alignment, have existing value $gap_structure[$ref_pos], want to set $running_gaps" unless($gap_structure[$ref_pos] == $running_gaps);
						}
					}
					$ref_pos++;
					$running_gaps = 0;
					$examine_gaps_n_alignment++;
				}
			}	
		}
	}
			
	close(D);
	
	print "Have loaded $n_alignments alignments -- $fn_for_BAM2VCF.\n";

	my $cmd = qq($bin_BAM2VCF --input $fn_for_BAM2VCF --referenceSequenceID $referenceSequenceID);
	
	print "Now executing: $cmd\n";	
	
	my $output_file = $fn_for_BAM2VCF . '.VCF';
	
	# todo
	#if(system($cmd))
	#{
	#	die "Command $cmd failed";
	#}
	
	unless(-e $output_file)
	{
		die "File $output_file not existing";
	}
	
	open(CHROUT, '<', $output_file) or die "Cannot open $output_file";
	while(<CHROUT>)
	{
		print OUT $_;
	}
}

close(OUT);

print "\n\nOK\n\n";


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
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}


### I don't see this function ever used....

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
}

### I don't see this function ever used either...


sub printHaplotypesAroundPosition
{
	my $posI = shift;
	my $alignments_starting_at_href = shift;
	
	my @aroundPositions = grep {$_ >= 0} (($posI - 2) .. ($posI + 2));
	foreach my $startPos (keys %{$alignments_starting_at_href})
	{
		foreach my $alignment_aref (@{$alignments_starting_at_href->{$startPos}})
		{
			my $apparent_startPos = $alignment_aref->[3];
			if($apparent_startPos != $startPos)
			{
				warn Dumper("Apparent / startpos mismatch", $startPos, $apparent_startPos, $alignment_aref->[2]);
			}
			my $stopPos = $alignment_aref->[4];
			die unless($stopPos >= $startPos);
			my $interesting = 0;
			foreach my $interestingPos (@aroundPositions)
			{
				if(($interestingPos >= $startPos) and ($interesting <= $stopPos))
				{
					$interesting = 1;
				}
			}
			
			my $ref = $alignment_aref->[0];
			my $query = $alignment_aref->[1];
			
			if($interesting)
			{
				my %gt_per_position;
				my $ref_pos = $startPos - 1;
				my $running_allele = '';
				for(my $i = 0; $i < length($ref); $i++)
				{
					my $c_ref = substr($ref, $i, 1);
					my $c_query = substr($query, $i, 1);
					
					if(($c_ref eq '-') or ($c_ref eq '*'))
					{
						$running_allele .= $c_query;
					}
					else
					{
						if($running_allele)
						{
							$gt_per_position{$ref_pos} = $running_allele;
						}
						
						$running_allele = $c_query;
						$ref_pos++;						
					}
				}
				if($running_allele)
				{
					$gt_per_position{$ref_pos} = $running_allele;
				}
						
				print "Positions $alignment_aref->[2] \n";
				foreach my $interestingPos (@aroundPositions)
				{
					if(exists $gt_per_position{$interestingPos})
					{
						print "\t", $interestingPos, "\t", $gt_per_position{$interestingPos}, "\n";
					}
				}
			}
		}
	}
}

