#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

## Usage:
## checkMAFFT_input_and_output.pl --MAFFTdir <path to MAFFT directory used in BAM2MAFFT.pl and CALLMAFFT.pl>
##                                --contigLengths <path to file with contig names/lengths>
##                                --preMAFFTBAM <path to output BAM from outputFile, FIND_GLOBAL_ALIGNMENTS.pl>
##                                --finalOutputCRAM <path to CRAM>
##                                --fas2bam_path <path to script 'fas2bam.pl'>
##                                --samtools_path <path to SAMtools executable for fas2bam.pl>
##                                --bamheader <path to file containing header for BAM file for fas2bam.pl>
##
##
## Example command:
## ./checkMAFFT_input_and_output.pl --MAFFTdir /intermediate_files/forMAFFT/
##                                  --contigLengths /intermediate_files/postGlobalAlignment_readLengths
##                                  --preMAFFTBAM /intermediate_files/forMAFFT.bam
##                                  --finalOutputCRAM /intermediate_files/combined.cram
##                                  --fas2bam_path /intermediate_files/fas2bam.pl
##                                  --samtools_path /usr/local/bin/samtools
##                                  --bamheader windowbam.header.txt


my $contigLengths;
my $MAFFTdir;
my $preMAFFTBAM;
my $finalOutputCRAM;
my $fas2bam_path;
my $samtools_path;
my $bamheader;
GetOptions (
	'MAFFTdir:s' => \$MAFFTdir, 
	'contigLengths:s' => \$contigLengths, 
	'preMAFFTBAM:s' => \$preMAFFTBAM,
	'finalOutputCRAM:s' => \$finalOutputCRAM, 
	'fas2bam_path:s' => \$fas2bam_path,
	'samtools_path:s' => \$samtools_path,
	'bamheader:s' => \$bamheader,
);

die unless($contigLengths);
die unless($MAFFTdir);
die unless($preMAFFTBAM);
die unless($finalOutputCRAM);
die unless($fas2bam_path);
die unless($samtools_path);
die unless($bamheader);

## Keep track of all errors
my $number_total_errors = 0;

my $readWindowsInfo = read_windowbams_info($MAFFTdir);

my %expectedLengts;
open(L, '<', $contigLengths) or die "Cannot open --contigLengths $contigLengths";
while(<L>)
{
	my $l = $_; 
	chomp($l);
	my @f = split(/\t/, $l);
	$expectedLengts{$f[0]} = $f[1];
}
close(L);

my %contigs_fromFASTA;
my %contigs_fromFASTA_fileSources;
my %contigs_fromMFA;
my %contigs_fromMFA_fileSources;
my %contigs_from_singleBAMs;
foreach my $chr (keys %$readWindowsInfo)
{
	print "Reading $chr... - ", scalar(@{$readWindowsInfo->{$chr}}), " windows.\n";
	my $windowI = 0;
	foreach my $window (@{$readWindowsInfo->{$chr}})
	{
		$windowI++;
		if(($windowI % 10000) == 0)
		{
			print "\tWindow $windowI ...\n";
		}
		my $faFile = $window->{fasta};
		my $mfaFile = $window->{mfa};
		die unless($faFile);
		die unless($mfaFile);
		my $fa_href = readFASTA($faFile);
		my %fa_href_seq2;
		my %lengths_fa;
		foreach my $seqID (keys %$fa_href)
		{
			next if($seqID =~ /^ref/);
			my $seqID2 = $seqID;
			$seqID2 =~ s/_FIRST//;
			
			my $seq_noGaps = $fa_href->{$seqID};
			$seq_noGaps =~ s/[_\-]//g;
					
			$fa_href->{$seqID} =~ s/-//g;
			$contigs_fromFASTA{$seqID2} .= $seq_noGaps;
			$fa_href_seq2{$seqID2} = $seq_noGaps;
			$contigs_fromFASTA_fileSources{$seqID2}{$faFile}++;
			$lengths_fa{$seqID2} = length($seq_noGaps);
		}
		
		my $mfa_href = readFASTA($mfaFile);
		my $mfaLength;
		foreach my $seqID (keys %$mfa_href)
		{
			next if($seqID =~ /^ref/);
			my $seqID2 = $seqID;
			$seqID2 =~ s/_FIRST//;
			
			my $seq = $mfa_href->{$seqID};
			unless(defined $mfaLength)
			{
				$mfaLength = length($seq);
			}
			warn Dumper("Problem with sequence lengths", $mfaFile, $mfaLength, length($seq), $seqID) unless(length($seq) == $mfaLength);
			unless(length($seq) == $mfaLength)
			{
				$number_total_errors++; ## keep record of total number of errors
			}
			
			my $seq_noGaps = $seq;
			$seq_noGaps =~ s/[_\-]//g;
		
			$contigs_fromMFA{$seqID2} .= $seq_noGaps;
			$contigs_fromMFA_fileSources{$seqID2}{$mfaFile}++;
		}

		my $BAM = $window->{bam};
		die unless(-e $BAM);
		
		my %lengths_BAM;
		my %BAM_href_seq2;
		open(SINGLEBAM, $samtools_path, " view $BAM |") or die "Cannot view $BAM $!";
		while(<SINGLEBAM>)
		{
			my $l = $_;
			chomp($l);
			next unless($l);
			my @f = split(/\t/, $l);
			die unless(scalar(@f) > 3);
			my $seq = $f[9];
			$seq = '' if($seq eq '*');
			my $readID = $f[0];
			my $readID2 = $readID;
			$readID2 =~ s/_FIRST//;
			
			die "Invalid sequence for $readID $readID2 in $BAM" if($seq =~ /[\-_]/);
			
			$contigs_from_singleBAMs{$readID2} .= $seq;
			$BAM_href_seq2{$readID2} = $seq;
			$lengths_BAM{$readID2} = length($seq);
		}
		close(SINGLEBAM);
		  
		my %joint_keys = map {$_ => 1} ((keys %lengths_fa), (keys %lengths_BAM));
		my $printInfo = sub {
			print "Print error report - differing lengths|\n";
			foreach my $seqID (keys %joint_keys)
			{
				print "\t", $seqID, "\n";
				print "\t\tFASTA: ", $lengths_fa{$seqID}, " ", $fa_href_seq2{$seqID}, "\n";
				print "\t\tBAM:   ", $lengths_BAM{$seqID}, " ", $BAM_href_seq2{$seqID}, "\n";
				unless($lengths_fa{$seqID} == $lengths_BAM{$seqID})
				{
					print "\t\t!!!!\n";
				}
			}
			print Dumper($window), "\n";			
		};	
		
		my $redo_BAM = 0;
		foreach my $seqID (keys %joint_keys)
		{
			unless((defined $lengths_fa{$seqID}) and (defined $lengths_BAM{$seqID}) and ($lengths_fa{$seqID} == $lengths_BAM{$seqID}))
			{
				$redo_BAM = 1;
				$printInfo->();
				$number_total_errors++;  ## keep record of total number of errors
			}
		}	

		$redo_BAM = 0;
		if($redo_BAM)
		{
			my $cmd_re_execute = qq(perl $fas2bam_path --input $window->{mfa} --ref "ref" --output $window->{bam} --bamheader $bamheader --samtools_path $samtools_path);
			print $cmd_re_execute, "\n";
			system($cmd_re_execute) and die "Command $cmd_re_execute failed!";		
		}
	}
}

my %preMAFFT_BAM_lengths;
my %preMAFFT_BAM_sequence;
open(PREMAFFT, $samtools_path, " view $preMAFFTBAM |") or die "Cannot view $preMAFFTBAM";
while(<PREMAFFT>)
{
	my $l = $_;
	chomp($l);
	next unless($l);
	my @f = split(/\t/, $l); 
	die unless(scalar(@f) > 3);
	my $seq = $f[9];
	my $readID = $f[0];
	$preMAFFT_BAM_lengths{$readID} = length($seq);
	$preMAFFT_BAM_sequence{$readID} = $seq;
}
close(PREMAFFT);


my %finalOutput_BAM_lengths;
my %finalOutput_BAM_sequence;
if($finalOutputCRAM =~ /\.cram$/)
{
	open(OUTPUT, $samtools_path, " view $finalOutputCRAM |") or die "Cannot view $finalOutputCRAM";
	while(<OUTPUT>)
	{
		my $l = $_;
		chomp($l);
		next unless($l);
		my @f = split(/\t/, $l);
		die unless(scalar(@f) > 3);
		my $seq = $f[9];
		my $readID = $f[0];
		$finalOutput_BAM_lengths{$readID} = length($seq);
		$finalOutput_BAM_sequence{$readID} = $seq;
	}
	close(OUTPUT);
}
else
{
	die unless($finalOutputCRAM =~ /\.sam$/);
	open(OUTPUT, '<', $finalOutputCRAM) or die "Cannot open $finalOutputCRAM";
	while(<OUTPUT>)
	{
		my $l = $_;
		chomp($l);
		next unless($l);
		my @f = split(/\t/, $l);
		die unless(scalar(@f) > 3);
		my $seq = $f[9];
		my $readID = $f[0];
		$finalOutput_BAM_lengths{$readID} = length($seq);
		$finalOutput_BAM_sequence{$readID} = $seq;
	}
	close(OUTPUT);
}

my %combinedKeys = map {$_ => 1} ((keys %preMAFFT_BAM_lengths), (keys %contigs_fromFASTA));

foreach my $k (keys %combinedKeys)
{
	my $l_expected = (exists $expectedLengts{$k}) ? $expectedLengts{$k} : "NOT_PRESENT";
	my $l_observed = (exists $contigs_fromFASTA{$k}) ? length($contigs_fromFASTA{$k}) : "NOT_PRESENT";
	my $l_BAM_preMAFFT = (exists $preMAFFT_BAM_lengths{$k}) ? $preMAFFT_BAM_lengths{$k} : "NOT_PRESENT";
	my $l_MFA_postMafft = (exists $contigs_fromMFA{$k}) ? length($contigs_fromMFA{$k}) : "NOT_PRESENT";
	my $l_singleBAM_postMafft = (exists $contigs_from_singleBAMs{$k}) ? length($contigs_from_singleBAMs{$k}) : "NOT_PRESENT";
	my $l_finalOutput = (exists $finalOutput_BAM_sequence{$k}) ? length($finalOutput_BAM_sequence{$k}) : "NOT_PRESENT";
	print join("\t", $k, $l_expected, $l_BAM_preMAFFT, $l_MFA_postMafft, $l_singleBAM_postMafft, $l_finalOutput), "\n";
	if($l_expected ne $l_finalOutput)
	{
	        $number_total_errors++;  ## keep record of total number of errors
		my $contigFromBAM = $preMAFFT_BAM_sequence{$k};
		print "Discrepancy: $k \n";
		print "\tExpected       : $l_expected \n";
		print "\tPre-MAFFT, BAM : $l_BAM_preMAFFT \n";
		print "\tPost-MAFFT, mfa: $l_MFA_postMafft \n";
		print "\tPost-MAFFT, BAM: $l_singleBAM_postMafft \n";
		print "\tCombined final : $l_finalOutput \n";
		print "\tFile sources:\n";
		foreach my $k (keys %{$contigs_fromFASTA_fileSources{$k}})
		{
			print "\t - ", $k, "\n";
		}
	}
}

if($n_printed_errors)
{
	die "A total of $number_total_errors issues were detected. Please see output above.";
}

print "Done. No errors found!\n";


sub read_windowbams_info {
    my $dir = shift;

    my %entry_hash = ();
    open WINDOWS, "$dir/_windowsInfo"
        or die "Couldn\'t open $dir/_windowsInfo: $!\n";
    my ($lastentry, $laststart); # make sure sorted!
    while (<WINDOWS>) {
        next if (/^referenceContigID/);
        if (/^(\S+)\s(\S+)\s(\d+)\s(\d+)\s(\d+)/) {
            my ($entry, $chrDir, $windowid, $start, $end) = ($1, $2, $3, $4);
			my $bamFile = "$dir/$chrDir/$entry\_$windowid.bam";
			my $faFile = "$dir/$chrDir/$entry\_$windowid.fa";
			my $mfaFile = "$dir/$chrDir/$entry\_$windowid.mfa";
			die "BAM file $bamFile not present" unless(-e $bamFile);
			die "FASTA file $faFile not present" unless(-e $faFile);
			die "MFA file $mfaFile not present" unless(-e $mfaFile);
            push @{$entry_hash{$entry}}, {bam => $bamFile,
										  fasta => $faFile,
										  mfa => $mfaFile,
                                          offset => $start};
            if ($lastentry && ($entry eq $lastentry) && ($start < $laststart)) {
                die "Windows in _windowInfo are unsorted! ($entry:$start is less than $laststart)\n";
            }
            $lastentry = $entry;
            $laststart = $start;
        }
        else {
            die "Illegal format in _windowsInfo file:\n$_";
        }
    }
    close WINDOWS;

    return {%entry_hash};
}

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

