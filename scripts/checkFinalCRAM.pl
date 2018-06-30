#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::Sam;
$| = 1;

my $contigLengths;
my $MAFFTdir;
my $preMAFFTBAM;
my $finalOutputCRAM;
GetOptions (
	'MAFFTdir:s' => \$MAFFTdir, 
	'contigLengths:s' => \$contigLengths, 
	'preMAFFTBAM:s' => \$preMAFFTBAM, 
	'finalOutputCRAM:s' => \$finalOutputCRAM, 
);

die unless($contigLengths);
die unless($MAFFTdir);
die unless($preMAFFTBAM);
die unless($finalOutputCRAM);

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

my %preMAFFT_BAM_lengths;
my %preMAFFT_BAM_sequence;
open(PREMAFFT, "/data/projects/phillippy/software/samtools-1.4/samtools view $preMAFFTBAM |") or die "Cannot view $preMAFFTBAM";
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
open(OUTPUT, "/data/projects/phillippy/software/samtools-1.4/samtools view $finalOutputCRAM |") or die "Cannot view $finalOutputCRAM";
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


# my %combinedKeys = map {$_ => 1} ((keys %expectedLengts), (keys %contigs_fromFASTA), (keys %preMAFFT_BAM_lengths));
# my %combinedKeys = map {$_ => 1} ((keys %preMAFFT_BAM_lengths), (keys %contigs_fromFASTA));
my %combinedKeys = map {$_ => 1} ((keys %preMAFFT_BAM_lengths));

my $perfect = 0;
my $issues =  0;
foreach my $k (keys %combinedKeys)
{
	my $l_expected = (exists $expectedLengts{$k}) ? $expectedLengts{$k} : "NOT_PRESENT";
	# my $l_observed = (exists $contigs_fromFASTA{$k}) ? length($contigs_fromFASTA{$k}) : "NOT_PRESENT";
	my $l_BAM_preMAFFT = (exists $preMAFFT_BAM_lengths{$k}) ? $preMAFFT_BAM_lengths{$k} : "NOT_PRESENT";
	# my $l_MFA_postMafft = (exists $contigs_fromMFA{$k}) ? length($contigs_fromMFA{$k}) : "NOT_PRESENT";
	# my $l_singleBAM_postMafft = (exists $contigs_from_singleBAMs{$k}) ? length($contigs_from_singleBAMs{$k}) : "NOT_PRESENT";
	my $l_finalOutput = (exists $finalOutput_BAM_sequence{$k}) ? length($finalOutput_BAM_sequence{$k}) : "NOT_PRESENT";
	# print join("\t", $k, $l_expected, $l_BAM_preMAFFT, $l_MFA_postMafft, $l_singleBAM_postMafft, $l_finalOutput), "\n";
	print join("\t", $k, $l_expected, $l_BAM_preMAFFT, $l_finalOutput), "\n";
	if($l_expected ne $l_finalOutput)
	{
		my $contigFromBAM = $preMAFFT_BAM_sequence{$k};
		print "Discrepancy: $k \n";
		print "\tExpected       : $l_expected \n";
		print "\tPre-MAFFT, BAM : $l_BAM_preMAFFT \n";
		# print "\tPost-MAFFT, mfa: $l_MFA_postMafft \n";
		# print "\tPost-MAFFT, BAM: $l_singleBAM_postMafft \n";
		print "\tCombined final : $l_finalOutput \n";
		#print "\tFile sources:\n";
		#foreach my $k (keys %{$contigs_fromFASTA_fileSources{$k}})
		#{
		#	print "\t - ", $k, "\n";
		#}
		$issues++;
	}
	else
	{
		$perfect++;
	}
}

print "Done.\n";
print "\tPerfect: $perfect \n";
print "\tIssues : $issues \n";

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

