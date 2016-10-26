#!/usr/bin/perl

use strict;
use Data::Dumper;
use Bio::DB::Sam;
use Getopt::Long;   
use List::MoreUtils qw/mesh/;

$| = 1;

# Example command:
# ./checkWindowCorrectness.pl --contigsFasta /home/data/contigs/AllContigs.fa --fileWithAllGapWindows ../forMafft2/_alignments_inWindow_onlyGaps --windowFilesDirectory ../forMafft2 --windowsInfo ../forMafft2/_windowsInfo

my $contigsFasta;
my $fileWithAllGapWindows;
my $windowFilesDirectory;
my $windowsInfo;

GetOptions (
	'contigsFasta:s' => \$contigsFasta, 
	'fileWithAllGapWindows:s' => \$fileWithAllGapWindows, 
	'windowsInfo:s' => \$windowsInfo, 
	'windowFilesDirectory:s' => \$windowFilesDirectory,	
);

die "Please specify --contigsFasta" unless($contigsFasta);
die "Please specify --fileWithAllGapWindows" unless($fileWithAllGapWindows);
die "Please specify --windowFilesDirectory" unless($windowFilesDirectory);
die "Please specify --windowsInfo" unless($windowsInfo);

die "--contigsFasta $contigsFasta not existing" unless(-e $contigsFasta);
die "--fileWithAllGapWindows $fileWithAllGapWindows not existing" unless(-e $fileWithAllGapWindows);
die "--windowFilesDirectory $windowFilesDirectory not existing" unless(-e $windowFilesDirectory);
die "--windowsInfo $windowsInfo not existing" unless(-e $windowsInfo);

print "Read $contigsFasta\n";
my $contigs_href = readFASTA($contigsFasta, 0);
print "\tdone.\n";

my @fasta_files_in_dir = glob($windowFilesDirectory . '/*');

my %windows_all_gaps;
{
	my $n_entries_all_gaps = 0;
	open(WINDOWSALLGAPS, '<', $fileWithAllGapWindows) or die "Cannot open $fileWithAllGapWindows";
	my $headerLine = <WINDOWSALLGAPS>;
	chomp($headerLine);
	my @header_fields = split(/\t/, $headerLine);
	while(<WINDOWSALLGAPS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless($#line_fields == $#header_fields);
		my %line_hash = (mesh @header_fields, @line_fields);
		die unless($line_hash{alignedSequenceID});
		die unless($line_hash{referenceContigID});
		die unless(defined $line_hash{windowI});
		$windows_all_gaps{$line_hash{referenceContigID}}{$line_hash{windowI}}{$line_hash{alignedSequenceID}}++;  
		$n_entries_all_gaps++;
	}
	close(WINDOWSALLGAPS);
	print "Have $n_entries_all_gaps window / contig combinations that are expected as all-gap.\n";
}

open(WINDOWS, '<', $windowsInfo) or die "Cannot open $windowsInfo";
my $headerLine = <WINDOWS>;
chomp($headerLine);
my @header_fields = split(/\t/, $headerLine);
my %runningContigSequences;
my $n_agree = 0;
my $n_disagree = 0;
while(<WINDOWS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	die unless($#line_fields == $#header_fields);
	my %line_hash = (mesh @header_fields, @line_fields);
	die unless($line_hash{referenceContigID});
	die unless(defined $line_hash{windowI});
	
	my $expectedFile = $windowFilesDirectory . '/' . $line_hash{referenceContigID} . '_' . $line_hash{windowI} . '.fa';
	unless(-e $expectedFile)
	{
		die "Did not found expected file: $expectedFile";
	}
	
	my $fasta_window = readFASTA($expectedFile);
	
	my %running_ids_not_seen = map {$_ => 1} keys %runningContigSequences;
	foreach my $gapSequenceID (keys %{$windows_all_gaps{$line_hash{referenceContigID}}{$line_hash{windowI}}})
	{
		if($windows_all_gaps{$line_hash{referenceContigID}}{$line_hash{windowI}}{$gapSequenceID})
		{
			$running_ids_not_seen{$gapSequenceID} = 0;
		}
	}
	
	foreach my $sequenceID (keys %$fasta_window)
	{
		if($windows_all_gaps{$line_hash{referenceContigID}}{$line_hash{windowI}}{$sequenceID})
		{
			die "Sequence ID $sequenceID is supposed to be gap in window $line_hash{referenceContigID} , $line_hash{windowI}";
		}
		
		$running_ids_not_seen{$sequenceID} = 0;
		$runningContigSequences{$sequenceID} .= $fasta_window->{$sequenceID};
	}
	
	foreach my $finishedSequenceID (keys %running_ids_not_seen)
	{
		next unless($running_ids_not_seen{$finishedSequenceID});
		
		unless(exists $contigs_href->{$finishedSequenceID})
		{
			warn "Don't have truth sequence for $finishedSequenceID";
			next;
		}
		my $trueSequence = $contigs_href->{$finishedSequenceID};
		my $trueSequence_revCmp = reverseComplement($trueSequence);
		
		my $sequenceFromWindows = $runningContigSequences{$finishedSequenceID};
		if(($trueSequence eq $sequenceFromWindows) or ($trueSequence_revCmp eq $sequenceFromWindows))
		{
			# print "OK\n";
		}
		else
		{
			print "Disagreement!\n";
			print "\tlen trueSequence: ", length($trueSequence), "\n";
			print "\tlen sequenceFromWindows: ", length($sequenceFromWindows), "\n";
			print "\ttrueSequence:        ", $trueSequence, "\n";
			print "\tsequenceFromWindows: ", $sequenceFromWindows, "\n";
			$n_disagree++;
			exit if($n_disagree > 10);
		}
		
		delete $runningContigSequences{$finishedSequenceID};
	}	
}
close(WINDOWS);

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

