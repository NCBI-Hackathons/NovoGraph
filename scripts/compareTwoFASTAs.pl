#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
use Set::IntervalTree;
$| = 1;

# Example command
# ./compareTwoFASTAs.pl --f1 /home/devsci7/globalize_windowbams/global_multiple_alignments.frombam.fasta --f2 /home/data/contigs/AllContigs.fa

my $f1;
my $f2;
GetOptions (
	'f1:s' => \$f1, 
	'f2:s' => \$f2, 
);

die "Please specify --f1" unless($f1);
die "Please specify --f2" unless($f2);
die "--f1 $f1 not existing" unless(-e $f1);
die "--f2 $f2 not existing" unless(-e $f2);

print "Loading f1 into memory, then iteratively process f2...\n";

my $f1_href = readFASTA($f1);

my $f1_exclusive = 0;
my $f2_exclusive = 0;
my $shared_identical = 0;
my $shared_different = 0;

my %saw_f2;
my $process_sequence_from_f2 = sub {
	my $sequenceID = shift;
	my $sequence = shift;
	
	die unless($sequenceID);
	die unless(defined $sequence);
	
	$saw_f2{$sequenceID}++;
	
	if(exists $f1_href->{$sequenceID})
	{
		if(($f1_href->{$sequenceID} eq $sequence) or ($f1_href->{$sequenceID} eq reverseComplement($sequence)))
		{
			$shared_identical++;
		}
		else
		{
			$shared_different++;
		}
	}
	else
	{
		$f2_exclusive++;
	}
};

my $currentSequenceID;
my $runningSequence;
open(F2, '<', $f2) or die "Cannot open $f2";
while(<F2>)
{
	my $line = $_;
	chomp($line);
	$line =~ s/[\n\r]//g;
	if(substr($line, 0, 1) eq '>')
	{
		if(defined $runningSequence)
		{
			$process_sequence_from_f2->($currentSequenceID, $runningSequence);
		}
		$currentSequenceID = substr($line, 1);
		$currentSequenceID =~ s/\s.+//;
		$runningSequence = '';
	}
	else
	{
		$runningSequence .= $line;
	}
}	
close(F2);

if(defined $runningSequence)
{
	$process_sequence_from_f2->($currentSequenceID, $runningSequence);
}

foreach my $f1ID (keys %$f1_href)
{
	unless($saw_f2{$f1ID})
	{
		$f1_exclusive++;
	}
}

print "Exclusive to f1 [${f1}]: $f1_exclusive", "\n";
print "Exclusive to f2 [${f2}]: $f2_exclusive", "\n";
print "Shared: ", ($shared_identical + $shared_different), "\n";
print "\tIdentical: ", $shared_identical, "\n";
print "\tDifferent: ", $shared_different, "\n";


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

sub reverseComplement
{
	my $kMer = shift;
	$kMer =~ tr/ACGT/TGCA/;
	return reverse($kMer);
	return $kMer;
}




