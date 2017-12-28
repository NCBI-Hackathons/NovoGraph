#!/usr/bin/perl

use strict;

my $inBAM = $ARGV[0];
unless($inBAM and ($inBAM =~ /\.bam/))
{
	die "Usage: totalCoverage.pl BAMFILE";
}

my $n_alignments = 0;
my $totalBases = 0;
open(SAMTOOLS, '-|', "samtools view $inBAM") or die "Cannot open $inBAM with samtools view";
while(<SAMTOOLS>)
{
	my $line = $_;
	my @fields = split(/\t/, $line);
	my $sequence = $fields[9];
	$n_alignments++;
	$totalBases += length($sequence);
}

print "\nCoverage summary $inBAM\n";
print "\tTotal alignments   : $n_alignments \n";
print "\tTotal bases        : $totalBases \n";
print "\tAvg whole-gen. cov.: ", ($totalBases/3.2e9), "\n";