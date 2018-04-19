#!/usr/bin/env perl

use strict;

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

$| = 1;

my $output;

GetOptions (
	'output:s' => \$output,
);

my $chr = 'chr20';
my $fn_output_distribution = $output . '.CRAM2VCF_INDELLengths';
if($chr)
{
	$fn_output_distribution .= '_' . $chr;
}

my $files_done = 0;
my @commands;
my @inputFiles;
my $fn_cmds = $output . '_CRAM2VCF_commands.txt';
open(CMDS, '<', $fn_cmds) or die "Cannot open $fn_cmds";
while(<CMDS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	$line =~ s/\&//g;
	
	die unless($line =~ /--input (\S+?) --referenceSequenceID/);	
	my $inputFile = $1;	
	
	next unless($inputFile =~ /chr(\d+|X|Y)$/);
	if(defined $chr)
	{
		next unless($inputFile =~ /$chr/);
	}
	push(@inputFiles, $inputFile);
}
close(CMDS);

my $n_matches = 0;
my $n_mismatches = 0;

my %INDELs;
foreach my $inputFile (@inputFiles)
{
	print "Processing $inputFile ...\n";
	open(F, '<', $inputFile) or die "Cannot open $inputFile";
	<F>;
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @files = split(/\t/, $line);
		my $ref = $files[0];
		my $query = $files[1];
		die unless(length($ref) == length($query));
		die unless((substr($ref, 0, 1) ne '-') and (substr($query, 0, 1) ne '-'));
		my $running_insertion = 0;
		my $running_deletion = 0;
		for(my $i = 0; $i < length($ref); $i++)
		{
			my $C_ref = substr($ref, $i, 1);
			my $C_query = substr($query, $i, 1);
			next if(($C_ref eq '-') and ($C_query eq '-'));
			if($C_ref eq '-')
			{
				if($running_deletion)
				{
					$INDELs{-1 * $running_deletion}++;
					$running_deletion = 0;
				}
				$running_insertion++;
			}
			elsif($C_query eq '-')
			{
				if($running_insertion)
				{
					$INDELs{$running_insertion}++;
					$running_insertion = 0;
				}			
				$running_deletion++;
			}
			else
			{
				die Dumper("Double gap?", $C_ref, $C_query) unless(($C_ref ne '-') and ($C_query ne '-'));
				if($running_deletion)
				{
					$INDELs{-1 * $running_deletion}++;
					$running_deletion = 0;
				}	
				if($running_insertion)
				{
					$INDELs{$running_insertion}++;
					$running_insertion = 0;
				}	
				if($C_ref eq $C_query)
				{
					$n_matches++;
				}
				else
				{
					$n_mismatches++;
				}
			}
		}
	}
	close(F);
}

print "Matches   : $n_matches \n";
print "Mismatches: $n_mismatches \n";
print "Indel lengths: ", scalar(keys %INDELs), "\n";

open(F, '>', $fn_output_distribution) or die "Cannot open $fn_output_distribution";
foreach my $length (sort {$a <=> $b} keys %INDELs)
{
	print F join("\t", $length, $INDELs{$length}), "\n";
}
close(F);

print "\nProduced file $fn_output_distribution\n\n";

