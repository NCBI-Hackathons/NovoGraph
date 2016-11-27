#!/usr/bin/perl

# ./dealWithTooManyCIGAROperations.pl --input ../../intermediate_files/forMAFFT.bam.sam --output ../../intermediate_files/forMAFFT.bam.sam.filtered

use strict;
use warnings;
use Getopt::Long;   

# example command: 
my $input;
my $output;

GetOptions (
	'input:s' => \$input, 
	'output:s' => \$output, 
);

die "Please specify --input" unless($input);
die "--input $input not existing" unless(-e $input);
die "Please specify --output" unless($output);

open(IN, '<', $input) or die "Cannot open $input";
open(OUT, '>', $output) or die "Cannot open $output";
my $n_removed = 0;
while(<IN>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	if(substr($line, 0, 1) eq '@')
	{
		print OUT $line, "\n";
	}
	else
	{
		my @fields = split(/\t/, $line);
		die unless(scalar(@fields) >= 11);
		my $CIGAR = $fields[5];
		my $CIGARoperations = 0;
		while($CIGAR)
		{
			die unless($CIGAR =~ /^(\d+)(\w)/);
			my $n = $1;
			my $op = $2;
			die unless(($op eq 'M') or ($op eq 'D') or ($op eq 'I'));
			$CIGARoperations++;
			$CIGAR =~ s/^(\d+)(\w)//;
		}
		if($CIGARoperations <= 65530)
		{
			print OUT $line, "\n";
		}
		else
		{
			$n_removed++;
		}
	}
}

print "\nRemoved $n_removed alignments.\n";


sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
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
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= uc($line);
		}
	}	
	close(F);
		
	return \%R;
}

sub writeFASTA
{
	my $file = shift;
	# print "Writing $file\n";
	my $href = shift;
	open(F, '>', $file) or die "Cannot open $file";
	foreach my $key (keys %$href)
	{
		my $seq = $href->{$key};
		print F '>', $key, "\n";
		# print "\t", $key, "\t", length($seq), "\n";
		while($seq)
		{
			my $toPrint;
			if(length($seq) > 50)
			{
				$toPrint = substr($seq, 0, 50);
				substr($seq, 0, 50) = '';
			}
			else
			{
				$toPrint = $seq;
				$seq = '';
			}	
			print F $toPrint, "\n";
		}
	}
	close(F);	
}
