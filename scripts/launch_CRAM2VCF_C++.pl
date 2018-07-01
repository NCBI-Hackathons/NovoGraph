#!/usr/bin/env perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all/;
use List::MoreUtils qw/mesh/;
use Bio::DB::HTS;

## Usage:
## launch_CRAM2VCF_C++.pl --output <path to VCF created by CRAM2VCF.pl>
##
## Example command:
## 	./launch_CRAM2VCF_C++.pl --output VCF/graph_v2.vcf



$| = 1;

my $output;

GetOptions (
	'output:s' => \$output,
);

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
	
	my $VCF = $inputFile . '.VCF';
	my $doneFile = $VCF . '.done';
	if(-e $doneFile)
	{
		open(DONE, '<', $doneFile) or die "Cannot open $doneFile";
		my $done = <DONE>;
		chomp($done);
		$done = (length($done) > 0) ? substr($done, 0, 1) : 0;
		close(DONE);
		
		if($done)
		{
			$files_done++;
			next;
		} 
		else 
		{
			unlink($doneFile) or die "Cannot delete $doneFile";		
		}
	}
	
	push(@inputFiles, $inputFile);
	push(@commands, $line);
	
	
}
close(CMDS);

print "Files done already: $files_done -- delete $output*.done if you want to redo these!\n";

my @command_batches = ([]);
my $runningSize = 0;
my $runningCommands = 0;
my $totalCommands = 0;
for(my $i = 0; $i <= $#commands; $i++)
{
	if(($runningCommands >= 20) or ($runningSize >= 1e6))
	{
		push(@command_batches, []);	
		$runningCommands = 0;
		$runningSize = 0;
	}
	push(@{$command_batches[$#command_batches]}, $commands[$i]);
	$runningCommands++;
	$totalCommands++;
	$runningSize += (-s $inputFiles[$i]);
}	

if($totalCommands == 0)
{
	print "\nAll done.\n\n";
}
else
{
	print "\nTotal command batches: ", scalar(@command_batches), "\n\n";

	foreach my $command_batch (@command_batches)
	{
		my $combined_command = join(';', @$command_batch);
		my $pid = fork;
		die "fork failed" unless defined $pid;
		if ($pid == 0) {
			system($combined_command) and die "Could not execute command: $combined_command";
			exit;
		}	
	}

	print "\n\nProcesses launched.\n";
}
