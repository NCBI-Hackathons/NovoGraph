#!/usr/bin/env perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   
use List::Util qw/max all shuffle/;
use List::MoreUtils qw/mesh/;
use Errno;
require POSIX;    # provides WNOHANG

$SIG{CHLD} = \&reap_kids;

## Usage:
## launch_CRAM2VCF_C++.pl --prefix <path to output VCF without .VCF at the end> --number_processes <integer of total number of processes to run simultaneously>
##
## Example command:
## ./launch_CRAM2VCF_C++.pl --prefix VCF/graph_v2.vcf --number_processes 3

$| = 1;

my $output;
my $number_processes;

GetOptions (
	'prefix:s' => \$output,
	'number_processes:s' => \$number_processes
);

die "Please specify --prefix" unless($output);
## die "Please specify the total number of processes to be simultaneously run --processes" unless($processes);

## if undefined, set 
if (not defined $number_processes){
	$number_processes = 2;
}

die "Please specify a positive integer for input parameter `--number_processes` " unless ( $number_processes=~ /^[+]?\d+$/ && $number_processes != 0  );


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
	
	die "Weird line - $line - in file $fn_cmds" unless($line =~ /--input (\S+?) --referenceSequenceID/);	
	my $inputFile = $1;
	
	my $VCF = $inputFile . '.VCF';
	my $doneFile = $inputFile . '.done';
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

print "Files done already: $files_done -- delete $output*.done if you want to redo these! Remaining: ", scalar(@commands), "\n";

my @command_batches = ([]);
my $runningSize = 0;
my $runningCommands = 0;
my $totalCommands = 0;
my @indices = (0 .. $#commands);
@indices = shuffle(@indices);
die unless($#indices == $#commands);

for(my $iI = 0; $iI <= $#indices; $iI++)
{
	my $i = $indices[$iI];
	# if(($runningCommands >= 10) or ($runningSize >= 1e6))
	if(($runningCommands >= 10))
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
	foreach my $iI (@indices)
	{
		my $command = $commands[$iI];
		while(scalar(keys %still_running) >= $number_processes)
		{
			sleep 10;
		}
		
		my $pid = fork;
		die "fork failed" unless defined $pid;
		if ($pid == 0) {
			system($command) and die "Could not execute command: $command";
			exit;
		}		
		$still_running{$pid} = 1;			
	}
	
	# print "\nTotal command batches: ", scalar(@command_batches), "\n\n";

	# foreach my $command_batch (@command_batches)
	# {		
		# print "\nNow starting: ", scalar(@$command_batch), " commands.\n";
		# foreach my $cmd (@$command_batch)
		# {
			# my $pid = fork;
			# die "fork failed" unless defined $pid;
			# if ($pid == 0) {
				# system($cmd) and die "Could not execute command: $cmd";
				# exit;
			# }		
			# $still_running{$pid} = 1;
		# }
		
				
		# while( scalar(keys %still_running) > 0 ) {
			# print "Waiting... PIDs " . join(", ", keys %still_running), " still running.\n";		
			# sleep 10;
		# }		
	# }

	print "\n\nProcesses launched.\n";
}

sub reap_kids {
    local $!; # good practice. avoids changing errno.
    while (1) {
        my $kid = waitpid( -1, POSIX->WNOHANG );
        last unless ($kid > 0); # No more to reap.
        delete $still_running{$kid}; # untrack kid.
    }
}
