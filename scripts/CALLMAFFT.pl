#!/usr/bin/perl

use strict;
use Getopt::Long;   
use Data::Dumper;
use File::Path;
use POSIX qw(ceil);
use File::Copy "cp";
use FindBin;
use File::Spec;
use Cwd;
my $current_dir = getcwd;

$| = 1;

# Example command
# ./CALLMAFFT.pl --action kickOff --mafftDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT --qsub 1
# ./CALLMAFFT.pl --action processChunk --mafftDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT --chunkI 0
#samtools sort -o /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.bam; samtools index /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam 
#samtools sort -o /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_2.bam /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_2.bam; samtools index /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_2.bam 
# ./fas2bam.pl --input /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.mfa --output /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.bam --ref "ref" --bamheader ../config/windowbam.header.txt; samtools sort -o /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.bam; samtools index /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam; ./validate_BAM_MSA.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.fa

# ./validate_BAM_MSA.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_1.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_1.fa
# ./validate_BAM_MSA.pl --BAM /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/sorted_chr1_2.bam --referenceFasta /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT/chr1/chr1_2.fa

my $mafft_bin = find_present_alternative('/data/projects/phillippy/projects/rDNA/mafft/mafft-7.273-with-extensions/install/bin/mafft');
my $FAS2BAM_bin = find_present_alternative('./fas2bam.pl');
my $validateMSABAM = find_present_alternative('./validate_BAM_MSA.pl');
my $GRCh38_header_file = find_present_alternative('../config/windowbam.header.txt');

my $mafftDirectory;
my $action;
my $inputTruncatedReads;
my $readsFasta;
my $paranoid = 1;
my $chunkSize = 1000;
my $processChunk;
my $chunkI;
my $qsub = 1;
my $path_to_script = $FindBin::Bin.'/'.$FindBin::Script;
my $temp_qsub = 'temp_qsub';

GetOptions (
	'action:s' => \$action,
	'mafftDirectory:s' => \$mafftDirectory, 
	'chunkSize:s' => \$chunkSize, 
	'processChunk:s' => \$processChunk, 
	'chunkI:s' => \$chunkI, 
	'qsub:s' => \$qsub, 
);

unless($action)
{
	die "Please specify --action";
}
my $fn_files_to_process = $mafftDirectory . '/files_to_process';
if($action eq 'kickOff')
{
	die unless($chunkSize);
	
	unless($mafftDirectory)
	{
		die "Please specify --mafftDirectory";
	}
	unless((-e $mafftDirectory) and (-e $mafftDirectory))
	{
		die "Please specify valid directory as --mafftDirectory -- $mafftDirectory is not";
	}
	
	open(FILES_TO_PROCESS, '>', $fn_files_to_process) or die "Cannot open $fn_files_to_process";
	
	my $dirHandle;
	opendir ($dirHandle, $mafftDirectory) or die "$0: opendir: $!";
	
	
	my @sub_directories = grep {-d "$mafftDirectory/$_" && ! /^\.{1,2}$/} readdir($dirHandle);
	
	closedir $dirHandle;
	
	print "Identified ", scalar(@sub_directories), " directories.\n";
	print join("\n", map {' - '.$_} @sub_directories), "\n";
	
	my $n_files = 0;
	foreach my $sub_directory (@sub_directories)
	{
		my $subDir_fullPath = $mafftDirectory . '/' . $sub_directory;
		my $subDirHandle;
		opendir ($subDirHandle, $subDir_fullPath) or die "Cannot open '$subDir_fullPath'";
	
		my @fasta_files = sort grep {(-f "$subDir_fullPath/$_") &&  ($_ =~ /\.fa$/)} readdir($subDirHandle);
							
		closedir $subDirHandle;
		
		foreach my $fasta_file (@fasta_files)
		{
			my $fasta_file_full_path = $subDir_fullPath . '/' . $fasta_file;
			$fasta_file_full_path = File::Spec->rel2abs($fasta_file_full_path);
			die unless(-e $fasta_file_full_path);
			print FILES_TO_PROCESS $fasta_file_full_path, "\n";
			$n_files++;
			# last if($n_files > 100);
		}
	}
	
	my $n_chunks = ceil($n_files / $chunkSize);
	close(FILES_TO_PROCESS);
	
	invoke_self_array($n_chunks-1);
}
elsif($action eq 'processChunk')
{
	die unless(defined $chunkI);
	die unless($chunkSize);
	my $firstLine = $chunkI * $chunkSize;
	my $lastLine = ($chunkI + 1) * $chunkSize - 1;
	print "Go from line $firstLine to $lastLine\n";

	my @files_to_process;
	open(FILES_TO_PROCESS, '<', $fn_files_to_process) or die "Cannot open $fn_files_to_process";
	my $lineI = -1;
	while(<FILES_TO_PROCESS>)
	{
		$lineI++;
		if(($lineI >= $firstLine) and ($lineI <= $lastLine))
		{
			my $file = $_;
			chomp($file);
			next unless($file);
			die "File $file not existing" unless(-e $file);
			push(@files_to_process, $file);
		}
	}
	close(FILES_TO_PROCESS);
	
	foreach my $file (@files_to_process)
	{
		print "Processing $file \n";
		die "File weird name: $file" unless($file=~ /\.fa$/);
		my $msaFile = $file;
		$msaFile=~ s/\.fa$/.mfa/;
		
		my $bamFile = $file;
		$bamFile=~ s/\.fa$/.bam/;
		
		makeMSA($file, $msaFile);
		makeBAM($msaFile, $bamFile);
		checkBAM($bamFile, $file);
		unlink($msaFile);
	}
}
else
{
	die "Specified unknown --action $action";
}

sub checkBAM
{
	my $BAM = shift;
	my $fasta = shift;
	
	my $BAM_sorted = $BAM . '.sorted.bam';
	my $BAM_sorted_bai = $BAM_sorted . '.bai';
	
	my $cmd_sort = qq(samtools sort -o $BAM_sorted $BAM);
	if(system($cmd_sort))
	{
		die "Command $cmd_sort failed";
	}
	
	my $cmd_index = qq(samtools index $BAM_sorted);
	if(system($cmd_index))
	{
		die "Command $cmd_index failed";
	}	
	
	unless(-e $BAM_sorted_bai)
	{
		die "Expected file $BAM_sorted_bai not present!";
	}
	
	my $command_check = qq($validateMSABAM --BAM $BAM_sorted --referenceFasta $fasta);
	if(system($command_check))
	{
		die "Command $command_check failed";
	}
	
	unlink($BAM_sorted);
	unlink($BAM_sorted_bai);
}

sub makeBAM
{
	my $inputFile = shift;
	my $outputFile = shift;
	die unless($inputFile and $outputFile);
	
	validate_as_alignment($inputFile);
	
	my $cmd_makeBAM = qq($FAS2BAM_bin --input $inputFile --output $outputFile --ref "ref" --bamheader $GRCh38_header_file);
	
	# print $cmd_makeBAM, "\n";
	
	my $attempt = 0;
	my $ret;
	for($attempt = 0; $attempt < 5; $attempt++)
	{
		$ret = system($cmd_makeBAM);
		if($ret == 0)
		{
			if(-e $outputFile)
			{
				last;
			}
		}
		else
		{
			validate_as_alignment($inputFile);
			#if($attempt > 5)
			#{
			#	die "Five attempts at $cmd_makeBAM failed!";
			#}
		}
	}
	
	unless(($ret == 0) and (-e $outputFile))
	{
		die "File $outputFile not there, but after $cmd_makeBAM it should be - attempts $attempt - last exit status $ret!";
	}
}
sub invoke_self_array
{
	my $maxChunk_0based = shift;
	
	my $mafftDirectory_abs = File::Spec->rel2abs($mafftDirectory);
	if($qsub)
	{	
		open(QSUB, '>', $temp_qsub) or die "Cannot open $temp_qsub";
		my $minJobID = 1;
		my $maxJobID = $maxChunk_0based + 1;
			
		print QSUB qq(#!/bin/bash
#\$ -t ${minJobID}-${maxJobID}
#\$ -l mem_free=8G
#\$ -N 'CALLMAFFT'
jobID=\$(expr \$SGE_TASK_ID - 1)
cd $current_dir
perl $path_to_script --mafftDirectory $mafftDirectory_abs --action processChunk --chunkI \$jobID --chunkSize $chunkSize
);
		close(QSUB);
		my $qsub_cmd = "qsub $temp_qsub";
		if(system($qsub_cmd))
		{
			die "qsub $qsub_cmd failed";
		}
	}
	else
	{
		for(my $chunkI = 0; $chunkI <= $maxChunk_0based; $chunkI++)
		{
			print "Call myself for chunk $chunkI\n";
			my $cmd = qq(perl $path_to_script --mafftDirectory $mafftDirectory_abs --action processChunk --chunkI $chunkI --chunkSize $chunkSize);
			if(system($cmd))
			{
				die "Command $cmd failed";
			}
		}
	}
}


sub makeMSA
{
	my $inputFile = shift;
	my $outputFile = shift;
	die unless($inputFile and $outputFile);
	
	my $input_href = readFASTA($inputFile);
	if(scalar(keys %$input_href) >= 2)
	{

	my $cmd_mafft = qq($mafft_bin --auto $inputFile > $outputFile);
	my $ret = system($cmd_mafft);

	unless(($ret == 0) and (-e $outputFile))
	{
		die "File $outputFile not there, but after $cmd_mafft it should be";
	}
	}
	else
	{
		cp($inputFile, $outputFile) or die "Cannot cp $inputFile $outputFile";
	}
	validate_as_alignment($outputFile);
}

sub validate_as_alignment
{
	my $inputFile = shift;
	my $alignment_href = readFASTA($inputFile);
	
	if(scalar(keys %$alignment_href))
	{
		my $l;
		foreach my $key (keys %$alignment_href)
		{
			my $seq_l = length($alignment_href->{$key});
			if(not defined $l)
			{
				$l = $seq_l;
			}
			unless($seq_l == $l)
			{
				die "File $inputFile is not an alignment, length mismatch - $l vs $seq_l";
			}
		}
	}
	else
	{
		die "File $inputFile does not contain any sequences.";
	}
	
}

sub count_FASTA_sequences
{
	my $file = shift;
	
	my $n_sequences = 0;
	open(F, '<', $file) or die "Cannot open $file";
	while(<F>)
	{
		my $line = $_;
		if($line and(substr($line, 0, 1) eq '>'))
		{
			$n_sequences++;
		}
	}
	close(F);
	
	return $n_sequences;
}

sub find_present_alternative
{
        foreach my $a (@_)
        {
                if(-e $a)
                {
                        return $a;
                }
        }
        die "Could not find a present alternative from list:\n".join("\n", @_);
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
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}
