#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

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

## Usage:
## CALLMAFFT.pl --action <'kickOff', 'check', 'reprocess', 'processChunk'> 
##              --mafftDirectory <path to subdirectory forMAFFT>
##              --qsub <only used with action 'kickOff', default set to 1>
##              --chunkI <only used with action 'processChunk', JobID>
##
## Example command
## ./CALLMAFFT.pl --action kickOff --mafftDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT --qsub 1
## ./CALLMAFFT.pl --action check --mafftDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT
## ./CALLMAFFT.pl --action reprocess --mafftDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT
## ./CALLMAFFT.pl --action processChunk --mafftDirectory /data/projects/phillippy/projects/hackathon/intermediate_files/forMAFFT --chunkI 0

my $mafft_bin = find_present_alternative('/data/projects/phillippy/projects/rDNA/mafft/mafft-7.273-with-extensions/install/bin/mafft');
my $FAS2BAM_bin = find_present_alternative('./fas2bam.pl');
my $GRCh38_header_file = find_present_alternative('../config/windowbam.header.txt');

my $mafftDirectory;
my $action;
my $inputTruncatedReads;
my $readsFasta;
my $paranoid = 1;
my $chunkSize = 5;
my $processChunk;
my $chunkI;
my $qsub = 1;
my $path_to_script = $FindBin::Bin.'/'.$FindBin::Script;
my $temp_qsub = 'temp_qsub';
my $reprocess;

GetOptions (
	'action:s' => \$action,
	'mafftDirectory:s' => \$mafftDirectory, 
	'chunkSize:s' => \$chunkSize, 
	'processChunk:s' => \$processChunk, 
	'chunkI:s' => \$chunkI, 
	'qsub:s' => \$qsub, 
	'reprocess:s' => \$reprocess, 
);

unless($action)
{
	die "Please specify --action";
}
my $fn_files_to_process = $mafftDirectory . '/files_to_process';
my $fn_files_to_reprocess = $mafftDirectory . '/files_to_reprocess';
if($action eq 'reprocess')
{
	print "Reprocess....\n\n";
	
	my $dirHandle;
	opendir ($dirHandle, $mafftDirectory) or die "$0: opendir: $!";
	my @sub_directories = grep {-d "$mafftDirectory/$_" && ! /^\.{1,2}$/} readdir($dirHandle);
	closedir $dirHandle;
	
	print "Identified ", scalar(@sub_directories), " directories.\n";
	print join("\n", map {' - '.$_} @sub_directories), "\n";
	
	open(REPROCESS, '>', $fn_files_to_reprocess) or die "Cannot open $fn_files_to_reprocess";
	
	my $n_files = 0;
	my $n_files_found = 0;
	my $n_files_found_withMFA = 0;
	my $n_files_found_withBAM = 0;
	foreach my $sub_directory (@sub_directories)
	{
		my $subDir_fullPath = $mafftDirectory . '/' . $sub_directory;
		my $subDirHandle;
		opendir ($subDirHandle, $subDir_fullPath) or die "Cannot open '$subDir_fullPath'";
	
		my @fasta_files = sort grep {(-f "$subDir_fullPath/$_") &&  ($_ =~ /\.fa$/)} readdir($subDirHandle);
							
		closedir $subDirHandle;
		
		foreach my $fasta_file (@fasta_files)
		{
			$n_files_found++;
			
			my $mfa_file = $subDir_fullPath . '/' . $fasta_file;
			my $bam_file = $subDir_fullPath . '/' . $fasta_file;
			$mfa_file =~ s/\.fa$/.mfa/;
			$bam_file =~ s/\.fa$/.bam/;
		
			if(-e $mfa_file)
			{
				$n_files_found_withMFA++;
			}
			if(-e $bam_file)
			{
				$n_files_found_withBAM++;
			}
			else
			{
				# print "Not existing: $bam_file\n";
			}
			
			unless(-e $bam_file)
			{
				my $fasta_file_full_path = $subDir_fullPath . '/' . $fasta_file;
				$fasta_file_full_path = File::Spec->rel2abs($fasta_file_full_path);
				die unless(-e $fasta_file_full_path);
				print REPROCESS $fasta_file_full_path, "\n";
				$n_files++;
			}
		}
	}
	
	close(REPROCESS);
	
	print "Total found files: $n_files_found\n";
	print "With MFA: $n_files_found_withMFA \n";
	print "With BAM: $n_files_found_withBAM \n\n";
	print "Now redo: $n_files \n";
	
	my $n_chunks = ceil($n_files / $chunkSize);
	
	invoke_self_array($n_chunks-1, 1);
	
}
elsif($action eq 'check')
{
	print "Check....\n\n";
	
	my $dirHandle;
	opendir ($dirHandle, $mafftDirectory) or die "$0: opendir: $!";
	my @sub_directories = grep {-d "$mafftDirectory/$_" && ! /^\.{1,2}$/} readdir($dirHandle);
	closedir $dirHandle;
	
	print "Identified ", scalar(@sub_directories), " directories.\n";
	print join("\n", map {' - '.$_} @sub_directories), "\n";
		
	my $n_files = 0;
	my $n_files_found = 0;
	my $n_files_found_withMFA = 0;
	my $n_files_found_withBAM = 0;
	foreach my $sub_directory (@sub_directories)
	{
		my $subDir_fullPath = $mafftDirectory . '/' . $sub_directory;
		my $subDirHandle;
		opendir ($subDirHandle, $subDir_fullPath) or die "Cannot open '$subDir_fullPath'";
	
		my @fasta_files = sort grep {(-f "$subDir_fullPath/$_") &&  ($_ =~ /\.fa$/)} readdir($subDirHandle);
							
		closedir $subDirHandle;
		
		foreach my $fasta_file (@fasta_files)
		{
			$n_files_found++;
			
			my $mfa_file = $subDir_fullPath . '/' . $fasta_file;
			my $bam_file = $subDir_fullPath . '/' . $fasta_file;
			$mfa_file =~ s/\.fa$/.mfa/;
			$bam_file =~ s/\.fa$/.bam/;
		
			if(-e $mfa_file)
			{
				$n_files_found_withMFA++;
			}
			if(-e $bam_file)
			{
				$n_files_found_withBAM++;
			}
			else
			{
				# print "Not existing: $bam_file\n";
			}
			
			unless(-e $bam_file)
			{
				my $fasta_file_full_path = $subDir_fullPath . '/' . $fasta_file;
				$fasta_file_full_path = File::Spec->rel2abs($fasta_file_full_path);
				die unless(-e $fasta_file_full_path);
				$n_files++;
			}
		}
	}
	
	close(REPROCESS);
	
	print "Total found files: $n_files_found\n";
	print "With MFA: $n_files_found_withMFA \n";
	print "With BAM: $n_files_found_withBAM \n\n";
	print "Would now redo: $n_files \n";
}
elsif($action eq 'kickOff')
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

	my $useFile = ($reprocess) ? $fn_files_to_reprocess : $fn_files_to_process;
	my @files_to_process;
	open(FILES_TO_PROCESS, '<', $useFile) or die "Cannot open $useFile";
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
		# todo
		#unlink($msaFile);
	}
}
else
{
	die "Specified unknown --action $action";
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
	my $reprocess = shift;
	
	my $reprocess_string = ($reprocess) ? " --reprocess 1 " : '';
	my $mafftDirectory_abs = File::Spec->rel2abs($mafftDirectory);
	if($qsub)
	{	
		open(QSUB, '>', $temp_qsub) or die "Cannot open $temp_qsub";
		my $minJobID = 1;
		my $maxJobID = $maxChunk_0based + 1;
			
		print QSUB qq(#!/bin/bash
#\$ -t ${minJobID}-${maxJobID}
#\$ -q public.q
#\$ -tc 100
#\$ -l mem_free=8G
#\$ -N 'CALLMAFFT'
jobID=\$(expr \$SGE_TASK_ID - 1)
cd $current_dir
perl $path_to_script --mafftDirectory $mafftDirectory_abs --action processChunk --chunkI \$jobID --chunkSize $chunkSize $reprocess_string
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
			my $cmd = qq(perl $path_to_script --mafftDirectory $mafftDirectory_abs --action processChunk --chunkI $chunkI --chunkSize $chunkSize $reprocess_string);
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
	my $n_sequences_originallyIn = (scalar keys %$input_href);
	
	my $temp_file_in = $inputFile . ".tmp_in";
	my $temp_file_out = $inputFile . ".tmp_out";
	
	my %empty_sequences;
	foreach my $key (keys %$input_href)
	{
		my $seq = $input_href->{$key};
		die unless(length($seq));
		$seq =~ s/[\-_]//g;
		if(length($seq) == 0)
		{
			$empty_sequences{$key} = 1;
		}
		else
		{
			$input_href->{$key} = $seq;
		}
	}
	
	die if(scalar(keys %$input_href) == scalar(keys %empty_sequences));
	foreach my $emptyKey (keys %empty_sequences)
	{
		die unless(exists $input_href->{$emptyKey});
		delete $input_href->{$emptyKey};
	}	
	
	if(scalar(keys %$input_href) >= 2)
	{
		writeFASTA($temp_file_in, $input_href);
		#my $cmd_mafft = qq($mafft_bin --auto --quiet $temp_file_in > $temp_file_out);
		my $cmd_mafft = qq($mafft_bin --retree 1 --maxiterate 0 --quiet $temp_file_in > $temp_file_out);
		print "Executing $cmd_mafft \n";
		
		my $ret = system($cmd_mafft);

		# print "Return code $ret\n";
		
		unless(($ret == 0) and (-e $temp_file_out))
		{
			die "File $temp_file_out not there, but after $cmd_mafft it should be";
		}
		
		validate_as_alignment($temp_file_out);
	}
	else
	{
		writeFASTA($temp_file_out, $input_href);
		#cp($inputFile, $temp_file_out) or die "Cannot cp $inputFile $temp_file_out";
	}
	
	if(keys %empty_sequences)
	{
		my $temp_out_href = readFASTA($temp_file_out);
		my $l = length((values %$temp_out_href)[0]);
		foreach my $k (keys %empty_sequences)
		{
			my $gaps = ('-' x $l);
			die unless(length($gaps) == $l);
			$temp_out_href->{$k} = $gaps;
		}
		writeFASTA($outputFile, $temp_out_href);
		
	}
	else
	{
		cp($temp_file_out, $outputFile) or die "Cannot cp $inputFile $temp_file_out";	
	}
	
	validate_as_alignment($outputFile);
	my $output_href = readFASTA($outputFile);
	die unless(scalar(keys %$output_href) == $n_sequences_originallyIn);
	
	# todo
	#unlink($temp_file_in);
	#unlink($temp_file_put);
	
	
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
