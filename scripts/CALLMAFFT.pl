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
use List::Util qw/shuffle/;
use Cwd;
my $current_dir = getcwd;

$| = 1;

## Usage:
## CALLMAFFT.pl --action <'kickOff', 'check', 'reprocess', 'processChunk'> 
##              --mafftDirectory <path to subdirectory forMAFFT>
##              --qsub <only used with action 'kickOff', default set to 1>
##              --chunkI <only used with action 'processChunk', JobID>
##              --mafft_executable <path to MAFFT executable in /bin, required>
##              --fas2bam_path <path to script 'fas2bam.pl', required>
##              --samtools_path <path to SAMtools executable for fas2bam.pl, required>
##              --bamheader <path to file containing header for BAM file for fas2bam.pl, required>
##
## Example command
## ./CALLMAFFT.pl --action kickOff --mafftDirectory ../intermediate_files/forMAFFT --qsub 1 
##                --mafft_executable /mafft/mafft-7.273-with-extensions/install/bin/mafft --fas2bam_path /intermediate_files/fas2bam.pl --samtools_path /usr/local/bin/samtools --bamheader windowbam.header.txt
## ./CALLMAFFT.pl --action check --mafftDirectory ../intermediate_files/forMAFFT
##                --mafft_executable /mafft/mafft-7.273-with-extensions/install/bin/mafft --fas2bam_path /intermediate_files/fas2bam.pl --samtools_path /usr/local/bin/samtools --bamheader windowbam.header.txt
## ./CALLMAFFT.pl --action reprocess --mafftDirectory ../intermediate_files/forMAFFT
##                --mafft_executable /mafft/mafft-7.273-with-extensions/install/bin/mafft --fas2bam_path /intermediate_files/fas2bam.pl --samtools_path /usr/local/bin/samtools --bamheader windowbam.header.txt
## ./CALLMAFFT.pl --action processChunk --mafftDirectory ../intermediate_files/forMAFFT --chunkI 0
##                --mafft_executable /mafft/mafft-7.273-with-extensions/install/bin/mafft --fas2bam_path /intermediate_files/fas2bam.pl --samtools_path /usr/local/bin/samtools --bamheader windowbam.header.txt
##

## Note: For the majority of use cases, this file 'windowbam.header.txt' should remain untouched. 
## Users should use the file 'windowbam.header.txt' as provided, unless there are assemblies with contigs 
## longer than chr1 in hg38, 248956422 bp. In this case, please change this value to be the size of the largest contig. 

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
my $mafft_executable;
my $fas2bam_path;
my $samtools_path;
my $bamheader;
my $PBSPro = 0;
my $PBSPro_select;
my $PBSPro_A;
my $torque = 0;
my $torque_select;
my $torque_A;
my $preExec;
my $useGinsi;
my $usePreClustering;
my $preCluster_k = 21;
my $preCluster_jaccard_threshold = 0.2; # 0.95**21/(2 - 0.95**21)

GetOptions (
	'action:s' => \$action,
	'mafftDirectory:s' => \$mafftDirectory, 
	'chunkSize:s' => \$chunkSize, 
	'processChunk:s' => \$processChunk, 
	'chunkI:s' => \$chunkI, 
	'qsub:s' => \$qsub, 
	'reprocess:s' => \$reprocess, 
	'mafft_executable:s' => \$mafft_executable,
	'fas2bam_path:s' => \$fas2bam_path,
	'samtools_path:s' => \$samtools_path,
	'bamheader:s' => \$bamheader,
	'PBSPro:s' => \$PBSPro,
	'PBSPro_select:s' => \$PBSPro_select,
	'PBSPro_A:s' => \$PBSPro_A,
	'torque:s' => \$torque,
	'torque_select:s' => \$torque_select,
	'torque_A:s' => \$torque_A,
	'useGinsi:s' => \$useGinsi,
	'preExec:s' => \$preExec,
	'usePreClustering:s' => \$usePreClustering,
	'preCluster_k:s' => \$preCluster_k,
	'preCluster_jaccard_threshold:s' => \$preCluster_jaccard_threshold,
);

# die jaccardSimilary('ACGTAAAT', 'ACGTAAAAAAT', 4);
# my %sequences = (
 # A => 'ACGT',
 # B => 'ACGTAATCGAGA',
 # C => 'ACGTATTCGAGA',
 # D => 'ACTACGAGCGCAGCGACCGAGCGACACTA',
 # E => 'ACTACGTGCGCAGCGACCGAGCGACACTA',
# );
# my $t_n = '_t';
# writeFASTA($t_n, \%sequences);
# createMSA_preCluster($t_n, 'out.txt', 4, 0.3);


die unless($mafft_executable);
die unless($fas2bam_path);
die unless($samtools_path);
die unless($bamheader);
die unless($preCluster_k =~ /^\d+$/);
die unless(($preCluster_jaccard_threshold >= 0) and ($preCluster_jaccard_threshold <= 1));

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
	
	my @command_lines_for_print;
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
			die "File $fasta_file_full_path not existing" unless(-e $fasta_file_full_path);
			push(@command_lines_for_print, $fasta_file_full_path);
			$n_files++;
			# last if($n_files > 100);
		}
	}
	
	my $n_chunks = ceil($n_files / $chunkSize);
	
	@command_lines_for_print = shuffle(@command_lines_for_print);
	print FILES_TO_PROCESS join("\n", @command_lines_for_print), "\n";
	
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
	
	my $cmd_makeBAM = qq(perl $fas2bam_path --input $inputFile --output $outputFile --ref "ref" --bamheader $bamheader --samtools_path $samtools_path);
	
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
	my $ginsi = (defined $useGinsi) ? " --useGinsi $useGinsi " : '';	
	my $preClustering_string = (defined $usePreClustering) ? " --usePreClustering $usePreClustering " : '';	
	
	my $mafftDirectory_abs = File::Spec->rel2abs($mafftDirectory);
	if($qsub)
	{	
		open(QSUB, '>', $temp_qsub) or die "Cannot open $temp_qsub";
		my $minJobID = 1;
		my $maxJobID = $maxChunk_0based + 1;
		if($PBSPro)
		{
		print QSUB qq(#!/bin/bash
#PBS -l $PBSPro_select
#PBS -l walltime=23:00:00
#PBS -A '$PBSPro_A'
#PBS -N CALLMAFFT
#PBS -J ${minJobID}-${maxJobID}
#PBS -r y

jobID=\$(expr \$PBS_ARRAY_INDEX - 1)
);		
		}
		elsif($torque)
		{
		print QSUB qq(#!/bin/bash
#PBS -l $torque_select
#PBS -l walltime=23:00:00
#PBS -A '$torque_A'
#PBS -N CALLMAFFT
#PBS -t ${minJobID}-${maxJobID}
#PBS -r y

jobID=\$(expr \$PBS_ARRAYID - 1)
);
		}
		else
		{
		print QSUB qq(#!/bin/bash
#\$ -t ${minJobID}-${maxJobID}
#\$ -q public.q
#\$ -tc 100
#\$ -l mem_free=8G
#\$ -N 'CALLMAFFT'

jobID=\$(expr \$SGE_TASK_ID - 1)
);
		}
		
		print QSUB qq(
$preExec
cd $current_dir
perl $path_to_script --mafftDirectory $mafftDirectory_abs --action processChunk --chunkI \$jobID --chunkSize $chunkSize --mafft_executable $mafft_executable --fas2bam_path $fas2bam_path --samtools_path $samtools_path --bamheader $bamheader $reprocess_string $ginsi $preClustering_string
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
			my $cmd = qq(perl $path_to_script --mafftDirectory $mafftDirectory_abs --action processChunk --chunkI $chunkI --chunkSize $chunkSize $reprocess_string --mafft_executable $mafft_executable --fas2bam_path $fas2bam_path --samtools_path $samtools_path --bamheader $bamheader $reprocess_string $ginsi $preClustering_string);
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
		
		if($usePreClustering)
		{
			die unless(defined $preCluster_k);
			die unless(defined $preCluster_jaccard_threshold);
			createMSA_preCluster($temp_file_in, $temp_file_out, $preCluster_k, $preCluster_jaccard_threshold);
		}
		else
		{
			my $cmd_mafft = qq(${mafft_executable} --retree 1 --maxiterate 0 --quiet $temp_file_in > $temp_file_out);
			my $cmd_mafft_ginsi = qq(${mafft_executable}-ginsi --quiet $temp_file_in > $temp_file_out);
			my $cmd_mafft_use = $useGinsi ? $cmd_mafft_ginsi : $cmd_mafft;
			
			print "Executing $cmd_mafft_use \n";
			
			my $ret = system($cmd_mafft_use);

			# print "Return code $ret\n";
			
			unless($ret == 0)
			{
				warn "Original mafft command failed, retry with fast settings";
				$ret = system($cmd_mafft);
			}
			
			unless(($ret == 0) and (-e $temp_file_out))
			{
				die "File $temp_file_out not there, but after $cmd_mafft it should be";
			}
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

sub createMSA_preCluster
{
	my $temp_file_in = shift;
	my $temp_file_out = shift;
	my $k = shift;
	my $jaccardT = shift;
	die unless(defined $mafft_executable);
	
	die unless(defined $temp_file_in);
	die unless(defined $temp_file_out);
	die unless(defined $k);
	die unless(defined $jaccardT);

	my $sequences_href = readFASTA($temp_file_in);
	die unless(scalar(keys %$sequences_href));
	
	my @sequences_keys_tooShort = grep {length($sequences_href->{$_}) <= $k} sort keys %$sequences_href;
	my @sequences_keys_longEnough = grep {length($sequences_href->{$_}) > $k} sort keys %$sequences_href;

	my %sequence_to_cluster = map {$_ => 0} @sequences_keys_tooShort;
	my %cluster_to_sequence = (0 => {map {$_ => 1} @sequences_keys_tooShort});
		
	my $verbose = 0;
	
	for(my $i1 = 0; $i1 <= $#sequences_keys_longEnough; $i1++)
	{
		my $k1 = $sequences_keys_longEnough[$i1];	
		my $s1 = $sequences_href->{$k1};
		
		my $existingCluster;
		die if(defined $sequence_to_cluster{$k1});
		
		print "Analyzing $k1...\n" if($verbose);	
		
		for(my $i2 = 0; $i2 <= $#sequences_keys_longEnough; $i2++)
		{
			my $k2 = $sequences_keys_longEnough[$i2];
			
			my $s2 = $sequences_href->{$k2};
			my $jS = jaccardSimilary($s1, $s2, $k);
			
			print "\tComparing $k1 and $k2, lengths " . length($s1) . " and " . length($s2) . "\n" if($verbose);
			if(defined $sequence_to_cluster{$k2})
			{
				print "\t\t$k2 has cluster and Jaccard $jS\n" if($verbose);
			
				if($jS >= $jaccardT)
				{
					print "\t\t\t$k1 and $k2 similar\n" if($verbose);
					$existingCluster = $sequence_to_cluster{$k2};
					last;
				}
				else
				{
					print "\t\t\t$k1 and $k2 not similar\n" if($verbose);
				}
			}
			else
			{
				print "\t\t$k2 has no cluster, ignore.\n" if($verbose);
			}
		}
		
		if(defined $existingCluster)
		{
			$sequence_to_cluster{$k1} = $existingCluster;		
			$cluster_to_sequence{$existingCluster}{$k1} = 1;	
			print "\tAdding $k1 to cluster $existingCluster\n" if($verbose);		
		}
		else
		{
			print "\tNo suitable existing clusters for $k1 found, create new cluster\n" if($verbose);
			my $newClusterID = scalar(keys %sequence_to_cluster) + 1;
			$sequence_to_cluster{$k1} = $newClusterID;
			$cluster_to_sequence{$newClusterID}{$k1} = 1;			
		}
	}
	
	die unless(scalar(keys %sequence_to_cluster) == scalar(keys %$sequences_href));
	
	foreach my $sequenceID (keys %sequence_to_cluster)
	{
		my $clusterID = $sequence_to_cluster{$sequenceID};
		die unless($cluster_to_sequence{$clusterID}{$sequenceID});
	}
	
	print "Identified ", scalar(keys %cluster_to_sequence), " clusters \n" if($verbose);
	
	my @concat_files;
	foreach my $clusterID (sort keys %cluster_to_sequence)
	{
		next unless(scalar(keys %{$cluster_to_sequence{$clusterID}}));
		
		my $fn_rawSeq = $temp_file_out . '.c.' . $clusterID;
		my $fn_msaSeq = $temp_file_out . '.c.' . $clusterID . '.mfa';
		my %cluster_href;
		foreach my $clusterSeqID (keys %{$cluster_to_sequence{$clusterID}})
		{
			die unless(defined $sequences_href->{$clusterSeqID});
			$cluster_href{$clusterSeqID} = $sequences_href->{$clusterSeqID};
			unless($sequences_href->{$clusterSeqID} =~ /^[ACGTN]+$/i)
			{
				warn "Sequence $clusterSeqID in file $temp_file_in is not [ACGTN]+ - mafft may fail because of this.";
			}
		}
		die unless(scalar(keys %cluster_href));		
		writeFASTA($fn_rawSeq, \%cluster_href);
		
		if(scalar(keys %cluster_href) == 1)
		{
			cp($fn_rawSeq, $fn_msaSeq) or die "Copy from $fn_rawSeq to $fn_msaSeq failed";
		}
		else
		{
			my $cmd_mafft = qq(${mafft_executable} --retree 1 --maxiterate 0 --quiet $fn_rawSeq > $fn_msaSeq);
			my $mafft_retcode = system($cmd_mafft);
			if($mafft_retcode != 0)
			{
				die "Command $cmd_mafft failed with return code $mafft_retcode";
			}
		}
		
		validate_as_alignment($fn_msaSeq);	

		push(@concat_files, $fn_msaSeq);
	}
	
	my %joint_sequence_IDs;
	foreach my $f (@concat_files)
	{
		my $sequences_href = readFASTA($f);	
		foreach my $seqID (keys %$sequences_href)
		{
			$joint_sequence_IDs{$seqID}++;
		}
	}
	
	my %combined_sequences;
	foreach my $f (@concat_files)
	{
		my $sequences_href = readFASTA($f);	
		my $seqL = length((values %$sequences_href)[0]);
		foreach my $seqID (keys %joint_sequence_IDs)
		{
			if(exists $sequences_href->{$seqID})
			{
				my $seq_to_add = uc($sequences_href->{$seqID});
				die unless(length($seq_to_add) == $seqL);
				$combined_sequences{$seqID} .= $seq_to_add;
			}
			else
			{
				my $seq_to_add = ('-' x $seqL);
				die unless(length($seq_to_add) == $seqL);
				$combined_sequences{$seqID} .= $seq_to_add;				
			}
		}
	}
	
	writeFASTA($temp_file_out, \%combined_sequences);
	
	validate_as_alignment($temp_file_out);	
	
	print "Produced $temp_file_out\n" if($verbose);
	
	#die Dumper(\%sequence_to_cluster, \%cluster_to_sequence);
	
	#die Dumper(\@sequences_keys_tooShort, \@sequences_keys_longEnough);	
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
	foreach my $key (sort keys %$href)
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

sub jaccardSimilary
{
	my $s1 = shift;
	my $s2 = shift;
	my $k = shift;
	
	die unless(defined $s1);
	die unless(defined $s2);
	die unless(defined $k);
	
	my @kmers_s1 = kmers(uc($s1), $k);
	my @kmers_s2 = kmers(uc($s2), $k);
	
	if((scalar(@kmers_s1) == 0) and (scalar(@kmers_s2) == 0))
	{
		return 0;
	}
	else
	{
		my %_k1 = map {$_ => 1} @kmers_s1;
		my %_k2 = map {$_ => 1} @kmers_s2;
		my %k_union = map {$_ => 1} (@kmers_s1, @kmers_s2);
		my %k_intersection = map {$_ => 1} grep {$_k2{$_}} @kmers_s1;
		die if(scalar(keys %k_union) == 0);
		return scalar(keys %k_intersection)/scalar(keys %k_union);
	}
}	

sub kmers
{
	my $s = shift;
	my $k = shift;
	
	if(length($s) < $k)
	{
		return ();
	}
	my @forReturn;
	my $expectedMers = length($s) - $k + 1;
	for(my $i = 0; $i < $expectedMers; $i++)
	{	
		my $kMer = substr($s, $i, $k);
		die unless(length($kMer) == $k);
		push(@forReturn, $kMer);
	}
	
	die unless(scalar(@forReturn) == $expectedMers);
	return @forReturn;
}
