package MSAandBAM;

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL), Torsten Houwaart (HHU/UKD)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use warnings;
use File::Copy "cp";

sub makeBAM
{
	my $inputFile = shift;
	my $outputFile = shift;
    my $bamheader = shift;
    my $samtools_path = shift;
    my $fas2bam_path = shift;

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

sub makeMSA
{
	my $inputFile = shift;
	my $outputFile = shift;
	my $mafft_path = shift;
    my $useGinsi = shift;
	my $usePreClustering = shift;
    my $preCluster_k = shift;
    my $preCluster_jaccard_threshold = shift;

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
			createMSA_preCluster($temp_file_in, $temp_file_out, $mafft_path, $usePreClustering,  $preCluster_k, $preCluster_jaccard_threshold);
		}
		else
		{
			my $cmd_mafft = qq(${mafft_path} --retree 1 --maxiterate 0 --quiet $temp_file_in > $temp_file_out);
			my $cmd_mafft_ginsi = qq(${mafft_path}-ginsi --quiet $temp_file_in > $temp_file_out);
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

sub lastOption_linearMSA
{
	my $temp_file_in = shift;
	my $temp_file_out = shift;
	
	my $sequences_href = readFASTA($temp_file_in);
	die unless(scalar(keys %$sequences_href));

	my %combined_sequences;
	my $joint_length = 0;
	foreach my $seqID (sort keys %$sequences_href)
	{
		my $seqL = length($sequences_href->{$seqID});
		$joint_length += $seqL;
		foreach my $seqID2 (sort keys %$sequences_href)
		{
			if($seqID2 eq $seqID)
			{
				my $seq_to_add = uc($sequences_href->{$seqID2});
				die unless(length($seq_to_add) == $seqL);
				$combined_sequences{$seqID2} .= $seq_to_add;
			}
			else
			{
				my $seq_to_add = ('-' x $seqL);
				die unless(length($seq_to_add) == $seqL);
				$combined_sequences{$seqID2} .= $seq_to_add;				
			}
		}
	}
	
	# print "Joint sequence length $temp_file_out : $joint_length \n";
	
	foreach my $k (keys %combined_sequences)
	{
		die unless(length($combined_sequences{$k}) == $joint_length);
	}
	
	writeFASTA($temp_file_out, \%combined_sequences);
}

sub createMSA_preCluster
{
	my $temp_file_in = shift;
	my $temp_file_out = shift;
    my $mafft_path = shift;
	my $usePreClustering = shift;
	my $k = shift;
	my $jaccardT = shift;
	die unless(defined $mafft_path);
	
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
			my $cmd_mafft = qq(${mafft_path} --retree 1 --maxiterate 0 --quiet $fn_rawSeq > $fn_msaSeq);
			my $mafft_retcode;
			if($usePreClustering ne 'ultra')
			{
				$mafft_retcode = system($cmd_mafft);
			}
			if(($usePreClustering eq 'ultra') or ($mafft_retcode != 0))
			{
				warn "Command $cmd_mafft failed with return code $mafft_retcode or --usePreClustering set to 'ultra' - will linearly concatenate the input sequences";
				lastOption_linearMSA($fn_rawSeq, $fn_msaSeq);
			}
			else
			{
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

1;
