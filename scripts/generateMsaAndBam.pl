#!/usr/bin/perl

## Author: Alexander Dilthey (HHU/UKD, NHGRI-NIH), Evan Biederstedt (NYGC), Nathan Dunn (LBNL), Nancy Hansen (NIH), Aarti Jajoo (Baylor), Jeff Oliver (Arizona), Andrew Olsen (CSHL), Torsten Houwaart (HHU/UKD)
## License: The MIT License, https://github.com/NCBI-Hackathons/Graph_Genomes/blob/master/LICENSE

use strict;
use Getopt::Long;   
use Data::Dumper;
use File::Path;
use POSIX qw(ceil);
use File::Copy "cp";
use File::Spec;
use MSAandBAM;

$| = 1;

## Usage:
## generateMsaAndBam.pl <'makeMSA', 'makeBAM'> --fas2bam_path f2bp [--samtools_path stp] [--mafft_path mp] [--useGinsi] [--usePreClustering] [--preCluster_k pck] [--preCluster_jaccard_threshold pcjt] inputfile outputfile

my $mafft_path;
my $samtools_path;
my $bamheader;
my $useGinsi=0;
my $usePreClustering = 0;
my $preCluster_k = 21;
my $preCluster_jaccard_threshold = 0.2; # 0.95**21/(2 - 0.95**21)
my $fas2bam_path;

GetOptions (
	'mafft_path:s' => \$mafft_path,
	'useGinsi:s' => \$useGinsi,
	'fas2bam_path:s' => \$fas2bam_path,
	'samtools_path:s' => \$samtools_path,
	'bamheader:s' => \$bamheader,
	'usePreClustering:s' => \$usePreClustering,
	'preCluster_k:s' => \$preCluster_k,
	'preCluster_jaccard_threshold:s' => \$preCluster_jaccard_threshold,
);

my $action = shift @ARGV;
my $inputPath = shift @ARGV;
my $outputPath = shift @ARGV;


unless($action)
{
	die "Please specify either makeBAM or makeMSA";
}
unless($inputPath)
{
    die "No input file given";
}
unless($outputPath)
{
    die "No output file given";
}

# by default executables should be in $PATH
unless($mafft_path)
{
    $mafft_path = 'mafft';
}
unless($samtools_path)
{
    $samtools_path = 'samtools';
}


if($action eq 'makeMSA')
{
	print "Making MSA....\n\n";
    MSAandBAM::makeMSA($inputPath, $outputPath, $mafft_path, $useGinsi, $usePreClustering, $preCluster_k, $preCluster_jaccard_threshold )
	
}
elsif($action eq 'makeBAM')
{
    unless($fas2bam_path)
    {
        die "Specify --fas2bam_path"
    }
    unless($bamheader){
        die "Specify --bamheader";
    }
	print "Making BAM....\n\n";
    MSAandBAM::makeBAM($inputPath, $outputPath, $bamheader, $samtools_path, $fas2bam_path)

}
else
{
	die "Specified action unknown $action. Should be makeMSA or makeBAM";
}
