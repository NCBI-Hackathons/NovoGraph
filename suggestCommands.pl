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
## suggestCommands.pl --referenceGenome <path to reference in FASTA format> --inputContigs <path to input contigs in FASTA format>
## Example:
## suggestCommands.pl --referenceGenome Jason/pgf_grch37.fa --inputContigs Jason/JasonLongContigs_0525.fa --outputDirectory Jason --prefix 0525 --mafft_executable /software/mafft/7.407/login/bin/mafft --samtools_path /software/samtools/1.6/ivybridge/bin/samtools --limitToPGF 1

my $referenceGenome;
my $inputContigs;
my $qsub = 1;
my $outputDirectory;
my $prefix;
my $mafft_executable;
my $samtools_path;
my $limitToPGF = 0;
GetOptions (
	'referenceGenome:s' => \$referenceGenome,
	'inputContigs:s' => \$inputContigs, 
	'outputDirectory:s' => \$outputDirectory, 
	'prefix:s' => \$prefix, 
	'qsub:s' => \$qsub, 
	'mafft_executable:s' => \$mafft_executable, 
	'samtools_path:s' => \$samtools_path, 
	'limitToPGF:s' => \$limitToPGF, 
);

unless(-e 'src/CRAM2VCF')
{
	warn "Could not find src/CRAM2VCF, either you don't call me from the NovoGraph root, or you have not compiled the C++ component yet.";
}

unless($mafft_executable)
{
	my $which_output = `which mafft`;
	$which_output =~ s/[\n\r]//g;
	if($which_output and (-e $which_output))
	{
		$mafft_executable = $which_output;
	}		
	else
	{
		die "Could not determine path to MAFFT by 'which mafft', please specify via --mafft_executable";
	}
}

unless($samtools_path)
{
	my $which_output = `which samtools`;
	$which_output =~ s/[\n\r]//g;
	if($which_output and (-e $which_output))
	{
		$samtools_path = $which_output;
	}		
	else
	{
		die "Could not determine path to samtools by 'which samtools', please specify via --samtools_path";
	}
}

unless(defined $outputDirectory)
{
	die "Please specify --outputDirectory";
}

unless(defined $referenceGenome)
{
	die "Please specify --reference";
}

unless(-e $referenceGenome)
{
	warn "File specified via --referenceGenome not existing";
}

unless(defined $inputContigs)
{
	die "Please specify --inputContigs";
}

unless(-e $inputContigs)
{
	warn "File specified via --inputContigs not existing";
}

unless(-d $outputDirectory)
{
	print "mkdir $outputDirectory\n";
	print "mkdir ${outputDirectory}/intermediate_files\n";
	
}

my $ginsiMAFFT = ($limitToPGF) ? ' --useGinsi 1 ' : '';

print qq(
bwa mem -t 4 $referenceGenome  $inputContigs | samtools view -Sb - > ${outputDirectory}/${prefix}_allContigs_unsorted.bam;\
$samtools_path sort -o  ${outputDirectory}/${prefix}_allContigs_sorted.bam  ${outputDirectory}/${prefix}_allContigs_unsorted.bam;\
$samtools_path index  ${outputDirectory}/${prefix}_allContigs_sorted.bam;\
$samtools_path view -c -f 0x4  ${outputDirectory}/${prefix}_allContigs_sorted.bam

# Make sure the last samtools view command returned 0!

perl scripts/checkBAM_SVs_and_INDELs.pl --BAM  ${outputDirectory}/${prefix}_allContigs_sorted.bam --referenceFasta $referenceGenome --readsFasta $inputContigs;\
perl scripts/BAM2ALIGNMENT.pl --BAM  ${outputDirectory}/${prefix}_allContigs_sorted.bam --referenceFasta $referenceGenome --readsFasta $inputContigs --outputFile  ${outputDirectory}/intermediate_files/${prefix}_AlignmentInput.txt;\
perl scripts/FIND_GLOBAL_ALIGNMENTS.pl --alignmentsFile  ${outputDirectory}/intermediate_files/${prefix}_AlignmentInput.txt.sortedWithHeader  --referenceFasta $referenceGenome --outputFile  ${outputDirectory}/${prefix}_forMAFFT.bam --outputTruncatedReads  ${outputDirectory}/${prefix}_truncatedReads --outputReadLengths  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --CIGARscript_path scripts/dealWithTooManyCIGAROperations.pl;\
perl scripts/countExpectedGlobalAlignments.pl --BAM  ${outputDirectory}/${prefix}_forMAFFT.bam

perl scripts/BAM2MAFFT.pl --BAM  ${outputDirectory}/${prefix}_forMAFFT.bam --referenceFasta $referenceGenome --readsFasta $inputContigs --outputDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --inputTruncatedReads  ${outputDirectory}/${prefix}_truncatedReads  --processNonChrReferenceContigs 1;\
perl scripts/CALLMAFFT.pl --action kickOff --mafftDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --mafft_executable $mafft_executable --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt --qsub $qsub $ginsiMAFFT
perl scripts/CALLMAFFT.pl --action reprocess --mafftDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --mafft_executable $mafft_executable --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt --qsub 0

perl scripts/globalize_windowbams.pl --fastadir  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT/  --msadir  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT/ --contigs  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --output  ${outputDirectory}/${prefix}_combined.sam

$samtools_path view -h -t pgf.GrCh38.headerfile.txt  ${outputDirectory}/${prefix}_combined.sam >  ${outputDirectory}/${prefix}_combined.sam_with_header.sam;\
$samtools_path sort  ${outputDirectory}/${prefix}_combined.sam_with_header.sam -o  ${outputDirectory}/${prefix}_combined.sam_with_header_sorted.sam;\
cat  ${outputDirectory}/${prefix}_combined.sam_with_header_sorted.sam | $samtools_path view -C -T $referenceGenome - >  ${outputDirectory}/${prefix}_combined.cram;\
$samtools_path index  ${outputDirectory}/${prefix}_combined.cram;\
perl scripts/checkMAFFT_input_and_output.pl --MAFFTdir  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT/  --contigLengths  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --preMAFFTBAM  ${outputDirectory}/${prefix}_forMAFFT.bam  --finalOutputCRAM  ${outputDirectory}/${prefix}_combined.cram --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt

perl scripts/CRAM2VCF.pl --CRAM  ${outputDirectory}/${prefix}_combined.cram  --referenceFasta $referenceGenome --prefix ${outputDirectory}/${prefix}_finalVCF  --contigLengths  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --CRAM2VCF_executable src/CRAM2VCF
);

if($limitToPGF)
{
	print qq(
# Generate VCF specifically for PGF, GIAB style:
src/CRAM2VCF --input ${outputDirectory}/${prefix}_finalVCF.part_pgf --referenceSequenceID pgf --max_gap_length 1000000 --doProduce_pseudoSample 1 --doProduce_separatedVCF 0	
);
}
else
{
print qq(
# Start VCF generation process for all reference contigs:
perl scripts/launch_CRAM2VCF_C++.pl --prefix ${outputDirectory}/${prefix}_finalVCF
);	
}




