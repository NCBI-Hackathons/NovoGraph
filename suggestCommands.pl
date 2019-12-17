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

my $ginsiMAFFT = ($limitToPGF) ? '--useGinsi 1 ' : ' ';

my $samtools_view_header = ($limitToPGF) ? 'pgf.GrCh38.headerfile.txt' : 'GRCh38.headerfile.txt';

# $samtools_path sort -o  ${outputDirectory}/${prefix}_allContigs_sorted.bam  ${outputDirectory}/${prefix}_allContigs_unsorted.bam;\
# $samtools_path index  ${outputDirectory}/${prefix}_allContigs_sorted.bam;\
# $samtools_path view -c -f 0x4  ${outputDirectory}/${prefix}_allContigs_sorted.bam

 
## Removed:
# # Produce some summary statistics:
# # (this will only work if you have the Bio::DB::HTS module installed)
# perl scripts/countExpectedGlobalAlignments.pl --BAM  ${outputDirectory}/${prefix}_forMAFFT.cram


print qq(

# If your contigs are interspersed with long stretches of unknown nucleotides (NNN...NNN) you should split them:
# needs: Python 3
python scripts/split_identify_fasta.py $inputContigs > ${outputDirectory}/split_${inputContigs}

# Map your input contigs with minimap2:
minimap2 -t 4 -a -x asm20 $referenceGenome ${outputDirectory}/split_${inputContigs} | $samtools_path view -Sb - > ${outputDirectory}/${prefix}_allContigs_unsorted.bam

# You could also use bwa:
# bwa mem -t 4 $referenceGenome  ${outputDirectory}/split_${inputContigs} | $samtools_path view -Sb - > ${outputDirectory}/${prefix}_allContigs_unsorted.bam

# Check that the following command returns 0 - otherwise remove entries of unmapped entries from BAM:
$samtools_path view -c -f 0x4 ${outputDirectory}/${prefix}_allContigs_unsorted.bam

# Check that all input data look OK:
perl scripts/checkBAM_SVs_and_INDELs.pl --BAM ${outputDirectory}/${prefix}_allContigs_unsorted.bam --referenceFasta $referenceGenome --readsFasta $inputContigs --sam2alignment_executable src/sam2alignment --samtools_path $samtools_path

# Convert BAM into a simple text format readable by the next step:
perl scripts/BAM2ALIGNMENT.pl --BAM  ${outputDirectory}/${prefix}_allContigs_unsorted.bam --referenceFasta $referenceGenome --readsFasta $inputContigs --outputFile  ${outputDirectory}/intermediate_files/${prefix}_AlignmentInput.txt --sam2alignment_executable src/sam2alignment --samtools_path $samtools_path

# Find globally best alignment for each input contig:
perl scripts/FIND_GLOBAL_ALIGNMENTS.pl --alignmentsFile  ${outputDirectory}/intermediate_files/${prefix}_AlignmentInput.txt.sortedWithHeader  --referenceFasta $referenceGenome --outputFile  ${outputDirectory}/${prefix}_forMAFFT.sam --outputTruncatedReads  ${outputDirectory}/${prefix}_truncatedReads --outputReadLengths  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --samtools_path $samtools_path

# Prepare the multiple sequence alignment step:
perl scripts/BAM2MAFFT.pl --SAM  ${outputDirectory}/${prefix}_forMAFFT.sam --referenceFasta $referenceGenome --readsFasta $inputContigs --outputDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --inputTruncatedReads  ${outputDirectory}/${prefix}_truncatedReads  --processNonChrReferenceContigs 1 --sam2alignment_executable src/sam2alignment --samtools_path $samtools_path

# Kick off multiple sequence alignment generation. This commands assumes that you are using an SGE cluster environment.
# If you have no cluster available, you can specify --qsub 0 to directly execute all required MSA commands, but be prepared for this to take a rather long time.
# If you are using PBSPro instead of SGE, you can add the following arguments (modified for your local environment):
#        --PBSPro 1 --PBSPro_select 'select=1:ncpus=16:mem=48GB' --PBSPro_A IMMGEN --preExec 'module load Perl; module load SamTools; module load Mafft/7.407' --chunkSize 500
# If you are using the open source scheduler TORQUE, you can add the following arguments (modified for your local environment):
#        --torque 1 --torque_select 'select=1:ncpus=16:mem=48GB' --torque_A IMMGEN --preExec 'module load Perl; module load SamTools; module load Mafft/7.407' --chunkSize 500
# The --chunkSize parameter determines how many alignment jobs are assigned to each submitted job on your cluster, i.e. as you increase --chunkSize, the total number
# of submitted jobs is reduced.
perl scripts/CALLMAFFT.pl --action kickOff --mafftDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --mafft_executable $mafft_executable --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt --qsub $qsub $ginsiMAFFT
 
# It often happens that individual alignment jobs fail for idiosyncratic reasons. Therefore, when all jobs have finished, try the following command - if all cluster jobs were
# executed successfully, the command will tell you; otherwise, it will try to create the missing alignments. Also supports --qsub 1 and PBSPro parameters like the preceding command.
perl scripts/CALLMAFFT.pl --action reprocess --mafftDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --mafft_executable $mafft_executable --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt --qsub 0

# Now check that all alignments were computed successfully:
perl scripts/CALLMAFFT.pl --action check --mafftDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --mafft_executable $mafft_executable --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt

# If there are still chunks without valid valignments, try adding --usePreClustering 1.
# ... when --usePreClustering 1 is active, the algorithm will try increasigly aggressive multiple sequence alignment strategies.
perl scripts/CALLMAFFT.pl --action reprocess --usePreClustering 1 --mafftDirectory  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT --mafft_executable $mafft_executable --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt --qsub 0

# Combine the multiple sequence alignments (MSAs) created during the previous step:
perl scripts/globalize_windowbams.pl --fastadir  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT/  --msadir  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT/ --contigs  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --output  ${outputDirectory}/${prefix}_combined.sam --samtools_path $samtools_path

# Convert the combined genome-wide MSAs into CRAM format:
$samtools_path view -h -t ${samtools_view_header}  ${outputDirectory}/${prefix}_combined.sam >  ${outputDirectory}/${prefix}_combined.sam_with_header.sam;\
$samtools_path sort  ${outputDirectory}/${prefix}_combined.sam_with_header.sam -o  ${outputDirectory}/${prefix}_combined.sam_with_header_sorted.sam;\
cat  ${outputDirectory}/${prefix}_combined.sam_with_header_sorted.sam | $samtools_path view -C -T $referenceGenome - >  ${outputDirectory}/${prefix}_combined.cram;\
$samtools_path index  ${outputDirectory}/${prefix}_combined.cram

# Check that the data still look OK:
perl scripts/checkMAFFT_input_and_output.pl --MAFFTdir  ${outputDirectory}/intermediate_files/${prefix}_forMAFFT/  --contigLengths  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --preMAFFTBAM  ${outputDirectory}/${prefix}_forMAFFT.sam  --finalOutputCRAM  ${outputDirectory}/${prefix}_combined.cram --fas2bam_path scripts/fas2bam.pl --samtools_path $samtools_path --bamheader windowbam.header.txt

# Prepare graph/VCF creation:
perl scripts/CRAM2VCF.pl --CRAM  ${outputDirectory}/${prefix}_combined.cram  --referenceFasta $referenceGenome --prefix ${outputDirectory}/${prefix}_finalVCF  --contigLengths  ${outputDirectory}/intermediate_files/${prefix}_postGlobalAlignment_readLengths --CRAM2VCF_executable src/CRAM2VCF --sam2alignment_executable src/sam2alignment --samtools_path $samtools_path 


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
# Launch graph/VCF generation process: (10 threads by default, executed in parallel on your local machine)
perl scripts/launch_CRAM2VCF_C++.pl --prefix ${outputDirectory}/${prefix}_finalVCF

# Create one combined graph/VCF:
perl scripts/CRAM2VCF_createFinalVCF.pl --CRAM ${outputDirectory}/${prefix}_combined.cram --referenceFasta $referenceGenome --prefix ${outputDirectory}/${prefix}_finalVCF

);	
}

								






