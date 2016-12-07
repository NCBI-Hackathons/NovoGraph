#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
$| = 1;

my $referenceFasta;
my $BAM;
my $window;
my $samples;
my $script;
my $outputDirectory = 'forVG';
my $threads = 8;
GetOptions (
  'outputDirectory:s' => \$outputDirectory,	
  'referenceFasta:s' => \$referenceFasta,
  'BAM:s' => \$BAM,
  'samples:s' => \$samples,
  'window:i' => \$window,
  'script:s' => \$script
);

die "Please specify --script <path to BAM2VCF.pl>" unless($script);
die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "Please specify --window 10000" unless($window);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);
die "--samples $samples not existing" unless(-e $samples);

unless($outputDirectory =~ /forVG/)
{
	die "Safety check failure - I want to rm -rf $outputDirectory, but I can't because it doesn't contain the string forVG";
} 
my $cmd_delete = qq(rm -rf $outputDirectory);
if(system($cmd_delete))
{
	die "Command $cmd_delete failed.";
}
unless((-e $outputDirectory) and (-d $outputDirectory))
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory directory";
}

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my @sequence_ids = $sam->seq_ids();
my @sequence_lengths = @{$sam->header->target_len};

my @commands;
my @output_files;
for(my $i=0;$i<@sequence_ids;$i++) {
  my $name = $sequence_ids[$i];
  my $length = $sequence_lengths[$i];
  
  # check if there are alignments
  my @alignments = $sam->features(-seq_id => $name);
  next unless (@alignments);

  my $pieceDirectory = "$outputDirectory/$name";
  unless((-e $pieceDirectory) and (-d $pieceDirectory)) {
    mkdir($pieceDirectory) or die "cannot mkdir $pieceDirectory";
  }
  for(my $j=1;$j<$length;$j+=$window) {
    my $end = $j+$window-1;
    if ($end > $length) {
      $end = $length;
    }
	
    @alignments = $sam->features(-seq_id => $name, -start => $j, -end => $end);
    next unless (@alignments);

    my $padded_j = sprintf("%09d", $j);
    my $outFile = "$pieceDirectory/$padded_j.vcf";
    push(@commands, "$script --BAM $BAM --referenceFasta $referenceFasta --samples $samples --name $name --start $j --end $end > $outFile\n");
	
	push(@output_files, $outFile);
  }
}

my $fn_commands = '_cmd_for_BAM2VCF';
open(CMD, '>', $fn_commands) or die "Can't open file $fn_commands";
print CMD join("\n", @commands), "\n";
close(CMD);

print "\nNow process ", scalar(@commands), " commands.\n\n";
my $parallel_cmd = qq(parallel -a $fn_commands -j $threads);
if(system($parallel_cmd))
{
	die "Command $parallel_cmd failed";
}

print "\nparallel exited with code 0\n";

my $output_fn = $outputDirectory . '/combined.vcf';
open(OUT, '>', $output_fn) or die "Cannot open $output_fn";
my $columnsString = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
print OUT qq{##fileformat=VCFv4.2
##fileDate=20161026
##source=BAM2VCF.pl
##reference=file://$referenceFasta
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n};
print OUT $columnsString, "\n";
close(OUT);

foreach my $expected_output_file (@output_files)
{
	die "Expected file $expected_output_file not present" unless(-e $expected_output_file);
	my $cmd_append = qq(cat $expected_output_file >> $output_fn);
	if(system($cmd_append))
	{
		die "Command $cmd_append failed.";
	}
}	

print "\nDone. Produced output file $output_fn\n\n";