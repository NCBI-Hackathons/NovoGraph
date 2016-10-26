#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
$| = 1;

my $referenceFasta;
my $BAM;
my $outputDirectory = 'forVG';
my $output;

GetOptions (
  'outputDirectory:s' => \$outputDirectory,	
  'output:s' => \$output,
  'referenceFasta:s' => \$referenceFasta,
	'BAM:s' => \$BAM
);

die "Please specify --BAM" unless($BAM);
die "Please specify --referenceFasta" unless($referenceFasta);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);

unless((-e $outputDirectory) and (-d $outputDirectory))
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory directory";
}

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my @pieces;
for my $seqName ($sam->seq_ids()) {

  # check if there are alignments
  my @alignments = $sam->features(-seq_id => $seqName);
  next unless (@alignments);

  my $pieceDirectory = "$outputDirectory/$seqName";
  push @pieces, "$pieceDirectory/*.vcf.gz";
}
#my $cmd = "vcf-concat @pieces | bgzip -c > $output; tabix -p $output\n";
my @outformat = 'vcf';
my $cmd = "vcf-concat @pieces | bgzip -c > $output; tabix -p @outformat\n";
system($cmd);
