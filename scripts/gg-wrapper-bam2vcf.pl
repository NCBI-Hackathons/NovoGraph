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

unless((-e $outputDirectory) and (-d $outputDirectory))
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory directory";
}

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my @sequence_ids = $sam->seq_ids();
my @sequence_lengths = @{$sam->header->target_len};


my %pieces;
for(my $i=0;$i<@sequence_ids;$i++) {
  my $pieceDirectory = "$outputDirectory/$sequence_ids[$i]";
  unless((-e $pieceDirectory) and (-d $pieceDirectory)) {
    mkdir($pieceDirectory) or die "cannot mkdir $pieceDirectory";
  }
  for(my $j=1;$j<$sequence_lengths[$i];$j+=$window) {
    my $end = $j+$window;
    if ($end > $sequence_lengths[$i]) {
      $end = $sequence_lengths[$i];
    }
    my $padded_j = sprintf("%09d", $j);
    my $outFile = "$pieceDirectory/$padded_j.vcf.gz";
    print "$script --BAM $BAM --referenceFasta $referenceFasta --samples $samples --name $sequence_ids[$i] --start $j --end $end | bgzip -c > $outFile\n";
  }
}

