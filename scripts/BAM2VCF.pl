#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Getopt::Long;   
use Data::Dumper;
$| = 1;

my $referenceFasta;
my $samples;
my $BAM;
my ($regionName,$regionStart,$regionEnd);

GetOptions (
  'referenceFasta:s' => \$referenceFasta,
  'samples:s' => \$samples,
	'BAM:s' => \$BAM, 
  'name:s' => \$regionName,
  'start:i' => \$regionStart,
  'end:i' => \$regionEnd
);

die "Please specify --BAM" unless($BAM);
die "Please specify --samples" unless($samples);
die "Please specify --referenceFasta" unless($referenceFasta);
die "--BAM $BAM not existing" unless(-e $BAM);
die "--samples $samples not existing" unless(-e $samples);
die "--referenceFasta $referenceFasta not existing" unless(-e $referenceFasta);
die "please specify a region with --name, --start, --end" unless ($regionName and $regionStart and $regionEnd);

my @samples;
open(my $fh, "<", $samples);
while (<$fh>) {
  chomp;
  push @samples, $_;
}
close $fh;

printHeader();

my $sam = Bio::DB::Sam->new(-fasta => $referenceFasta, -bam => $BAM);

my $snp_caller = sub {
  my ($seqid,$pos,$p) = @_; # $pos is 1-based
  ($seqid eq $regionName and $pos >= $regionStart and $pos <= $regionEnd) or return;
  my $refbase = $sam->segment($seqid,$pos,$pos)->dna;
  my %base2name;
  for my $pileup (@$p) {
    my $b = $pileup->alignment;
    my $qname = $b->query->name;
    my $indel = $pileup->indel;
    my $qbase  = substr($b->qseq,$pileup->qpos,1);
    # If this column is an indel, return a positive integer for an insertion
    # relative to the reference, a negative integer for a deletion relative
    # to the reference, or 0 for no indel at this column.

    my $key = $qbase;
    if ($indel>0) {
      $key = substr($b->qseq,$pileup->qpos,1+$indel);
    }
    elsif ($indel<0) {
      $key = $refbase;
      $refbase = $sam->segment($seqid,$pos,$pos-$indel)->dna;
    }
    my ($sample) = $qname =~ m/(.+?)\./;
    $base2name{$key}{$sample}=$qbase;
    # print join("\t", $seqid,$pos,$refbase,$sample,$qname,$indel,$key),"\n";
  }
  outputVCF($regionName,$pos,$refbase,\%base2name) if (scalar keys %base2name > 1);
};

$sam->pileup("$regionName:$regionStart-$regionEnd", $snp_caller);

sub outputVCF {
  my ($seqid,$pos,$ref,$variation) = @_;
  my @altAlleles;
  my %sample2alt;
  for my $alt (keys %$variation) {
    if ($alt eq $ref) {
      for my $sample (keys %{$variation->{$alt}}) {
        $sample2alt{$sample} = 0;
      }
    }
    else {
      push @altAlleles, $alt;
      for my $sample (keys %{$variation->{$alt}}) {
        $sample2alt{$sample} = scalar @altAlleles;
      }
    }
  }
  my $ns = scalar keys %sample2alt;
  # my @genotypes;
  # for my $sample (@samples) {
  #   push @genotypes, exists $sample2alt{$sample} ? $sample2alt{$sample} : '.';
  # }
  print join("\t", $seqid,$pos,'.',$ref,join(',',@altAlleles),10,'PASS',"NS=$ns"),"\n"; #,"GT",join("\t",@genotypes)),"\n";
}

sub printHeader {
  my $columnsString = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
  print qq{##fileformat=VCFv4.2
##fileDate=20161026
##source=BAM2VCF.pl
##reference=file://$referenceFasta
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
$columnsString
};
}
