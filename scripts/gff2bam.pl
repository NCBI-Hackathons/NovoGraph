#!/usr/bin/perl -w
# $Id: $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FileHandle;

our %Opt;

my $Id = q$Id:$;

=head1 NAME

gff2bam.pl - script to convert an NCBI gff file to SAM/BAM formatted files.

=head1 SYNOPSIS

This script reads in NCBI gff files containing alignments for the GRCh38 "alt" contigs to the primary reference and converts them to BAM format.

=head1 USAGE

gff2bam.pl --input <file of gff file locations> --altfasta <fasta file of alt sequences> --bamheader <path to file containing header for BAM file> --output <path to write output SAM or BAM file to>

=head1 DESCRIPTION

Each gff file listed in the input file specified with --input is parsed to create a single SAM-formatted alignment to the primary reference.  If given the path of a BAM header file with the --bamheader option, the script writes the alignments in BAM format.  Otherwise, it will write the alignments in SAM format.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $input_file = $Opt{'input'};
my $chrom_file = $Opt{'chrnames'};
my $altfasta_file = $Opt{'altfasta'};
my $output_file = $Opt{'output'};

my $ra_gff_files = read_input_file($input_file);
my $rh_chrom_names = read_chrom_names($chrom_file);


if (!$output_file) { # create a name from the input file
    my $filebase = $input_file;
    $filebase =~ s:.*/::; # put it in the current directory
    $output_file = ($Opt{bamheader}) ? "$filebase.bam" : "$filebase.sam";
}

my $sambam_file = ($Opt{bamheader}) ? " | samtools view -uS -t $Opt{bamheader} > $output_file"
                                    : "> $output_file";

my $sam_fh = FileHandle->new("$sambam_file");

my $no_files = @{$ra_gff_files};
print "$no_files gff files!\n";
foreach my $gff_file (@{$ra_gff_files}) {
    my $samstring = retrieve_alignment($gff_file, $rh_chrom_names);
    print $sam_fh "$samstring";
}
close $sam_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( input=s bamheader=s output=s chrnames=s altfasta=s
                manual help+ version
                verbose 
                )) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "gff2bam.pl, ", q$Revision: $, "\n"; }

    if (!$Opt{input}) { die "Must specify input filename with option --input.\n"; }
    if (!$Opt{chrnames}) { die "Must specify filename of file of chromosome names with option --chrnames.\nTry ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.28.assembly.txt"; }
    if (!$Opt{altfasta}) { die "Must specify location of indexed fasta file containing alt sequences with option --altfasta.\nTry ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.24_GRCh38.p9/GCA_000001405.24_GRCh38.p9_genomic.fna.gz\n"; }

}

sub read_input_file {
    my $file = shift;

    my @gff_list = ();
    open GFFS, $file
        or die "Couldn\'t open file $file: $!\n";
    while (<GFFS>) {
        if (/^(\S+)/) {
            push @gff_list, $1;
        }
        else { # assume sequence line
            die "Illegal format!\n";
        }
    }
    close GFFS;

    return [@gff_list];
}

sub read_chrom_names {
    my $file = shift;

    my $fh = FileHandle->new("$file");

    my %chrom_names = ();
    while (<$fh>) {
        next if (/^#/);
        chomp;
        my @fields = split /\t/, $_;
        $fields[$#fields] =~ s/\s+$//;
        $chrom_names{$fields[4]} = $fields[$#fields];
        $chrom_names{$fields[6]} = $fields[$#fields];
    }

    return {%chrom_names};
}

sub retrieve_alignment {
    my $gff_file = shift;
    my $rh_chrname = shift;

    my $gff_fh = FileHandle->new("$gff_file");

    while (<$gff_fh>) {
        next if (/^#/);
        chomp;
        my @fields = split /\t/, $_;
        my $contig = $fields[0];
        my $chrom = $rh_chrname->{$contig};
        if (!$chrom) {
            die "$contig does not have a chrom name!\n";
        }

        my $chromstart = $fields[3];
        my $chromend = $fields[4];
        my $field6 = $fields[5];
        my $info_field = $fields[8];
        #print "$info_field!\n";
        my @infotags = split /;/, $info_field;
        foreach my $infotag (@infotags) {
            print "INFO: $infotag\n";
        }
        my ($target_tag) = grep /Target=/, @infotags;
        my ($target, $targetstart, $targetend);
        if ($target_tag =~ /Target=(\S+)\s(\d+)\s(\d+)/) {
            ($target, $targetstart, $targetend) = ($1, $2, $3);
        }
        else {
            die "Unable to parse target info from tag $target_tag!\n";
        }

        #print STDERR "$gff_file\n$infotags[$#infotags]\n";
        my ($cigar_tag) = grep /Gap=/, @infotags;
        $cigar_tag = '*' if (!$cigar_tag);
        print "CIGAR: $cigar_tag\n";
        #if ($cigar_tag =~ /Gap=(.*)$/) {
        if ($info_field =~ /Gap=(.*)$/) {
            my $cigar = $1;
            $cigar =~ s/^\s+//;
            $cigar =~ s/\s+$//;
            my @cigar_ops = split /\s+/, $cigar;
            my $newcigar = '';
            my ($hapbases, $primarybases) = (0,0); # count bases in haplotype and in primary chromosome
            foreach my $op (@cigar_ops) {
                if ($op =~ /^([IMD])(\d+)$/) {
                    my ($type, $nbases) = ($1, $2);
                    if ($type =~ /^[IM]$/) {
                        $hapbases += $nbases;
                    }
                    if ($type =~ /^[DM]$/) {
                        $primarybases += $nbases;
                    }
                    $newcigar .= $nbases.$type; 
                }
            }
            my $targetseq = `samtools faidx $altfasta_file $target:$targetstart-$targetend`;;
            $targetseq =~ s/\>[^\n]*\n//;
            $targetseq =~ s/\n//g;
            $targetseq = uc($targetseq);
            #print "Retrieved seq for $target:\n$targetseq\n";
            print "$target $targetstart-$targetend matches $chrom $chromstart-$chromend\n";
            return "$target\t0\t$chrom\t$chromstart\t0\t$newcigar\t*\t0\t0\t$targetseq\t*\n";
        }
    }
    close $gff_fh;
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=over 4

=item B<--input>

The path to a file of gff file paths.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
