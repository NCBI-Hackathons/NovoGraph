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

globalize_windowbams.pl - Given a set of BAM files that are windowed to a particular genomic region, combine them into a global BAM file by changing positions to reflect offsets in a provided offset file.  In addition, combine entries for single contigs present in multiple windows. 

=head1 SYNOPSIS

Create a global BAM file from numerous local (windowed) BAM files, combining entry alignments as necessary.

=head1 USAGE

globalize_windowbams.pl --fastadir <path to directory with inputs for MAFFT> --msadir <path to directory with outputs from MAFFT> --contigs <path to file with contig names/lengths>

=head1 DESCRIPTION

This script reads in a ".windowInfo" file containing tab-delimited lines (one for each windowed BAM file) in which the first field is a reference entry, the second field is the a window id, and the third field is the offset along the entry for the window (i.e., the position in the output file = position in window BAM file + offset), as well as a ".contigs" file containing the names of all contigs in all the windowed BAM files, along with their original lengths (i.e., the length of the entire sequence that has been split to different windows).  It outputs a BAM-formatted file (with name specified with the --output option) in which alignments for single contigs that span multiple windows are combined into a single contiguous alignment.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $fastadir = $Opt{'fastadir'};
my $msadir = $Opt{'msadir'};
my $contigs_file = $Opt{'contigs'};
my $output_file = $Opt{'output'};

my $rh_windowinfo = read_windowbams_info($fastadir, $msadir); # hashed by entry, then sorted lists of hashes to bam paths ("bam") and offset positions ("offset")
my $rh_contiglength = read_contigs_file($contigs_file); # contig lengths hashed by contig name

my $sambam_file = ($Opt{bamheader}) ? " | samtools view -uS -t $Opt{bamheader} > $output_file"
                                    : "> $output_file";

my $sam_fh = FileHandle->new(">$sambam_file");
                   
foreach my $entry (sort keys %{$rh_windowinfo}) {
    my %current_contigs = (); # if contig lengths are right, this will never get too big

    foreach my $rh_baminfo (@{$rh_windowinfo->{$entry}}) {
        my $bam = $rh_baminfo->{bam};
        my $offset = $rh_baminfo->{offset};

        my $rh_window_contig_alignments = offset_windowbam($bam, $entry, $offset);

        foreach my $contig (keys %{$rh_window_contig_alignments}) {
            if ($current_contigs{$contig}) {
                $current_contigs{$contig} = combine_contig_alignments($current_contigs{$contig},
                    $rh_window_contig_alignments->{$contig});
            }
            else {
                $current_contigs{$contig} = $rh_window_contig_alignments->{$contig};
            }

            if (length($current_contigs{$contig}->{seq})==$rh_contiglength->{$contig}) {
                print_sam_entry($sam_fh, $current_contigs{$contig});
                delete $current_contigs{$contig};
            }
        }
    }
}

close $sam_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    
    # Set defaults here
    %Opt = ( output => 'uber.bam' );
    GetOptions(\%Opt, qw( fastadir=s msadir=s contigs=s bamheader=s output=s
                manual help+ version
                verbose 
                )) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "globalize_windowbams.pl, ", q$Revision: $, "\n"; }

    if (!$Opt{fastadir}) { die "Must specify input fasta directory with option --fastadir.\n"; }
    if (!$Opt{msadir}) { die "Must specify output msa directory with option --msadir.\n"; }
    if (!$Opt{contigs}) { die "Must specify input .contigs filename with option --contigs.\n"; }

}

# Note: This routine assumes the format of the file created by Alex D. at the CSHL Hackathon
sub read_windowbams_info {
    my $fastadir = shift;
    my $msadir = shift;

    my %entry_hash = ();
    open WINDOWS, "$fastadir/_windowsInfo"
        or die "Couldn\'t open $fastadir/_windowsInfo: $!\n";
    while (<WINDOWS>) {
        next if (/^referenceContigID/);
        if (/^(\S+)\s(\d+)\s(\d+)\s(\d+)/) {
            my ($entry, $windowid, $start, $end) = ($1, $2, $3, $4);
            push @{$entry_hash{$entry}}, {bam => "$msadir/$entry\_$windowid\_aligned.fa.bam", 
                                          offset => $start};
        }
        else {
            die "Illegal format in _windowsInfo file:\n$_";
        }
    }
    close WINDOWS;

    return {%entry_hash};
}

sub read_contigs_file {
    my $file = shift;

    my %contig_length = ();
    open CONTIGS, $file
        or die "Couldn\'t open .contigs file $file: $!\n";
    while (<CONTIGS>) {
        if (/^(\S+)\s(\d+)$/) {
            my ($contig, $contiglength) = ($1, $2);
            $contig_length{$contig} = $contiglength;
        }
        else {
            die "Illegal format in contigs file:\n$_";
        }
    }
    close CONTIGS;

    return {%contig_length};
}

sub offset_windowbam {
    my $bamfile = shift;
    my $entry = shift;
    my $offset = shift;

    my $samtools_view_cmd = "samtools view $bamfile | "; 
    my $fh = FileHandle->new($samtools_view_cmd); # assume SAM format for now

    my %aligns = ();
    while (<$fh>) {
        if (/^(\S+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)/) {
            my ($contig, $flag, $refentry, $position, $mapqual, $cigar, $seq, $qual) = 
                ($1, $2, $3, $4, $5, $6, $7, $8);
            my $globalpos = $position + $offset;
            my $final_covered_pos = calc_final_refpos($globalpos, $cigar);
            $aligns{$contig} = {contig => $contig,
                                flag => $flag,
                                refentry => $entry, # use chromosome, not REF in BAM file
                                pos => $globalpos,
                                mapqual => $mapqual,
                                cigar => $cigar,
                                seq => $seq,
                                qual => $qual,
                                final_refpos => $final_covered_pos };
        }
    }
    close $fh;

    return {%aligns};
}

sub combine_contig_alignments {
    my $rh_first_align = shift;
    my $rh_second_align = shift;

    my $contig = $rh_first_align->{contig};
    my $flag = $rh_first_align->{flag};
    my $refentry = $rh_first_align->{refentry};
    my $globalpos = $rh_first_align->{pos};
    my $finalpos1 = $rh_first_align->{final_refpos};
    my $mapqual = $rh_first_align->{mapqual};
    my $qual = $rh_first_align->{qual};
    my $cigar1 = $rh_first_align->{cigar};
    my $seq1 = $rh_first_align->{seq};
    my $cigar2 = $rh_second_align->{cigar};
    my $seq2 = $rh_second_align->{seq};
    my $pos2 = $rh_second_align->{pos};
    my $finalpos2 = $rh_second_align->{final_refpos};

    my $combined_cigar = $cigar1;
    my $deleted_length = $pos2 - $finalpos1 - 1;
    if ($deleted_length) {
        $combined_cigar .= $deleted_length.'D';
    }
    $combined_cigar .= $cigar2;
    my $combined_seq = $seq1.$seq2;

    return {contig => $contig,
            flag => $flag,
            refentry => $refentry,
            pos => $globalpos,
            mapqual => $mapqual,
            cigar => $combined_cigar,
            seq => $combined_seq,
            qual => $qual,
            final_refpos => $finalpos2 };
}

sub print_sam_entry {
    my $fh = shift;
    my $rh_align = shift;

    my $contig = $rh_align->{contig};
    my $flag = $rh_align->{flag};
    my $entry = $rh_align->{entry};
    my $pos = $rh_align->{pos};
    my $mapqual = $rh_align->{mapqual};
    my $cigar = $rh_align->{cigar};
    my $seq = $rh_align->{seq};
    my $qual = $rh_align->{qual};

    print $fh "$contig\t$flag\t$entry\t$pos\t$mapqual\t$cigar\t$seq\t$qual\n"; 
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=over 4

=item B<--fastadir>

The path to a directory containing all the input to MAFFT.  This should consist of fasta
files created for MAFFT with names like "A.XXXX" where "A" is the chromosome and XXXX is 
the starting position of the window.  This directory must also contain a "_windowsInfo"
file containing the offsets of each window.

=item B<--msadir>

The path to a directory containing all the output of MAFFT.  This should consist of fasta
files created by MAFFT with names like "A_XXXX_aligned.fa" where "A" is the chromosome 
and XXXX is the starting position of the window.

=item B<--output>

The name of the SAM or BAM file to be created (default = uber.bam).

=item B<--contigs>

The path to a file containing the length of each contig in the MAFFT alignments.

=item B<--bamheader>

Optional path to a header file to be used to create BAM format.  If not provided,
alignments (without a header) will be output in SAM format.

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
