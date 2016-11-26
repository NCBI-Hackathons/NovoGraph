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

fas2bam.pl - script to convert a multiple alignment (including a reference entry) to discover variants and write both SAM/BAM (and maybe VCF eventually) formatted files.

=head1 SYNOPSIS

This script reads in an MSA file with padded entries and a reference entry representing a multiple alignment (as produced by MAFFT, for example).  It then writes a BAM (by default) or SAM file with the alignments of each non-reference entry to the reference.

=head1 USAGE

fas2bam.pl --input <path to fasta file with multiple alignment> --ref <reference entry id> --bamheader <path to file containing header for BAM file> --output <path to write output SAM or BAM file to>

=head1 DESCRIPTION

When fed a multiple sequence alignment in FASTA format and a specification of which entry in this file is the "reference", this script will compare each non-reference entry to the reference entry to produce a SAM or BAM-formatted alignments for them.  If given the path of a BAM header file with the --bamheader option, the script writes the alignments in BAM format.  Otherwise, it will write the alignments in SAM format.

=cut

#------------
# Debug
#------------

my $last_call_read_fas_file_read_nucleotides = 0;

#------------
# Begin MAIN 
#------------

process_commandline();

my $input_file = $Opt{'input'};
my $ref_entry = $Opt{'ref'};
my $output_file = $Opt{'output'};

my $rh_entry_seqs = read_fas_file($input_file);

#check for ref entry
if (!$rh_entry_seqs->{$ref_entry}) {
	my $l_references = join("\n", map {"'" . $_ . "'"} keys %$rh_entry_seqs);
    die "No entry for reference $ref_entry in $input_file - read $last_call_read_fas_file_read_nucleotides sequence characters - have the following references:\n".$l_references;
}

my $ref_length = length($rh_entry_seqs->{$ref_entry});

# check that all entries are same length:
foreach my $entry (sort keys %{$rh_entry_seqs}) {
    next if ($entry eq $ref_entry);
    my $seq = $rh_entry_seqs->{$entry};
    my $seqlength = length($seq);
    if ($seqlength != $ref_length) {
        die "$entry and $ref_entry have unequal lengths ($seqlength, $ref_length)!\n";
    }
}

if (!$output_file) { # create a name from the input file
    my $filebase = $input_file;
    $filebase =~ s:.*/::; # put it in the current directory
    $output_file = ($Opt{bamheader}) ? "$filebase.bam" : "$filebase.sam";
}

my $sambam_file = ($Opt{bamheader}) ? " | samtools view -uS -t $Opt{bamheader} > $output_file"
                                    : "> $output_file";

my $sam_fh = FileHandle->new("$sambam_file");

                   
my ($flag, $score) = (0, 0);
foreach my $entry (sort keys %{$rh_entry_seqs}) {
    next if ($entry eq $ref_entry);

    my ($sam_start, $cigar_string, $rs_entryseq, $ra_vars) = parse_alignment($rh_entry_seqs, $entry, $ref_entry);
    print $sam_fh "$entry\t$flag\t$ref_entry\t$sam_start\t$score\t$cigar_string\t*\t0\t0\t$$rs_entryseq\t*\n";
}

#------------
# End MAIN
#------------

sub process_commandline {
    
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( ref=s input=s bamheader=s output=s
                manual help+ version
                verbose 
                )) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "fas2bam.pl, ", q$Revision: $, "\n"; }

    if (!$Opt{input}) { die "Must specify input filename with option --input.\n"; }
    if (!$Opt{ref}) { die "Must specify the name of the reference entry in the input file with option --ref.\n"; }

}

sub read_fas_file {
    my $file = shift;

	$last_call_read_fas_file_read_nucleotides = 0;
    my %seq_entry_hash = ();
    my ($current_entry, $current_seq);
    open FAS, $file
        or die "Couldn\'t open file $file: $!\n";
    while (<FAS>) {
        if (/^\>\s*(\S+)(\s+.*){0,1}$/) {
            my ($id, $desc) = ($1, $2);
            if ($current_entry && $current_seq) {
                $seq_entry_hash{$current_entry} = $current_seq;
                #print "$current_seq\n";
            }
            $current_entry = $id;
            $current_seq = '';
        }
        else { # assume sequence line
            chomp;
            $current_seq .= $_;
            if (!$current_entry || $_ =~ /([^ATGCatgcNn\-RrYySsWwKkMmBb])/) {
                die "Illegal format (no entry definition line - $current_entry) or illegal characters (not ATGCNatgcn or - -- $1)!  Aborting.\n";
            }
			
			$last_call_read_fas_file_read_nucleotides += length($_);
        }
    }
    close FAS;

    if ($current_entry && $current_seq) {
        $seq_entry_hash{$current_entry} = $current_seq;
        #print "$current_seq\n";
    }

    return {%seq_entry_hash};
}

sub parse_alignment {
    my $rh_entries = shift;
    my $this_entry = shift;
    my $ref_entry = shift;

    my $reverseentryseq = reverse($rh_entries->{$this_entry}); # reversed so we can use chop to step through
    my $reverserefseq = reverse($rh_entries->{$ref_entry});

    my $flat_cigar = '';
    my $start_pos = 1; # this will shift up or down if there are pads at the beginning of the entry or reference, resp.
    my $entry_seq = '';
    my @variants = (); # not dealing with VCF yet

    my $match_seen_yet = 0; # entry gaps ahead of aligned bases will not be in the cigar bases, but change start position
    my ($next_entrybase, $next_refbase) = (chop $reverseentryseq, chop $reverserefseq);
    while ($next_entrybase && $next_refbase) {
        if (($next_entrybase ne '-') && ($next_refbase ne '-')) { # aligned bases
            $match_seen_yet = 1;
            $flat_cigar .= 'M';
            $entry_seq .= uc($next_entrybase);
        }
        elsif (($next_entrybase eq '-')  && ($next_refbase ne '-')){
            if (!$match_seen_yet) {
                $start_pos++;
            }
            else {
                $flat_cigar .= 'D';
            }
        }
        elsif (($next_refbase eq '-') && ($next_entrybase ne '-')) {
            $flat_cigar .= 'I';
            $entry_seq .= uc($next_entrybase);
        }
        else {
            $flat_cigar .= 'P';
        }
        ($next_entrybase, $next_refbase) = (chop $reverseentryseq, chop $reverserefseq);
    }

    #print "$flat_cigar\n";
    my $cigar_string = condense_flat_cigar($flat_cigar);

    return ($start_pos, $cigar_string, \$entry_seq, [@variants]);
}

sub condense_flat_cigar {
    my $flat_string = shift;

    my $reverse_flatcigar = reverse($flat_string);

    my $cigar_string = '';
    my $current_op = '';
    my $no_currentops = 0;
    my $nextop;
    while ($nextop = chop $reverse_flatcigar) {
        if ($current_op && $nextop ne $current_op) { # write out count of old ops
            $cigar_string .= $no_currentops.$current_op;
            $current_op = $nextop;
            $no_currentops = 1;
        }
        elsif (!$current_op) { # must be first!
            $current_op = $nextop;
            $no_currentops = 1;
        }
        else { # continuation of same operator
            $no_currentops++;
        }
    }
    if ($current_op) {
        $cigar_string .= $no_currentops.$current_op;
    } 

    return $cigar_string;
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

The path to an MSA file with FASTA-formatted entries that are padded with "-" 
characters to form a multiple alignment.

=item B<--ref>

The id of one of the entries of the multiple alignment which serves as the 
reference entry (for purposes of describing variants and alignments).

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
