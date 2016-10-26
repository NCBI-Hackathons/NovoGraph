# Arguments for each step of pipeline

1. Local to Global
  + Coders: Biederstedt & Jajoo
  + scripts/NWAM_01.py
    + Arguments:
      + [parameter name/flag]
        + Description: FASTA file of reference genome
      + [parameter name/flag]
        + Description: Directory of directories with FASTA files of contigs;
        one directory corresponds to one assembly **TODO**: Check
      + [parameter name/flag]
        + Description: BAM alignment from BWA? **TODO**: Check
2. Multiple Sequence alignment
  1. Windowing
    + Coder: Dilthey
    + scripts/BAM2MAFFT.pl
      + Example call: `BAM2MAFFT.pl --referenceFasta data/fastas/GRCh38.fas
      --BAM data/bams/bam-from-bwa.bam`
      + Arguments:
        + --referenceFasta
          + Description: FASTA file of reference genome
        + [parameter name/flag] **TODO**: Check this - is it necessary?
          + Description: Directory of directories with FASTA files of contigs;
          one directory corresponds to one assembly **TODO**: Check
        + --BAM
          + Description: BAM alignment from BWA? **TODO**: Check
  2. Multiple Sequence Alignment via MAFFT
    1. MAFFT
      + Coder: Dunn
      + scripts/gg-align_mafft.sh
        + Example call: `./gg-align_mafft.sh -i data/fastas -o data/bams`
        + Arguments:
          + -i
            + Description: path to FASTA input files.
          + -o
            + Description: path to output aligned FASTA and BAMS.
    2. FASTA to BAM
      + Coder: Hansen
      + scripts/fas2bam.pl (called by scripts/gg-align_mafft.sh)
        + Example call: `fas2bam.pl --input data/fastas/fasta01.fas --output
        data/bams --ref "ref" --bamheader "./config/windowbam.header.txt"`
        + Arguments:
          + --input
            + Description: Name of input FASTA file produced by MAFFT
          + --output
            + Description: Name of output BAM file
          + --ref
            + Description: String identifier for reference sequence, defaults to
            "ref"
          + --bamheader
            + Description: BAM header file; currently one in
            config/windowbam.header.txt
  3. Reassemble aligned Genomes
    + Coder: Hansen
    + scripts/globalize_windowbams.pl
      + Example call: `globalize_windowbams.pl --fastadir data/fastas --msadir
      data/bams --contigs data/contiginfo.txt` --bamheader config/GRCh38_header.fai --output global_msa.bam
      + Arguments:
        + --fastadir
          + Description: path to directory with inputs for MAFFT (this is the
            same as argument -i fed to gg-align_mafft.sh - see 2.2.1)
        + --msadir
          + Description: path to directory with outputs from fas2bam.pl (files
            are assumed to be named "chr#_#.bam".  This directory is the same
            as the output directory fed to gg-align_mafft.sh as -o - see 2.2.1)
        + --contigs
          + Description: path to file with contig names/lengths (This file
            can be obtained from Alex's contig info file which is used by
            Evan and Aarti's software.)
        + --output
          + Description: path to where a bamfile should be written
        + --bamheader
          + Description: path to a header file for the entire genome
3. Graph genome
  1. BAM to VCF:
    + Coder: Olson
    + scripts/gg-wrapper-bam2vcf.pl
      + Example call: `gg-wrapper-bam2vcf.pl --script BAM2VCF.pl --BAM data/bams/uberbam.bam --referenceFasta
      data/reference/GRCh38.fa --window 100000 --samples samples.txt | parallel -j 8`
      + Arguments:
	    + --script
		  + The perl script that does the work
        + --BAM
          + Description: BAM file (the file produced by 2.3)
        + --referenceFasta
          + Description: Name of indexed reference genome fasta file
        + --window
          + Description: Size of window for parallelization (try 100000)
		+ --samples
		  + Description: Name of list of samples (one per line)
  2. VCF to VG
    + vg (no custom scripts)
      + Example call: `vg construct -r data/reference/GRCh38.fa -v data/vcf/final.vcf.gz -p >
      data/graph/graph.vg`
      + Arguments:
        + -r
          + Description: Reference genome fasta file
        + -v
          + Description: tabix indexed VCF file of aligned assemblies (product of 3.1)
    		+ -p
    		  + Description: Show progress on STDERR
        + [no flag]
          + Description: Name of output graph file
