# Arguments for each step of pipeline

1. Local to Global
  1. Pre-processing
     + Coder: Dilthey
     + scripts/BAM2MAFFT.pl
       + Example call: ./BAM2MAFFT.pl --BAM /home/data/alignments/statistics/SevenGenomesGlobalAligns.bam --referenceFasta /home/data/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /home/data/contigs/AllContigs.fa
       + Arguments:
         + --BAM
            + Description: input BAM file, produced by BWA
         + --referenceFasta
           + Description: Reference file used to produce the BAM. If GRCh37/38, no ALTs!
         + --paranoid 0/1
           + Description: Check that the emitted sequence coordinates and alignments are concordant with the input sequences that went into the BAM. Memory-intensive!
         + --readsFasta
            + Description: If --paranoid 1, path to FASTA with sequences that went into the BAM
  + 2. Local to global	   
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
        + --BAM
          + Description: BAM alignment from BWA
        + --outputDirectory
          + Description: Directory for output files (optional)
	+ --paranoid 0/1
	  + Description: Check that original input sequences (the sequences that went into the BAM) can be reconstructed from window sequences. Memory-intensive! Default: off.
	+ --readsFasta
	  + Description: If --paranoid 1, FASTA file with all sequences that went into the BAM.
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
    		  + Description: Name of list of samples [assemblies] (one per line) an
          example of this is in examples/list_of_assemblies.txt
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

# QC scripts
1. Comparing two FASTA files
  + Coder: Dilthey
  + scripts/compareTwoFASTAs.pl
    + Example call: `./compareTwoFASTAs.pl --f1 /home/devsci7/globalize_windowbams/global_multiple_alignments.frombam.fasta --f2 /home/data/contigs/AllContigs.fa`
    + Arguments:
      + --f1 File 1 for comparison (kept in memory)    
      + --f2 File 2 for comparison (read iteratively)

2. Making sure that per-window FASTAs (or MFAs) contain the original input sequences
  + Coder: Dilthey
  + scripts/checkWindowCorrectness.pl   
    + Example call: ./checkWindowCorrectness.pl --contigsFasta /home/data/contigs/AllContigs.fa --fileWithAllGapWindows ../forMafft2/_alignments_inWindow_onlyGaps --windowFilesDirectory ../forMafft2 --windowsInfo ../forMafft2/_windowsInfo
    + Arguments:
      + --contigsFasta File with original input sequences (that went into the original BAM)
      + --fileWithAllGapWindows File with windows that are all-gap (generated by BAM2MAFFT.pl)
      + --windowFilesDirectory Directory with per-window files (generated by BAM2MAFFT.pl)
      + --windowsInfo File with window info (generated by BAM2MAFFT.pl)
