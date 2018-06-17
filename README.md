# Graph Genomes

Initiated at the CSHL Hackathon 24-26 October 2016

## Pipeline

* reference file GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa (GRCh38 without ALTs)
* contigs file input_sequences.fa

```
## Requires samtools version >= 1.7

bwa index GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa
bwa mem GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa input_sequences.fa  | samtools view -Sb - > input_sequences_unsorted.bam
samtools sort -o SevenGenomesPlusGRCh38Alts.bam input_sequences_unsorted.bam
samtools index SevenGenomesPlusGRCh38Alts.bam

## Check that there are no unmapped reads in the input BAM because this might lead to unknown behaviour:
samtools view -c -f 0x4 SevenGenomesPlusGRCh38Alts.bam

perl BAM2ALIGNMENT.pl --BAM SevenGenomesPlusGRCh38Alts.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta AllContigs.fa --outputFile /intermediate_files/AlignmentInput

## Expect output '../intermediate_files/AlignmentInput.sortedWithHeader'

perl checkBAM_SVs_and_INDELs.pl --BAM /data/projects/phillippy/projects/hackathon/shared/alignments/SevenGenomesPlusGRCh38Alts.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /data/projects/phillippy/projects/hackathon/shared/contigs/AllContigs.fa

perl FIND_GLOBAL_ALIGNMENTS.pl --alignmentsFile ../intermediate_files/AlignmentInput.sortedWithHeader --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --outputFile forMAFFT.bam --outputTruncatedReads ../intermediate_files/truncatedReads --outputReadLengths ../intermediate_files/postGlobalAlignment_readLengths

perl BAM2MAFFT.pl --BAM ../intermediate_files/forMAFFT.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta AllContigs.fa --outputDirectory .../intermediate_files/forMAFFT --inputTruncatedReads .../intermediate_files/truncatedReads 

perl countExpectedGlobalAlignments.pl --BAM .../intermediate_files/forMAFFT.bam

## assumes you are using SGE to submit jobs
perl CALLMAFFT.pl --action kickOff --mafftDirectory .../intermediate_files/forMAFFT --qsub 1
perl CALLMAFFT.pl --action check --mafftDirectory .../intermediate_files/forMAFFT
perl CALLMAFFT.pl --action reprocess --mafftDirectory .../intermediate_files/forMAFFT

perl globalize_windowbams.pl --fastadir /intermediate_files/forMAFFT/ --msadir /intermediate_files/forMAFFT/ --contigs /intermediate_files/postGlobalAlignment_readLengths --output /intermediate_files/combined_2.sam

perl validate_BAM_MSA.pl --BAM /intermediate_files/combined_sorted.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa

perl BAM2VCF.pl --BAM /intermediate_files/combined_sorted.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output VCF/uber_vcf.vcf

# get CRAM, via SAM sort
samtools view -h -t GRCh38.headerfile.txt /intermediate_files/combined_2.sam > /intermediate_files/combined_2_with_header.sam;\
samtools sort /intermediate_files/combined_2_with_header.sam -o /intermediate_files/combined_2_with_header_sorted.sam

cat /intermediate_files/combined_2_with_header_sorted.sam | samtools view -C -T GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa - > /intermediate_files/combined_2.cram

samtools index /intermediate_files/combined_2.cram

perl checkMAFFT_input_and_output.pl --MAFFTdir /intermediate_files/forMAFFT_2/ --contigLengths /intermediate_files/postGlobalAlignment_readLengths_2 --preMAFFTBAM /intermediate_files/forMAFFT_2.bam --finalOutputCRAM /intermediate_files/combined_2.cram &> output_checkMAFFT_v2

perl CRAM2VCF.pl --CRAM /intermediate_files/combined_2.cram --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output VCF/graph_v2.vcf --contigLengths /intermediate_files/postGlobalAlignment_readLengths_2

perl CRAM2VCF_checkVariantDistribution.pl --output VCF/graph_v2.vcf

perl launch_CRAM2VCF_C++.pl --output VCF/graph_v2.vcf

perl CRAM2VCF_createFinalVCF.pl --CRAM /intermediate_files/combined_2.cram --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output VCF/graph_v2.vcf

```


## Workflow
### Step 0: Align all contigs to GRCh38
_This step is outside of our formalized pipeline. Users will be responsible for
running BWA (or equivalent tool) to produce BAM alignment from a reference
genome and additional assemblies_
+ Input:
  + Reference genome (GRCh38)
    + Format: _FASTA_ (text)
  + Individual contigs
    + Format: _FASTA_ (text)
+ Output:
  + Single file of local alignments
    + Format: _BAM_ (binary)
+ Tools:
  + BWA

### Step 1: Globally align all contigs
+ Input:
  + Reference genome (GRCh38)
    + Format: _FASTA_ (text)
  + Individual contigs
    + Format: _FASTA_ (text)
  + Local alignment (output of step 0)
    + Format: _BAM_ (binary)
+ Output:
  + Global alignment where each contig has a single alignment
    + Format: _BAM_ (binary)
+ Tools:
  + Local To Global Alignment
    + Assignees:
      + Evan
      + Aarti
    + Status:
      + In progress
      + scripts/NeedlemanWunsch.ph & scripts/NW_alignment_scratch.ph on GitHub
  + Wrappers
    + Assignees:
      + Nathan?
      + Jeff
    + Status: in progress
      + scripts/gg-01.1-local.sh on GitHub

### Step 2: Divide and conquer multiple sequence alignment
+ Input:
  + Global alignment of all contigs (output of step 1)
    + Format: _BAM_ (binary)
+ Output:
  + Single global multiple sequence alignment
    + Format: _BAM_ (binary)
+ Tools:
  + MAFFT
    + Status: _complete_ (installed)
  + AMC
    1. Identify window coordinates (start with windows of size ~10kb) in
    alignment file
      + Assignees:
        + Alex
      + Status: in progress
        + BAM2MAFFT.pl in top-level of GitHub repo
    2. Extract sequence data for each window
      + Assignees:
        + Alex
      + Status: in progress
        + BAM2MAFFT.pl in top-level of GitHub repo
    3. Convert window to FASTA format
      + Assignees:
        + Alex
      + Status: in progress
        + BAM2MAFFT.pl in top-level of GitHub repo
    4. Run MAFFT for each window [could be parallelized]
      + Assignees:
        + Nathan
      + Status: in progress
        + scripts/align_mafft.sh on GitHub
    5. Convert MAFFT-produced FASTA alignments to BAM files (one for each
    window)
      + Assignees:
        + Nancy
      + Status: _complete_
        + scripts/fas2bam.pl
    6. Reassemble the individual BAM MSAs into a single BAM alignment file
      + Assignees:
        + Nancy
      + Status: in progress
        + _Location for any code?_
  + Wrappers
    + Window size specification?
    + Call AMC
    + Assignees
      + Jeff?
    + Status: **not started**

### Step 3: Create graph genome!
+ Input:
  + Single global multiple sequence alignment (output of step 2)
    + Format: _BAM_ (binary)
+ Output:
  + Graph representation of variation in reference genome + input contigs
    + Format: _VG_, _GFA_?
+ Tools:
  + BAM to VCF
    + Assignees:
      + Andrew
    + Status: in progress
      + _Location for any code?_
  + Wrapper
    + Convert BAM alignment to VCF
    + Call vg to convert VCF to vg or gfa format
    + Assignees:
      + Jeff
    + Status: **not started**


## Contributors

Evan Biederstedt (NYGC, WCM)

Alexander Dilthey (NIH)

Nathan Dunn (LBNL)

Aarti Jajoo (Baylor)

Nancy Hansen (NIH)

Jeff Oliver (Arizona)

Andrew Olsen (CSHL)
