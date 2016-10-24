# Graph Genomes
## CSHL Hackathon 24-26 October 2016

## Preliminary description of tools, inputs, and outputs
### Problem: Whole genome graph creation from multiple assemblies

### Step 1: Align all contigs to GRCh38
+ Input:
  + Reference genome (GRCh38)
    + Format: _FASTA_ (text)
  + Individual contigs
    + Format: _FASTA_ (text)
+ Output:
  + Single file of local alignments
    + Format: _BAM_ (binary)
+ Tools:
  + BWA (exists)
  + **Wrappers** (do not exist)
    + Collect paths (input and output files) (validate?)
    + Call BWA

### Step 1.1: Globally align all contigs
+ Input:
  + Reference genome (GRCh38)
    + Format: _FASTA_ (text)
  + Individual contigs
    + Format: _FASTA_ (text)
  + Local alignment (output of step 1)
    + Format: _BAM_ (binary)
+ Output:
  + Global alignment where each contig has a single alignment
    + Format: _BAM_ (binary)
+ Tools:
  + **Local To Global Alignment** (does not exist)
  + **Wrappers**? (do not exist)

### Step 2: Divide and conquer multiple sequence alignment
+ Input:
  + Global alignment of all contigs (output of step 1.1)
    + Format: _BAM_ (binary)
+ Output:
  + Single global multiple sequence alignment
    + Format: _BAM_ (binary)
+ Tools:
  + **AMC** (does not exist; although some C++ code exists to accomplish part
    of this)
    1. Identify window coordinates (start with windows of size ~10kb) in
    alignment file
    2. Extract sequence data for each window
    3. Convert window to FASTA format
    4. Run MAFFT for each window [could be parallelized]
    5. Reassemble the individual multiple sequence alignments into a single
    alignment file
  + MAFFT
  + **Wrappers** (do not exist)
    + Window size specification?

### Step 3: Create graph genome!
+ Input:
  + Single global multiple sequence alignment (output of step 2)
    + Format: _BAM_ (binary)
+ Output:
  + Graph representation of variation in reference genome + input contigs
    + Format: _VG_, _GFA_?
+ Tools:
  + **BAM to VCF** (does not exist)
  + vg
  + **Wrapper** (does not exist)
    + Convert BAM alignment to VCF
    + Call vg to convert VCF to vg or gfa format
