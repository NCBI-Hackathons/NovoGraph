# Graph Genomes
## CSHL Hackathon 24-26 October 2016

## Preliminary description of tools, inputs, and outputs
### Problem: Whole genome graph creation from multiple assemblies

[See Initial presentation Google doc here](https://docs.google.com/a/lbl.gov/presentation/d/17MTjobkF-wfgamiK2NDwoRGQzlEJ1nzw8VtPqNGBZDI/edit?usp=sharing)

[Tuesday presentation Google doc](https://docs.google.com/presentation/d/1sEx0Q0LdAuBQF0t-JJbwZuXHNvDDwxHrevyfoYrxi68/edit?usp=sharing)

## Workflow
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
+ Status:
  + BWA
    + Status: _complete_
  + Wrappers
    + Assignees:
      + Nathan?
      + Jeff
    + Status: in progress
      + See scripts/gg-01-local.sh in GitHub repo

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
+ Status:
  + Local to Global Alignment
    + Assignees:
      + Evan
      + Aarti
    + Status:
      + In progress, see NeedlemanWunsch.ph & NW_alignment_scratch.ph in scripts
      directory on GitHub
  + Wrappers
    + Assignees:
      + Nathan?
      + Jeff
    + Status: in progress
      + See scripts/gg-01.1-local.sh in GitHub repo

### Step 2: Divide and conquer multiple sequence alignment
+ Input:
  + Global alignment of all contigs (output of step 1.1)
    + Format: _BAM_ (binary)
+ Output:
  + Single global multiple sequence alignment
    + Format: _BAM_ (binary)
+ Tools:
  + MAFFT
  + **AMC** (does not exist; although some C++ code exists to accomplish part
    of this)
    1. Identify window coordinates (start with windows of size ~10kb) in
    alignment file
    2. Extract sequence data for each window
    3. Convert window to FASTA format
    4. Run MAFFT for each window [could be parallelized]
    5. Convert MAFFT-produced FASTA alignments to BAM files (one for each
    window)
    6. Reassemble the individual BAM MSAs into a single BAM alignment file
  + **Wrappers** (do not exist)
    + Window size specification?
    + Call AMC
+ Status:
  + MAFFT
    + Status: complete
  + AMC
    1. Identify window coordinates
      + Assignees:
        + Alex
      + Status: in progress
        + BAM2MAFFT.pl in top-level of GitHub repo
    2. Extract sequence data from BAM file
      + Assignees:
        + Alex
      + Status: in progress
        + BAM2MAFFT.pl in top-level of GitHub repo
    3. Convert window to FASTA format
      + Assignees:
        + Alex
      + Status: in progress
        + BAM2MAFFT.pl in top-level of GitHub repo
    4. Run MAFFT
      + Assignees:
        + Nathan
      + Status: in progress
        + See scripts/align_mafft.sh in GitHub repo
    5. Convert each FASTA alignment to BAM
      + Assignees:
        + Nancy
      + Status: _complete_
        + _Location for any code?_
    6. Concatenate each window-BAM to single BAM
      + Assignees:
        + Nancy
      + Status: in progress
        + _Location for any code?_
  + Wrappers:
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
  + **BAM to VCF** (does not exist)
  + **Wrapper** (does not exist)
    + Convert BAM alignment to VCF
    + Call vg to convert VCF to vg or gfa format
+ Status:
  + BAM to VCF:
    + Assignees:
      + Andrew
    + Status: in progress
      + _Location for any code?_
  + Wrapper:
    + Assignees:
      + Jeff
    + Status: **not started**

### Step 4: Generate Docker container et al for distribution
+ Assignees:
  + Nathan
+ Status: in progress
  + See docker-compose.yml & Dockerfile in top-level of GitHub repo