# Graph Genomes
## CSHL Hackathon 24-26 October 2016

## Preliminary description of tools, inputs, and outputs
### Problem: Whole genome graph creation from multiple assemblies

### Step 1: Align all contigs to GRCh38
1. Input:
  + Reference genome (GRCh38)
    + Format: _FASTA_ (text)
  + Individual contigs
    + Format: _FASTA_ (text)
2. Output:
  + Single file of local alignments
    + Format: _BAM_ (binary)
3. Tools:
  + BWA (exists)
  + **Wrappers** (do not exist)
    + Collect paths (input and output files) (validate?)
    + Call BWA