#!/bin/bash

# Wrapper for entire pipeline, largely just calls other scripts for heavy
# lifting
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-10-25

usage() {
  echo "Usage: $0 [-g <path-to-reference-genome>] [-c <path-to-contigs>]"
  exit 1
}

while getopts ":g:c:b:" opt; do
  case "${opt}" in
    g)
      INFILE=${OPTARG}
      if [ ! -r $REFFASTA ]; then
        echo "Input file $REFFASTA does not exist or it is not readable."
        usage
      fi
      ;;
    c)
      CONTIGDIR=${OPTARG}
      if [ ! -r $CONTIGDIR ]; then
        echo "Contig directory $CONTIGDIR does not exist or is not readable."
        usage
      fi
      ;;
    b)
      BAMFILE=${OPTARG}
      if [ ! -r $BAMFILE]; then
        echo "Bath to $BAMFILE does not exist or is not readable."
        usage
      fi
      ;;
    *)
      echo -e "\nUnrecognized option: -${OPTARG}"
      usage
      ;;
  esac
done

### Step 1: Align all contigs to GRCh38
# TODO: This is moving outside of our pipeline
# + Call BWA
LOCALBAMFILE="localbam"
#./gg-01-local.sh -g $REFFASTA -c $CONTIGDIR

### Step 1: Globally align all contigs
# + Call Local to Global Alignment
./gg-01.1-local.sh -g $REFFASTA -c $CONTIGDIR -b $LOCALBAMFILE

### Step 2: Divide and conquer multiple sequence alignment
# Call BAM2MAFFT for window identification, extraction, and FASTA conversion
BAM2MAFFT.pl --referenceFasta $REFFASTA --BAM $LOCALBAMFILE

# MAFFT alignment
# TODO: check the arguments for this path; it might just take two arguments:
# The output destination for FASTA files, and input directory for FASTA
./align_mafft.sh -o $OUTFILE $FASTA1 $FASTA2 ... $FASTAN

# Many FASTA -> Many BAM
# TODO: The perl script fas2bam.pl needs to be called for each window (ideally
# in parallel; could just be added to updated version of align_mafft.sh)
# Almost certainly will be called from within align_mafft.sh (above)

# Many BAM -> single BAM
$FINALBAM="uberbam.bam"
$BAMDIRS=""
# TODO: Insert call to Nancy's code here
# PSEUDOCALL:
concatbam.pl -o $FINALBAM -i $BAMDIRS

### Step 3: Create graph genome!
$VCFFILE="myvcf.vcf"
# BAM to VCF
# TODO: Insert call to Andrew's code here
# PSEUDOCALL:
bam2vcf.pl -b $FINALBAM -o $VCFFILE

# Call vg to convert VCF to vg or gfa format
$VGOUTFILE="vgoutput.vg"
# TODO: What is the small/x.fa argument? a reference genome?
vg construct -r small/x.fa -v $VCFFILE > $VGOUTFILE
