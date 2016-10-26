#!/bin/bash

# Wrapper for entire pipeline, largely just calls other scripts for heavy
# lifting
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-10-25

usage() {
  echo "Usage: $0 [-g <path-to-reference-genome>] [-c <path-to-contigs>] [-s <path-to-list-of-samples>] [-j <number of parallel jobs>]"
  exit 1
}

while getopts ":g:c:s:j:" opt; do
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
    s)
      SAMPLESFILE=${OPTARG}
      if [ ! -r $SAMPLESFILE]; then
        echo "Path to $SAMPLESFILE does not exist or is not readable."
        usage
      fi
      ;;
    j)
      JOBS=${OPTARG}
      if [ ! $JOBS]; then
        echo "jobs param not specified"
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
# TODO: Update with a real Call
# PSEUDOCALL
scripts/NWAM_01.py -g $REFFASTA -c $CONTIGDIR -b $LOCALBAMFILE
#./gg-01.1-local.sh -g $REFFASTA -c $CONTIGDIR -b $LOCALBAMFILE

### Step 2: Divide and conquer multiple sequence alignment
# Call BAM2MAFFT for window identification, extraction, and FASTA conversion
BAM2MAFFT.pl --referenceFasta $REFFASTA --BAM $LOCALBAMFILE

# MAFFT alignment
# Many FASTA -> Many BAM
FASTADIR=""
BAMDIR=""
# TODO: Will ultimately also pass number of processors & number of threads
./gg-align_mafft.sh -i $FASTADIR -o $BAMDIR

# Many BAM -> single BAM
FINALBAM="uberbam.bam"
CONTIGINFO="" # TODO: Formalize this location
globalize_windowbams.pl --fastadir $FASTADIR --msadir $BAMDIR --contigs $CONTIGINFO
# globalize_windowbams.pl --fastadir <path to directory with inputs for MAFFT>
#                         --msadir <path to directory with outputs from MAFFT>
#                         --contigs <path to file with contig names/lengths>

### Step 3: Create graph genome!
VCFFILE="vgInput.vcf.gz"
# BAM to bgzipped VCF in parallel windows of 100000
scripts/gg-wrapper-bam2vcf.pl --script scripts/BAM2VCF.pl --BAM $FINALBAM --referenceFasta $REFFASTA --samples $SAMPLESFILE --window 100000 | parallel -j $JOBS
# concatenate vcf files and index with tabix
scripts/gg-wrapper-vcf-concat.pl --BAM $FINALBAM --referenceFasta $REFFASTA --output $VCFFILE

# Call vg to convert VCF to vg or gfa format
VGOUTFILE="vgoutput.vg"
vg construct -r $REFFASTA -v $VCFFILE > $VGOUTFILE
