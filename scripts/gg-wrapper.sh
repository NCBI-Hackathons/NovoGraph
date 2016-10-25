#!/bin/bash

usage() {
  echo "Usage: $0 [-g <path-to-reference-genome>] [-c <path-to-contigs>]"
  exit 1
}

while getopts ":g:c:" opt; do
  case "${opt}" in
    g)
      INFILE=${OPTARG}
      if [ ! -r $INFILE ]; then
        echo "Input file $INFILE does not exist or it is not readable."
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
    *)
      echo -e "\nUnrecognized option: -${OPTARG}"
      usage
      ;;
  esac
done

### Step 1: Align all contigs to GRCh38
# + Call BWA
# ./gg-01-local.sh -g $INFILE -c $CONTIGDIR

### Step 1.1: Globally align all contigs
# + Call Local to Global Alignment
# ./gg-01.1-local.sh -g $INFILE -c $CONTIGDIR -b <path-to-bamfile>

### Step 2: Divide and conquer multiple sequence alignment
# + Window size specification?
# + Call AMC
# + MAFFT alignment

### Step 3: Create graph genome!
# + Call BAM to VCF
# + Call vg to convert VCF to vg or gfa format
