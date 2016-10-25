#!/bin/bash

usage() {
  echo "Usage: $0 [-g <path-to-reference-genome>] [-c <path-to-contigs>]"
  exit 1
}

INFILE=""
CONTIGDIR="."

while getopts ":g:c:" opt; do
  case "${opt}" in
    g)
      INFILE=${OPTARG}
      # do test for file exists
      ;;
    c)
      CONTIGDIR=${OPTARG}
      # do test for directory exists
      ;;
    *)
      echo -e "\nUnrecognized option: ${OPTARG}"
      usage
      ;;
  esac
done

echo "$INFILE"
echo "$CONTIGDIR"

### Step 1: Align all contigs to GRCh38
# + Collect paths (input and output files) (validate?)
# + Call BWA

### Step 1.1: Globally align all contigs
# + Call Local to Global Alignment

### Step 2: Divide and conquer multiple sequence alignment
# + Window size specification?
# + Call AMC

### Step 3: Create graph genome!
# + Call BAM to VCF
# + Call vg to convert VCF to vg or gfa format
