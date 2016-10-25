#!/bin/bash

# Wrapper for calling perl script to convert a FASTA file to BAM
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-10-25

usage() {
  echo "Usage: $0 -f <path-to-fasta-file> [-b <BAM header file>] [-r <reference identifier>]"
  exit 1
}

FASTAFILE=""
BAMHEADER="config/windowbamheader.txt"
REFID="ref"

while getopts ":f:br" opt; do
  case "${opt}" in
    f)
      FASTAFILE=${OPTARG}
      if [ ! -r $FASTAFILE ]; then
        echo "Input file $FASTAFILE does not exist or it is not readable."
        usage
      fi
      ;;
    b)
      BAMHEADER=${OPTARG}
      if [ ! -r $BAMHEADER ]; then
        echo "Contig directory $CONTIGDIR does not exist or is not readable."
        usage
      fi
      ;;
    r)
      REFID=${OPTARG}
      ;;
    *)
      echo -e "\nUnrecognized option: -${OPTARG}"
      usage
      ;;
  esac
done

scripts/fas2bam.pl --fas $FASTAFILE --ref "${REFID}" --bamheader $BAMHEADER
