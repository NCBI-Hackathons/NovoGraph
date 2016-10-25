#!/bin/bash

usage() {
  echo "Usage: $0 [-o <path-to-output-file>] -f [input fasta file] -f [input fasta file 2] "
  exit 1
}

INPUTFILES=();

while getopts ":o:f:" opt; do
  case "${opt}" in
    o)
      OUTFILE=${OPTARG}
      if [ ! -r $OUTFILE ]; then
        echo "Output file $OUTFILE does not exist or it is not readable."
        usage
      fi
      ;;
    f)
	  INPUTFILE=$OPTARG
      if [ ! -r $INPUTFILE ]; then
        echo "Input file $INPUTFILE does not exist or is not readable."
        usage
      fi
      INPUTFILES+=("$INPUTFILE")
      ;;
    *)
      echo -e "\nUnrecognized option: -${OPTARG}"
      usage
      ;;
  esac
done

for file in "${INPUTFILES[@]}"; do
	# ## -op and -ep 
	# mafft --auto --thread 12 --reorder aligned --maxiterate 5 input_file.fa > output_file.fa 

	# mafft --auto --thread 12 --reorder aligned --maxiterate 5 $file >> $outfile
done



