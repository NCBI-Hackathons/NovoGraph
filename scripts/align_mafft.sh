#!/bin/bash

usage() {
  echo "Usage: $0 [-o <path-to-output-file>] -f [input fasta file] -f [input fasta file 2] "
  exit 1
}


while getopts ":o:i:" opt; do
  case "${opt}" in
    o)
      OUTFILE=${OPTARG}
      ;;
    i)
	  INPUTDIRECTORY=$OPTARG
      if [ ! -r $INPUTDIRECTORY ]; then
        echo "Input file $INPUTDIRECTORY does not exist or is not readable."
        usage
      fi
      ;;
    *)
      echo -e "\nUnrecognized option: -${OPTARG}"
      usage
      ;;
  esac
done




INPUTFILES=`ls $INPUTDIRECTORY`


for fullfile in "${INPUTFILES[@]}"; do
	filename=$(basename "$fullfile")
	extension="${filename##*.}"
	filename="${filename%.*}"

	outputfile="${filename}_aligned.fa"
#    echo "$filename from $fullfile" 
	# ## -op and -ep 
	# mafft --auto --thread 12 --reorder aligned --maxiterate 5 input_file.fa > output_file.fa 
#   echo $file
	mafft --globalpair --thread 1 --reorder aligned --maxiterate 10 $fullfile > $outfile
done



