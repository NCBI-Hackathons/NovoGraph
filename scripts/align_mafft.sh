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

if [ -z ${INPUTDIRECTORY+x} ]; then 
	echo "INPUTDIRECTORY is unset"; 
	exit 1; 
fi

if [ -z ${OUTFILE+x} ]; then 
	echo "OUTFILE is unset"; 
	exit 1; 
fi


i=0
while read line
do
    INPUTFILES[ $i ]="$line"        
    (( i++ ))
done < <(ls $INPUTDIRECTORY)


for fullfile in "${INPUTFILES[@]}"; do
  filename=$(basename "$fullfile");
  extension="${filename##*.}";
  filename="${filename%.*}";
  
  outputfile="${INPUTDIRECTORY}/${filename}_aligned.fa";
  mafft --auto --maxiterate 10 ${INPUTDIRECTORY}/$fullfile > $outputfile
done





