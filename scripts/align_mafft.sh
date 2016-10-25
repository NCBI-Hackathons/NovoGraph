#!/bin/bash

usage() {
  echo "Usage: $0 [-o <path-to-output-file>] -f [input fasta file] -f [input fasta file 2] "
  exit 1
}


while getopts ":o:i:" opt; do
  case "${opt}" in
    o)
      OUTPUT_DIRECTORY=${OPTARG}
      ;;
    i)
	  INPUT_DIRECTORY=$OPTARG
      if [ ! -r $INPUT_DIRECTORY ]; then
        echo "Input file $INPUT_DIRECTORY does not exist or is not readable."
        usage
      fi
      ;;
    *)
      echo -e "\nUnrecognized option: -${OPTARG}"
      usage
      ;;
  esac
done

if [ -z ${INPUT_DIRECTORY+x} ]; then 
	echo "INPUT_DIRECTORY is unset"; 
	exit 1; 
fi

if [ -z ${OUTPUT_DIRECTORY+x} ]; then 
	echo "OUTPUT_DIRECTORY is unset"; 
	exit 1; 
fi


i=0
while read line
do
    INPUTFILES[ $i ]="$line"        
    (( i++ ))
done < <(ls $INPUT_DIRECTORY)


for fullfile in "${INPUTFILES[@]}"; do
  filename=$(basename "$fullfile");
  extension="${filename##*.}";
  filename="${filename%.*}";

  fullinputfile="${INPUT_DIRECTORY}/$fullfile"  
  outputfile="${OUTPUT_DIRECTORY}/${filename}_aligned.fa";
  mafft --auto  ${fullinputfile} > $outputfile

  FILESIZE=$(du -sb $outputfile| awk '{ print $1 }')

  if (($FILESIZE==0)) ; then 
     cp -f $fullinputfile $outputfile
  fi 
done





