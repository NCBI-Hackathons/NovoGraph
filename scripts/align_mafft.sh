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

echo "Files to process"
INPUTFILES=($( grep -r -c "^>" $INPUT_DIRECTORY | grep -v ":1" | awk -F':' '{ print $1 }' ))

echo "Reference Only"
FILES_TO_IGNORE=($( grep -r -c "^>" $INPUT_DIRECTORY | grep  ":1" | awk -F':' '{ print $1 }' ))


echo "Files to process: ${#INPUTFILES[@]}"
echo "Files to ignore: ${#FILES_TO_IGNORE[@]}"


function copyEmpty(){
	i=0
	for line in "${FILES_TO_IGNORE[@]}"
	do
	  ((++i))
	  filename=$(basename "$line");
	  extension="${filename##*.}";
	  filename="${filename%.*}";

	  outputfile="${OUTPUT_DIRECTORY}/${filename}_aligned.fa";
#      echo "copying $line to $outputfile"
	  cp -f $line $outputfile
	done 
}

copyEmpty 

echo "Copied $i files."

function align(){
  fullfile=$1 
  filename=$(basename "$fullfile");
  extension="${filename##*.}";
  filename="${filename%.*}";

  fullinputfile="${INPUT_DIRECTORY}/$fullfile"  
  outputfile="${OUTPUT_DIRECTORY}/${filename}_aligned.fa";
  mafft --auto  ${fullinputfile} > $outputfile

}

#for fullfile in "${INPUTFILES[@]}"; do
#  parallel -j 4 align {} ::: $fullfile
#done





