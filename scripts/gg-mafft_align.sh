#!/bin/bash

# TODO: pass in number of procs and threads
usage() {
  echo "Usage: $0 [ -i <path to input files>  -o <path to output file>  ]  "
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
INPUTFILES=($( grep -r -c "^>" $INPUT_DIRECTORY | grep -v ":1"  | awk -F':' '{ print $1 }' | grep ".fa"  ))

echo "Reference Only"
FILES_TO_IGNORE=($( grep -r -c "^>" $INPUT_DIRECTORY | grep ":1"  | awk -F':' '{ print $1 }'| grep ".fa"   ))


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
	echo "Copied $i files."
}


function align(){

 echo $1 $2 $3

  fullfile=$1
  filename=$(basename "$fullfile");
  extension="${filename##*.}";
  filename="${filename%.*}";

  outputfile="$3/${filename}_aligned.fa";

  echo "from $1 to $outputfile"

  mafft --auto  $1 > $outputfile

#  FILESIZE=$(du -sb $outputfile| awk '{ print $1 }')
#  if (($FILESIZE==0)) ; then
#     cp -f $1 $outputfile
#  fi

  if [ $? -eq 1 ]; then
     cp -f $1 $outputfile
  fi

#  gg-wrapper-fas2bam
  ./fas2bam.pl --fas $outputfile --ref "ref" --bamheader "config/windowbamheader.txt"

}

export -f align


copyEmpty

parallel -j 8  align ::: ${INPUTFILES[@]} ::: $INPUT_DIRECTORY ::: $OUTPUT_DIRECTORY





