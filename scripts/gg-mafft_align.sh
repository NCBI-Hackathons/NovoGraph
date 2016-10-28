#!/bin/bash

# TODO: pass in number of procs and threads
usage() {
  echo "Usage: $0 [ -i <path to input files to grep>  -o <path to output file>  ]  "
  exit 1
}


while getopts ":o:i:" opt; do
  case "${opt}" in
     o)
      OUTPUT_DIRECTORY=$OPTARG
      if [ ! -r $OUTPUT_DIRECTORY ]; then
        echo "Output file $OUTPUT_DIRECTORY does not exist or is not readable."
        usage
      fi
      ;;
    i)
      INPUT_DIRECTORY=$OPTARG
#      if [ ! -r $INPUT_DIRECTORY ]; then
#        echo "Input file $INPUT_DIRECTORY does not exist or is not readable."
#        usage
#      fi
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


echo "INPUT DIRECTORY: $INPUT_DIRECTORY"
echo "OUTPUT DIRECTORY: $OUTPUT_DIRECTORY"

INPUTFILES=($( grep -r -c "^>" $INPUT_DIRECTORY | grep -v ":1$"  | awk -F':' '{ print $1 }' ))


echo "Files to process: ${#INPUTFILES[@]}"

function align(){

 echo input_file $1 input_directory $2 output_directory $3 and 4 $4

  fullfile=$1
  filename=$(basename "$fullfile");
  extension="${filename##*.}";
  filename="${filename%.*}";

  outputfile="$4/${filename}_aligned.fas";

  echo "from $1 to $outputfile"

  mafft --thread -1 --reorder --auto  $1 > $outputfile

  bamoutput="$4/$filename.bam"

  echo "BAM input $outputfile and output $bamoutput"

  ./scripts/fas2bam.pl --input $outputfile --output $bamoutput --ref "ref" --bamheader "./config/windowbam.header.txt" 

}

export -f align

# copyEmpty

# parallel -j 8  align ::: ${INPUTFILES[@]} ::: $INPUT_DIRECTORY ::: $OUTPUT_DIRECTORY

for INPUT_FILE in ${INPUTFILES[@]}
do 
  align $INPUT_FILE $INPUT_DIRECTORY $OUTPUT_DIRECTORY ;
done 






