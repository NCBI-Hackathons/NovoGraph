
usage() {
  echo "Usage: $0 [-g <path-to-reference-genome>] [-c <path-to-contigs>]"
  exit 1
}

while getopts ":g:c:" opt; do
  case "${opt}" in
    g)
      INFILE=${OPTARG}
      if [ ! -r $REFFASTA ]; then
        echo "Input file $REFFASTA does not exist or it is not readable."
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

# TODO: Call to BWA to create a BAM file
