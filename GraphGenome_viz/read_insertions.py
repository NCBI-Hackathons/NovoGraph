#!/usr/bin/python3

import pysam
import sys

## open file

bamFile = sys.argv[1];

## large insertion parameter value; set by user
## not using argparse because laziness

insertionsize = 200

bamfile = pysam.Samfile(bamFile, "rb")


for read in bamfile:
    if not read.is_unmapped:
        cigarLine = read.cigar
        pretty_big_insertion= False
        for cigarType, cigarLength in cigarLine:
            if cigarType == 1:   ## this is the insertion
                if cigarLength > insertionsize:  
                    pretty_big_insertion = True
                    break
                    
        if pretty_big_insertion:   ### print the positions
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))
