#!/usr/bin/python3

import pysam
import sys

## open file

bamFile = sys.argv[1];


## match/insertion/deletion/skipped/soft clipping/hard clipping/padding parameter value
## Set by user
## not using argparse because laziness

matchsize = 55500
insertionsize = 200
deletionsize = 200
skippedsize = 500
softclipping_size = 5550000
hardclipping_size = 5550000
padding_size = 500
sequence_match_size = 2000
sequence_mismatch_size = 2000
backwards_size = 2000

## M 0 alignment match (can be a sequence match or mismatch)
## I 1 insertion to the reference
## D 2 deletion from the reference
## N 3 skipped region from the reference
## S 4 soft clipping (clipped sequences present in SEQ)
## H 5 hard clipping (clipped sequences NOT present in SEQ)
## P 6 padding (silent deletion from padded reference)
## = 7 sequence match
## X 8 sequence mismatch
## B 9 "backwards", http://seqanswers.com/forums/showthread.php?t=34440
## cf. http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

bamfile = pysam.Samfile(bamFile, "rb")


for read in bamfile:
    if not read.is_unmapped:
        cigarLine = read.cigar
        pretty_big_match = False
        pretty_big_insertion = False
        pretty_big_deletion = False
        pretty_big_skipping = False
        pretty_big_softclipping = False
        pretty_big_hardclipping = False
        pretty_big_padding = False
        pretty_big_sequence_match = False
        pretty_big_sequence_mismatch = False
        pretty_big_backwards = False             
        for cigarType, cigarLength in cigarLine:
            if cigarType == 0:   ## alignment match 
                if cigarLength > matchsize:  
                    pretty_big_match = True
                    break
            if cigarType == 1:   ## insertion
                if cigarLength > insertionsize:  
                    pretty_big_insertion = True
                    break
            if cigarType == 2:   ## deletion
                if cigarLength > deletionsize:  
                    pretty_big_deletion = True
                    break
            if cigarType == 3:   ## skipped
                if cigarLength > skippedsize:  
                    pretty_big_skipping = True
                    break
            if cigarType == 4:   ## soft clipping 
                if cigarLength > softclipping_size:  
                    pretty_big_softclipping = True
                    break
            if cigarType == 5:   ## hard clipping 
                if cigarLength > hardclipping_size:  
                    pretty_big_hardclipping = True
                    break      
            if cigarType == 6:   ## padding
                if cigarLength > padding_size:  
                    pretty_big_padding = True
                    break       
            if cigarType == 7:   ## sequence match
                if cigarLength > sequence_match_size:  
                    pretty_big_sequence_match = True
                    break
            if cigarType == 8:   ## sequence mismatch
                if cigarLength > sequence_mismatch_size:  
                    pretty_big_sequence_mismatch = True
                    break      
            if cigarType == 9:   ## backwards
                if cigarLength > backwards_size:  
                    pretty_big_backwards = True
                    break       
        ### print the positions
        if pretty_big_match:   ### M 0 alignment match
            print('MATCH:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_insertion:  ### I 1 insertion 
            print('INSERTION:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_deletion:   ### D 2 deletion
            print('DELETION:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_skipping:   ### N 3 skipped
            print('SKIPPING:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_softclipping:   ### S 4 soft clipping
            print('SOFT_CLIPPING:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_hardclipping:   ### H 5 hard clipping 
            print('HARD_CLIPPING:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_padding:   ### P 6 padding 
            print('PADDING:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_sequence_match:   ### = 7 sequence match 
            print('SEQUENCE MATCH:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if sequence_mismatch_size:   ### X 8 sequence mismatch
            print('SEQUENCE MISMATCH:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))

        if pretty_big_backwards:   ### B 9 "backwards" 
            print('BACKWARDS:')
            print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))





