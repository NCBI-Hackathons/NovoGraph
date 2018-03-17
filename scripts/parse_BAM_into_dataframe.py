#!/usr/bin/python3

import sys
import pysam
import pandas as pd

## name output *csv file
output_filename = 'output_parsedBAM.csv'


## open file
bam = sys.argv[1];

bamfile = pysam.Samfile(bam, "rb")


locations = []
cigar_tups = []
    
## function to parse CIGAR tuples from pysam
from itertools import groupby
def group_sum(input):
    """input list, returns dictionary of summation based on initial key"""
    result_dict = {k: sum(v_[1] for v_ in v) for k, v in groupby(sorted(input, key=lambda x: x[0]), lambda x: x[0])}
    return result_dict


## standard pysam for loop
for read in bamfile:
    if not read.is_unmapped:
        cigartuple = read.cigar
        chrom_start_end = "{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end)
        cigar_tups.append(cigartuple)
        locations.append(chrom_start_end)


data = [{f[0]:[] for f in e} for e in cigar_tups]

## add values to the dict as list which can capture multiple values
## creates lots of noise
[[data[k][e[0]].append(e[1]) for e in v] for k,v in enumerate(cigar_tups)]

#sum values for each key for each row.
cigar_dict = [{k:sum(v) for k,v in e.items()} for e in data]

## M 0 alignment match (can be a sequence match or mismatch)
## I 1 insertion to the reference
## D 2 deletion from the reference
## N 3 skipped region from the reference
## S 4 soft clipping (clipped sequences present in SEQ)
## H 5 hard clipping (clipped sequences NOT present in SEQ)
## P 6 padding (silent deletion from padded reference)
## cf. http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

df = pd.DataFrame(list(map(group_sum, cigar_tups))).fillna(0)

df['position'] = locations

## renaming
df = df.rename(index=str, columns={0: "match", 1: "insertion", 2: "deletion", 3: "skipped", 4: "soft_clipping", 5: "hard_clipping", 6: "padding"})

## re-order columns
df = df[['position', 'match', 'insertion', 'deletion', 'soft_clipping', 'hard_clipping', 'padding']]

### save to csv
df.to_csv(output_filename, index=False)

