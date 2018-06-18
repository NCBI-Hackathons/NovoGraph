#!/usr/bin/python3

## 19, March 2018

import sys
import pysam
import pandas as pd

## using these BAM's with supplementary alignments removed

## samtools view -F 0x800 -b  Korean.contigs.bam > Korean.no800.bam
## samtools view -F 0x800 -b  CHM13.contigs.bam > CHM13.no800.bam
## samtools view -F 0x800 -b  CHM1.contigs.bam > CHM1.no800.bam
## samtools view -F 0x800 -b  HG003.contigs.bam > HG003.no800.bam
## samtools view -F 0x800 -b  HG004.contigs.bam > HG004.no800.bam
## samtools view -F 0x800 -b  HX1.contigs.bam > HX1.no800.bam
## samtools view -F 0x800 -b  NA19240.contigs.bam > NA19240.no800.bam


bam = 'HG003.no800.bam'

bamfile = pysam.Samfile(bam, "rb")


cigar_tups = []
nm_tags = []
aln_lengths = []
    
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
        cigar_tups.append(cigartuple)
        nm_tags.append(read.get_tag("NM"))
        aln_lengths.append(read.query_alignment_length)


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
## = 7 sequence match
## X 8 sequence mismatch
## B 9 "backwards", http://seqanswers.com/forums/showthread.php?t=34440
## cf. http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

df = pd.DataFrame(list(map(group_sum, cigar_tups))).fillna(0)

df = df.rename(index=str, columns={0: "match", 1: "insertion", 2: "deletion", 3: "skipped", 4: "soft_clipping", 5: "hard_clipping", 6: "padding"})

df['NM'] = nm_tags
df['lengths'] = aln_lengths

df['INS_readLength'] = df['insertion']/df['lengths']
df['DEL_readLength'] = df['deletion']/df['lengths']
##  (NM) – (#INs) – (#DELs)
df['NM_INS_DEL'] = df['NM'] - df['insertion'] - df['deletion']


df.to_csv('NM_INS_DEL_HG003.no800.csv', index=False)





