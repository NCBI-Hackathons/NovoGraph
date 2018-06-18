#!/usr/bin/python3


import pysam
import pandas as pd



## using the BAM's with supplementary alignments removed

## samtools view -F 0x800 -b  Korean.contigs.bam > Korean.no800.bam
## samtools view -F 0x800 -b  CHM13.contigs.bam > CHM13.no800.bam
## samtools view -F 0x800 -b  CHM1.contigs.bam > CHM1.no800.bam
## samtools view -F 0x800 -b  HG003.contigs.bam > HG003.no800.bam
## samtools view -F 0x800 -b  HG004.contigs.bam > HG004.no800.bam
## samtools view -F 0x800 -b  HX1.contigs.bam > HX1.no800.bam
## samtools view -F 0x800 -b  NA19240.contigs.bam > NA19240.no800.bam


## CIGAR key
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



bam = 'CHM13.no800.bam'

bamfile = pysam.AlignmentFile(bam, "rb")

### insertions

nm_tags = []
lengths_tags =[] 

for read in bamfile:
    if not read.is_unmapped:
        cigarLine = read.cigar           
        for cigarType, cigarLength in cigarLine:
            if cigarType == 1:   ## insertion 
                nm_tags.append(read.get_tag("NM"))
                lengths_tags.append(read.query_alignment_length)



df = pd.DataFrame()
df['NM'] = nm_tags
df['length'] = lengths_tags

df['ratio'] = df['NM']/df['length']

df.to_csv('insertions_CHM13.no800.csv', index=False)


### deletions


bam = 'CHM13.no800.bam'

bamfile = pysam.AlignmentFile(bam, "rb")


nm_tags = []
lengths_tags =[] 

for read in bamfile:
    if not read.is_unmapped:
        cigarLine = read.cigar           
        for cigarType, cigarLength in cigarLine:
            if cigarType == 2:   ## deletions
                nm_tags.append(read.get_tag("NM"))
                lengths_tags.append(read.query_alignment_length)



df = pd.DataFrame()
df['NM'] = nm_tags
df['length'] = lengths_tags

df['ratio'] = df['NM']/df['length']

df.to_csv('deletions_CHM13.no800.csv', index=False)


### mismatch


bam = 'CHM13.no800.bam'

bamfile = pysam.AlignmentFile(bam, "rb")


nm_tags = []
lengths_tags =[] 

for read in bamfile:
    if not read.is_unmapped:
        cigarLine = read.cigar           
        for cigarType, cigarLength in cigarLine:
            if cigarType == 8:   ## mismatch
                nm_tags.append(read.get_tag("NM"))
                lengths_tags.append(read.query_alignment_length)



df = pd.DataFrame()
df['NM'] = nm_tags
df['length'] = lengths_tags

df['ratio'] = df['NM']/df['length']

df.to_csv('mismatch_CHM13.no800.csv', index=False)








