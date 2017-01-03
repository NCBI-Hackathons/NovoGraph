##############################################
##
## Simple script to create test contigs:
##     1) 10kb chrom1 with dozen SNPs
##     2) 10kb chrom1 with single 300bp deletion
##     3) 10kb chrom1 with single 200bp insertion of randomly created bp
##     4) 10kb chrom1 with SV:
##         a) 150bp, 500bp, 1000bp inversions
##         b) tandem duplications:
##             i)   150bp = 3 x 50bp
##             ii)  500bp = 2 x 250bp
##             iii) 1000bp = 10 x 100bp
##     5) three contigs: one of each of the above, non-overlapping
##     6) three contigs: one of each of the above, overlapping
##
##
## source("https://bioconductor.org/biocLite.R")
## biocLite("BSgenome")
##
## first install hg38
## source("https://bioconductor.org/biocLite.R")
## biocLite("BSgenome.Hsapiens.UCSC.hg38")
##
## I run each of the outputs below in http://www.ebi.ac.uk/Tools/sfc/emboss_seqret/
## It's probably unnecessary
## set of DNA sequences
## INPUT FORMAT: Unknown format; OUTPUT FORMAT: FASTA format
##
library("BSgenome")
human_ref_genome = getBSgenome('BSgenome.Hsapiens.UCSC.hg38', masked=FALSE)
## selection chrom1
human_ref_genome$chr1
chrom1 = human_ref_genome$chr1
section_chr1_10kb = chrom1[6000000:6009999] # whimsically chose this interval on chrom1
print(length(section_chr1_10kb)) # prints 10000
## my SNPs generally follow this scheme:
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
print(section_chr1_10kb[4547])  # C -> T
print(section_chr1_10kb[8326])  # T -> C
print(section_chr1_10kb[4602])  # C -> A
print(section_chr1_10kb[6171])  # A -> G
print(section_chr1_10kb[8328])  # A -> G
print(section_chr1_10kb[6394])  # C -> T
print(section_chr1_10kb[8399])  # A -> C
print(section_chr1_10kb[9088])  # C -> A
print(section_chr1_10kb[9569])  # T -> G
print(section_chr1_10kb[1313])  # A -> G
print(section_chr1_10kb[4387])  # G -> A
print(section_chr1_10kb[6813])  # T -> C
## create 12 SNPs on this section; chose these positions randomly
section_chr1_10kb[4547] = "T"
section_chr1_10kb[8326] = "C"
section_chr1_10kb[4602] = "A"
section_chr1_10kb[6171] = "G"
section_chr1_10kb[8328] = "G"
section_chr1_10kb[6394] = "T"
section_chr1_10kb[8399] = "C"
section_chr1_10kb[9088] = "A"
section_chr1_10kb[9569] = "G"
section_chr1_10kb[1313] = "G"
section_chr1_10kb[4387] = "A"
section_chr1_10kb[6813] = "C"
export(section_chr1_10kb, "dozenSNPS_chr1_hg38.fasta.gz", compress=TRUE)
## reset these values
section_chr1_10kb = chrom1[6000000:6009999]
## chose deletion position randomly, e.g. # deletion_position <- sample(1:10001, 1)
## print(deletion_position)     7419
## therefore, choose 7400 to 7700
single_deletion = c(section_chr1_10kb[1:7400], section_chr1_10kb[7701:10000])
print(length(single_deletion))
deleted = section_chr1_10kb[7401:7700]  # double check
print(deleted)
##  300-letter "DNAString" instance
## seq: CTGGGATTTAAGAAGTTAAAACCCACCCCCTCTGCC...GAGACTCAGCATTTCCTTTATCTCCTTTTACTGTGG
print(section_chr1_10kb[7390:7400])
print(section_chr1_10kb[7390:7402]) # cut out 'CTGGGATT...'
print(section_chr1_10kb[7698:7702]) # cut out '..AT'
## 11-letter "DNAString" instance, seq: CTTCTCCCATC
## 13-letter "DNAString" instance, seq: CTTCTCCCATCCT
## 5-letter "DNAString" instance, seq: TGGAT
export(single_deletion, "deletion_chr1_hg38.fasta.gz", compress=TRUE)
##
## Did this in Python:
## ''.join(random.choice("ACTG") for i in range(200))
##
## GCAATACGCAATCACCAATCCGAGGTAATCACCCTGGGCCCGGTCCGGGGTCATAGCATCCGTGAGG
##  TTGGATTGGGTAACAGAGTGCGGGGTGGCATTATTGTATGGGAACGAACCTTCAATTCATCTGGTAGGTAAGA
##  CGATGTACTTAAACCATGGCGATGGGCTCTACAATAGGCGACTGGTTACGGGGAGGGTCA
##
insertion = "GCAATACGCAATCACCAATCCGAGGTAATCACCCTGGGCCCGGTCCGGGGTCATAGCATCCGTGAGGTTGGATTGGGTAACAGAGTGCGGGGTGGCATTATTGTATGGGAACGAACCTTCAATTCATCTGGTAGGTAAGACGATGTACTTAAACCATGGCGATGGGCTCTACAATAGGCGACTGGTTACGGGGAGGGTCA"

## insertion_position <- sample(1:10001, 1)
## print(insertion_position) # 2502 ---> choose 2500
section_chr1_10kb = chrom1[6000000:6009999]  # reset
library(Biostrings)
inserted_segment = DNAString(insertion) # create correct data structure
insertion_section = c(section_chr1_10kb[1:2500], inserted_segment, section_chr1_10kb[2501:10000])
print(length(inserted_segment)) # 200 bp
print(section_chr1_10kb[2495:2500])
print(section_chr1_10kb[2495:2503])
print(insertion_section[2495:2503])
## 6-letter "DNAString" instance, seq: CAGACT
## 9-letter "DNAString" instance, seq: CAGACTCCC
## 9-letter "DNAString" instance, seq: CAGACTGCA
export(insertion_section, "insertion_chr1_hg38.fasta.gz", compress=TRUE)
## sv_position <- sample(1:10001, 1)
## print(sv_position)   # 4519 --> 4500
## inversion/flipped SV
section_chr1_10kb = chrom1[6000000:6009999] # reset
## SV: 150 base pair
## 150 bp inversion
## reverse(section_chr1_10kb[4501:4650])
inversion150 = c(section_chr1_10kb[1:4500], reverse(section_chr1_10kb[4501:4650]), section_chr1_10kb[4651:10000])
export(inversion150, "inversion150_chr1_hg38.fasta.gz", compress=TRUE)
## 3100
inversion500 = c(section_chr1_10kb[1:3100], reverse(section_chr1_10kb[3101:3600]), section_chr1_10kb[3601:10000])
export(inversion500, "inversion500_chr1_hg38.fasta.gz", compress=TRUE)
## 5900
inversion1000 = c(section_chr1_10kb[1:5900], reverse(section_chr1_10kb[5901:6900]), section_chr1_10kb[6901:10000])
export(inversion1000, "inversion1000_chr1_hg38.fasta.gz", compress=TRUE)
## also tandem duplications at 150bp, 500bp, 1000bp
section_chr1_10kb = chrom1[6000000:6009999]
## 1700   50bp /times 3
tandem_repeat150 = c(section_chr1_10kb[1:1700], section_chr1_10kb[1701:1750], section_chr1_10kb[1701:1750], section_chr1_10kb[1701:1750], section_chr1_10kb[1851:10000])
export(tandem_repeat150, "tandem_repeat150_chr1_hg38.fasta.gz", compress=TRUE) # 3 x 50 bp
## 8500  250 bp /times 2
tandem_repeat500 = c(section_chr1_10kb[1:8500], section_chr1_10kb[8501:8750], section_chr1_10kb[8501:8750], section_chr1_10kb[9001:10000])
export(tandem_repeat500, "tandem_repeat500_chr1_hg38.fasta.gz", compress=TRUE) # 10 x 100 bp
# 3900 --  10 100 bp
tandem_repeat1000 = c(section_chr1_10kb[1:3900], section_chr1_10kb[3901:4000], section_chr1_10kb[3901:4000],
section_chr1_10kb[3901:4000], section_chr1_10kb[3901:4000], section_chr1_10kb[3901:4000],
section_chr1_10kb[3901:4000], section_chr1_10kb[3901:4000], section_chr1_10kb[3901:4000],
section_chr1_10kb[3901:4000], section_chr1_10kb[3901:4000], section_chr1_10kb[4901:10000])
export(tandem_repeat1000, "tandem_repeat1000_chr1_hg38.fasta.gz", compress=TRUE) # 10 x 100 bp
##
## all rearrangements within one contig, non-overlapping (in principle...for the most part...laziness)
##
section_chr1_10kb = chrom1[6000000:6009999]
## 12 SNPs
section_chr1_10kb[4547] = "T"
section_chr1_10kb[8326] = "C"
section_chr1_10kb[4602] = "A"
section_chr1_10kb[6171] = "G"
section_chr1_10kb[8328] = "G"
section_chr1_10kb[6394] = "T"
section_chr1_10kb[8399] = "C"
section_chr1_10kb[9088] = "A"
section_chr1_10kb[9569] = "G"
section_chr1_10kb[1313] = "G"
section_chr1_10kb[4387] = "A"
section_chr1_10kb[6813] = "C"
## insertion 200 bp
##
## Did this in Python:
## import random
## ''.join(random.choice("ACTG") for i in range(200))
##
## TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGT
##  CTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCC
##  GTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA
##
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))  # 10200
## deletion 300 bp, at 2000 bp
deleted = c(inserted[1:2000], inserted[2301:10200])
print(length(deleted))   # 9900
## SV, inversion 150 bp, at 7000 bp
inverted = c(deleted[1:7000], reverse(deleted[7001:7150]), deleted[7151:9000])
print(length(inverted))  # 9900
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp
final = c(inverted[1:5000], inverted[5001:5050], inverted[5001:5050], inverted[5001:5051], inverted[5151:9000])
export(final, "nonoverlap_chr1_hg38.fasta.gz", compress=TRUE)
##
## same thing, different chromosome section, chrom1, non-overlapping
##
section_chr1_10kb = chrom1[4000000:4009999]
## 12 SNPs
print(section_chr1_10kb[4547])
print(section_chr1_10kb[8326])
print(section_chr1_10kb[4602])
print(section_chr1_10kb[6171])
print(section_chr1_10kb[8328])
print(section_chr1_10kb[6394])
print(section_chr1_10kb[8399])
print(section_chr1_10kb[9088])
print(section_chr1_10kb[9569])
print(section_chr1_10kb[1313])
print(section_chr1_10kb[4387])
print(section_chr1_10kb[6813])
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## A T T T C T A C T T G C
## G G G G T C C A C C T T
section_chr1_10kb = chrom1[4000000:4009999]
## 12 SNPs
## G G G G T C C A C C T T
section_chr1_10kb[4547] = "G"
section_chr1_10kb[8326] = "G"
section_chr1_10kb[4602] = "G"
section_chr1_10kb[6171] = "G"
section_chr1_10kb[8328] = "T"
section_chr1_10kb[6394] = "C"
section_chr1_10kb[8399] = "C"
section_chr1_10kb[9088] = "A"
section_chr1_10kb[9569] = "C"
section_chr1_10kb[1313] = "C"
section_chr1_10kb[4387] = "T"
section_chr1_10kb[6813] = "T"
## insertion 200 bp
##
## Did this in Python:
## import random
## ''.join(random.choice("ACTG") for i in range(200))
##
## TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGT
##  CTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTA
##  CCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA
#
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 2000 bp
deleted = c(inserted[1:2000], inserted[2301:10200])
print(length(deleted))
## SV, inversion 150 bp, at 7000 bp
inverted = c(deleted[1:7000], reverse(deleted[7001:7150]), deleted[7151:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp
final = c(inverted[1:5000], inverted[5001:5050], inverted[5001:5050], inverted[5001:5051], inverted[5151:9000])
export(final, "nonoverlap2_chr1_hg38.fasta.gz", compress=TRUE)
##
## same thing, different chromosome section, chrom1, part 3, non-overlapping
##
chrom1 = human_ref_genome$chr1
section_chr1_10kb = chrom1[8000000:8009999]
print(section_chr1_10kb[2326])
print(section_chr1_10kb[2602])
print(section_chr1_10kb[3171])
print(section_chr1_10kb[3328])
print(section_chr1_10kb[3394])
print(section_chr1_10kb[3399])
print(section_chr1_10kb[4088])
print(section_chr1_10kb[4569])
print(section_chr1_10kb[5313])
print(section_chr1_10kb[5387])
print(section_chr1_10kb[5547])
print(section_chr1_10kb[5813])
## C G G G T A A G A A C A
chrom1 = human_ref_genome$chr1
section_chr1_10kb = chrom1[8000000:8009999]
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## C G G G T A A G A A C A
section_chr1_10kb[2326] = "A"  # C -> A
section_chr1_10kb[2602] = "A"  # G -> A
section_chr1_10kb[3171] = "A"  # G -> A
section_chr1_10kb[3328] = "A"  # G -> A
section_chr1_10kb[3394] = "C"  # T -> C
section_chr1_10kb[3399] = "C"  # A -> C
section_chr1_10kb[4088] = "C"  # A -> C
section_chr1_10kb[4569] = "A"  # G -> A
section_chr1_10kb[5313] = "C"  # A -> C
section_chr1_10kb[5387] = "C"  # A -> C
section_chr1_10kb[5547] = "A"  # C -> A
section_chr1_10kb[5813] = "C"  # A -> C
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp  [3000:3200] ---> total length 10200
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 3100 bp  [3100:2800] --> total length 9000
deleted = c(inserted[1:2800], inserted[3101:10200])
print(length(deleted))
## SV, inversion 150 bp, at 3050 bp  [3050:3200]
inverted = c(deleted[1:3050], reverse(deleted[3050:3200]), deleted[3201:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp [3150:3300]
final = c(inverted[1:3150], inverted[3151:3200], inverted[3151:3200], inverted[3151:3200], inverted[3301:9000])
export(final, "overlap3_chr1_hg38.fasta.gz", compress=TRUE)
##
## same idea, non-overlapping, different chromosome section, chrom2
##
chrom2 = human_ref_genome$chr2
section_chr1_10kb = chrom2[6000000:6009999]
## 12 SNPs
## A C G G C C C C A T A T
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
section_chr1_10kb[4547] = "G"  # A
section_chr1_10kb[8326] = "T"  # C
section_chr1_10kb[4602] = "T"  # G
section_chr1_10kb[6171] = "T"  # G
section_chr1_10kb[8328] = "T"  # C
section_chr1_10kb[6394] = "T"  # C
section_chr1_10kb[8399] = "A"  # C
section_chr1_10kb[9088] = "A"  # C
section_chr1_10kb[9569] = "C"  # A
section_chr1_10kb[1313] = "C"  # T
section_chr1_10kb[4387] = "C"  # A
section_chr1_10kb[6813] = "C"  # T
## insertion 200 bp
##
## Did this in Python:
## import random
## ''.join(random.choice("ACTG") for i in range(200))
##
## TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGT
##  CTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGG
##  ACCCTAACTCTCCTGCTTTAAAAAAGA
##
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 2000 bp
deleted = c(inserted[1:2000], inserted[2301:10200])
print(length(deleted))
## SV, inversion 150 bp, at 7000 bp
inverted = c(deleted[1:7000], reverse(deleted[7001:7150]), deleted[7151:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp
final = c(inverted[1:5000], inverted[5001:5050], inverted[5001:5050], inverted[5001:5051], inverted[5151:9000])
export(final, "nonoverlap_chr2_hg38.fasta.gz", compress=TRUE)
##
## same idea, non-overlapping, different chromosome section, chrom3
##
chrom3 = human_ref_genome$chr3
section_chr1_10kb = chrom3[6000000:6009999]
## 12 SNPs
## A T T G T C T A A T A A
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
section_chr1_10kb[4547] = "G" # A
section_chr1_10kb[8326] = "G" # T
section_chr1_10kb[4602] = "G" # T
section_chr1_10kb[6171] = "T" # G
section_chr1_10kb[8328] = "C" # T
section_chr1_10kb[6394] = "A" # C
section_chr1_10kb[8399] = "G" # T
section_chr1_10kb[9088] = "C" # A
section_chr1_10kb[9569] = "G" # A
section_chr1_10kb[1313] = "C" # T
section_chr1_10kb[4387] = "T" # A
section_chr1_10kb[6813] = "G" # A
## insertion 200 bp
##
## Did this in Python:
## import random
## ''.join(random.choice("ACTG") for i in range(200))
##
## TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGT
##  CTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGG
##  GACCCTAACTCTCCTGCTTTAAAAAAGA
##
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 2000 bp
deleted = c(inserted[1:2000], inserted[2301:10200])
print(length(deleted))
## SV, inversion 150 bp, at 7000 bp
inverted = c(deleted[1:7000], reverse(deleted[7001:7150]), deleted[7151:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp
final = c(inverted[1:5000], inverted[5001:5050], inverted[5001:5050], inverted[5001:5051], inverted[5151:9000])
export(final, "nonoverlap_chr3_hg38.fasta.gz", compress=TRUE)
##
##
##
## Now, overlapping
##
chrom1 = human_ref_genome$chr1
section_chr1_10kb = chrom1[6000000:6009999]
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## C A T C A C G T A G C C
section_chr1_10kb[2326] = "A"  # C -> A
section_chr1_10kb[2602] = "C"  # A -> C
section_chr1_10kb[3171] = "G"  # T -> G
section_chr1_10kb[3328] = "T"  # C -> T
section_chr1_10kb[3394] = "G"  # A -> G
section_chr1_10kb[3399] = "A"  # C -> A
section_chr1_10kb[4088] = "T"  # G -> T
section_chr1_10kb[4569] = "C"  # T -> C
section_chr1_10kb[5313] = "G"  # A -> G
section_chr1_10kb[5387] = "T"  # G -> T
section_chr1_10kb[5547] = "A"  # C -> A
section_chr1_10kb[5813] = "A"  # C -> A
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp  [3000:3200] ---> total length 10200
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 3100 bp  [3100:2800] --> total length 9000
deleted = c(inserted[1:2800], inserted[3101:10200])
print(length(deleted))
## SV, inversion 150 bp, at 3050 bp  [3050:3200]
inverted = c(deleted[1:3050], reverse(deleted[3050:3200]), deleted[3201:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp [3150:3300]
final = c(inverted[1:3150], inverted[3151:3200], inverted[3151:3200], inverted[3151:3200], inverted[3301:9000])
export(final, "overlap_chr1_hg38.fasta.gz", compress=TRUE)
##
## overlapping, different section of chrom1
##
chrom1 = human_ref_genome$chr1
section_chr1_10kb = chrom1[4000000:4009999]
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## G G G T G G T G G C A A
section_chr1_10kb[2326] = "T"  # G -> T
section_chr1_10kb[2602] = "T"  # G -> T
section_chr1_10kb[3171] = "T"  # G -> T
section_chr1_10kb[3328] = "G"  # T -> G
section_chr1_10kb[3394] = "T"  # G -> T
section_chr1_10kb[3399] = "T"  # G -> T
section_chr1_10kb[4088] = "G"  # T -> G
section_chr1_10kb[4569] = "T"  # G -> T
section_chr1_10kb[5313] = "T"  # G -> T
section_chr1_10kb[5387] = "T"  # C -> T
section_chr1_10kb[5547] = "G"  # A -> G
section_chr1_10kb[5813] = "G"  # A -> G
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp  [3000:3200] ---> total length 10200
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 3100 bp  [3100:2800] --> total length 9000
deleted = c(inserted[1:2800], inserted[3101:10200])
print(length(deleted))
## SV, inversion 150 bp, at 3050 bp  [3050:3200]
inverted = c(deleted[1:3050], reverse(deleted[3050:3200]), deleted[3201:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp [3150:3300]
final = c(inverted[1:3150], inverted[3151:3200], inverted[3151:3200], inverted[3151:3200], inverted[3301:9000])
export(final, "overlap2_chr1_hg38.fasta.gz", compress=TRUE)
##
## overlapping, different section of chrom1, part 3
##
chrom1 = human_ref_genome$chr1
section_chr1_10kb = chrom1[8000000:8009999]
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## C G G G T A A G A A C A
section_chr1_10kb[2326] = "A"  # C -> A
section_chr1_10kb[2602] = "A"  # G -> A
section_chr1_10kb[3171] = "A"  # G -> A
section_chr1_10kb[3328] = "A"  # G -> A
section_chr1_10kb[3394] = "C"  # T -> C
section_chr1_10kb[3399] = "C"  # A -> C
section_chr1_10kb[4088] = "C"  # A -> C
section_chr1_10kb[4569] = "A"  # G -> A
section_chr1_10kb[5313] = "C"  # A -> C
section_chr1_10kb[5387] = "C"  # A -> C
section_chr1_10kb[5547] = "A"  # C -> A
section_chr1_10kb[5813] = "C"  # A -> C
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp  [3000:3200] ---> total length 10200
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 3100 bp  [3100:2800] --> total length 9000
deleted = c(inserted[1:2800], inserted[3101:10200])
print(length(deleted))
## SV, inversion 150 bp, at 3050 bp  [3050:3200]
inverted = c(deleted[1:3050], reverse(deleted[3050:3200]), deleted[3201:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp [3150:3300]
final = c(inverted[1:3150], inverted[3151:3200], inverted[3151:3200], inverted[3151:3200], inverted[3301:9000])
export(final, "overlap3_chr1_hg38.fasta.gz", compress=TRUE)
##
## overlapping, chrom2
##
chrom2 = human_ref_genome$chr2
section_chr1_10kb = chrom2[6000000:6009999]
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## T T A C C A G T T G A T
section_chr1_10kb[2326] = "G"  # T -> G
section_chr1_10kb[2602] = "G"  # T -> G
section_chr1_10kb[3171] = "G"  # A -> G
section_chr1_10kb[3328] = "T"  # C -> T
section_chr1_10kb[3394] = "T"  # C -> T
section_chr1_10kb[3399] = "G"  # A -> G
section_chr1_10kb[4088] = "A"  # G -> A
section_chr1_10kb[4569] = "C"  # T -> C
section_chr1_10kb[5313] = "C"  # T -> C
section_chr1_10kb[5387] = "A"  # G -> A
section_chr1_10kb[5547] = "G"  # A -> G
section_chr1_10kb[5813] = "G"  # T -> G
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp  [3000:3200] ---> total length 10200
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 3100 bp  [3100:2800] --> total length 9000
deleted = c(inserted[1:2800], inserted[3101:10200])
print(length(deleted))
## SV, inversion 150 bp, at 3050 bp  [3050:3200]
inverted = c(deleted[1:3050], reverse(deleted[3050:3200]), deleted[3201:9000])
print(length(inverted))
# SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp [3150:3300]
final = c(inverted[1:3150], inverted[3151:3200], inverted[3151:3200], inverted[3151:3200], inverted[3301:9000])
export(final, "overlap_chr2_hg38.fasta.gz", compress=TRUE)
##
## chrom3, overlapping
##
chrom3 = human_ref_genome$chr3
section_chr1_10kb = chrom3[6000000:6009999]
## T -> G/C
## A -> G/C
## G -> T/A
## C -> T/A
## A A A C A T A C G A C T
section_chr1_10kb[2326] = "G"  # A -> G
section_chr1_10kb[2602] = "G"  # A -> G
section_chr1_10kb[3171] = "G"  # A -> G
section_chr1_10kb[3328] = "G"  # C -> G
section_chr1_10kb[3394] = "G"  # A -> G
section_chr1_10kb[3399] = "G"  # T -> G
section_chr1_10kb[4088] = "G"  # A -> G
section_chr1_10kb[4569] = "T"  # C -> T
section_chr1_10kb[5313] = "T"  # G -> T
section_chr1_10kb[5387] = "G"  # A -> G
section_chr1_10kb[5547] = "T"  # C -> T
section_chr1_10kb[5813] = "G"  # T -> G
library(Biostrings)
insert2 = DNAString("TAGACGAGAGCGTTGCTCATTTCTCATTGCCTACCCGGCCGGTAGTGACGATTTGTACGATACTTATATTGTCAGCCACCGCGGCGTCTAGCCACACGTACGCATTTCGATTGAAATTTCTGCTCGTGCTGTCCTGCCTGTAATTGCACATATTACTGTGCCGTACCGTGGGGACCCTAACTCTCCTGCTTTAAAAAAGA")
## length(insert2) -- 200
## insertion at 3000 bp  [3000:3200] ---> total length 10200
inserted = c(section_chr1_10kb[1:3000], insert2, section_chr1_10kb[3001:10000])
print(length(inserted))
## deletion 300 bp, at 3100 bp  [3100:2800] --> total length 9000
deleted = c(inserted[1:2800], inserted[3101:10200])
print(length(deleted))
## SV, inversion 150 bp, at 3050 bp  [3050:3200]
inverted = c(deleted[1:3050], reverse(deleted[3050:3200]), deleted[3201:9000])
print(length(inverted))
## SV, tandem repeat 150 bp at 5000 bp, 3 times 50bp [3150:3300]
final = c(inverted[1:3150], inverted[3151:3200], inverted[3151:3200], inverted[3151:3200], inverted[3301:9000])
export(final, "overlap_chr3_hg38.fasta.gz", compress=TRUE)































