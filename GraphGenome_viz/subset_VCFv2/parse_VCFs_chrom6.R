

## April 17, 2018

library(VariantAnnotation)
library(data.table)
library(ggplot2)


### http://www.internationalgenome.org/wiki/Analysis/vcf4.0/

## if length(ALT) > length(REF), then insertion
## if length(ALT) < length(REF), then deletion

## dt = gr2dt(foo) ## convert GRanges to data.table
##  INTERESTING BUG!!
## ALT is a <DNAStringSet>


chrom6 = readVcf('chr6_only.vcf') 

fixed_chrom6 = fixed(chrom6)

fixed_chrom6$QUAL = NULL

fixed_chrom6$FILTER = NULL

fixed_chrom6$ALT = CharacterList(fixed_chrom6$ALT)

dt = as.data.table(fixed_chrom6)

print(dim(dt))

## CHR 6
## 1247905       2 


dt$REF = as.character(dt$REF)

dt$ALT = as.character(dt$ALT)


counts = dt[, do.call(rbind, mapply(function(x, snd) {
    lens = nchar(x[x!=""])
    longer = lens[lens > snd]
    if (length(longer) == 0L){
        longer = 0L
    }
    shorter = lens[lens < snd]
    if (length(shorter) == 0L){
        shorter = 0L
    }
    list(list(longer), list(shorter))
}, strsplit(ALT, ",| "), nchar(REF), SIMPLIFY=FALSE)), by=seq_len(dt[,.N])]



counts[, INS_length := V1]
counts[, DEL_length := V2]
counts[, VCF_pos := seq_len]
counts[, V1:= NULL]
counts[, V2:=NULL]
counts[, seq_len:=NULL]

setcolorder(counts, c("VCF_pos", "INS_length", "DEL_length"))


## now, create vector of insertions and deletions


total_insertions = unlist(counts$INS_length)
total_deletions = unlist(counts$DEL_length)


print(length(total_insertions))
## 1431899 

print(length(total_deletions))
## 1360443 



### remove the zeros, as we won't need them

total_insertions = total_insertions[total_insertions != 0]

total_deletions = total_deletions[total_deletions != 0]


print(length(total_insertions))
## 458714 

print(length(total_deletions))
## 631668 



insdt = data.table(insertions = total_insertions)


deldt = data.table(deletions = total_deletions)


##  > head(table(deldt))
##  deldt
##  1      2      3      4      5      6
##  487216  10738   4115   3595  11716   3612 

## > head(table(insdt))
## insdt
## 2     3     4     5     6     7
## 91724  11349  10789 164363  78291  32563  


## first plot lengths from 2 to 50, and then the rest



print(max(insdt$insertions))
## 290657 

print(max(deldt$deletions))
## 177005

library(skitools)


plt1 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(1, 50)) + ylim(0, 165000) + 
    labs(title="Chrom6, Insertion lengths <= 50, binwidth=1") + labs(y="Count", x="Insertion Lengths")


ppdf(print(plt1))



plt2 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(1, 50)) + ylim(0, 165000) + 
    labs(title="Chrom6, Deletions lengths <= 50, binwidth=1") + labs(y="Count", x="Deletion Lengths")  



ppdf(print(plt2))


########### now the rest

print(dim(insdt)[1])
## 458714

print(dim(deldt)[1])  
## 631668 


plt3 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 229236)) +
    ylim(0, 1200) + 
    labs(title="Chrom6, Insertion lengths > 50, binwidth=5") + labs(y="Count", x="Insertion Lengths")  


ppdf(print(plt3))



plt4 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 229236)) +
    ylim(0, 1200) + 
    labs(title="Chrom6, Deletion lengths > 50, binwidth=5") + labs(y="Count", x="Deletion Lengths")

ppdf(print(plt4))


fwrite(insdt, file = "chrom6_insertions.csv", sep = ",")

fwrite(deldt, file = "chrom6_deletions.csv", sep = ",")

print('end')

