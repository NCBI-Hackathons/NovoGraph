

## April 17, 2018

library(VariantAnnotation)
library(data.table)
library(ggplot2)
library(skitools)

### http://www.internationalgenome.org/wiki/Analysis/vcf4.0/

## if length(ALT) > length(REF), then insertion
## if length(ALT) < length(REF), then deletion

## dt = gr2dt(foo) ## convert GRanges to data.table
##  INTERESTING BUG!!
## ALT is a <DNAStringSet>


chrom21 = readVcf('chr21_only.vcf') 

fixed_chrom21 = fixed(chrom21)

fixed_chrom21$QUAL = NULL

fixed_chrom21$FILTER = NULL

fixed_chrom21$ALT = CharacterList(fixed_chrom21$ALT)

dt = as.data.table(fixed_chrom21)

print(dim(dt))

## CHR 21
## 353591      2


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
## 431991

print(length(total_deletions))
## 443587



### remove the zeros, as we won't need them

total_insertions = total_insertions[total_insertions != 0]

total_deletions = total_deletions[total_deletions != 0]


print(length(total_insertions))
## 159752

print(length(total_deletions))
## 229236 



insdt = data.table(insertions = total_insertions)


deldt = data.table(deletions = total_deletions)


##  > head(table(deldt))
##  deldt
##  1      2      3      4      5      6
##  128340   3497   1354   1001   3444   1225 

## > head(table(insdt))
## insdt
## 2     3     4     5     6     7
## 24411  2919  3662 45442 25895 13152   


## first plot lengths from 2 to 50, and then the rest



print(max(insdt$insertions))
## 50006

print(max(deldt$deletions))
## 49442 




plt1 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(2, 50)) + ylim(0, 50010) + 
    labs(title="Chrom21, Insertion lengths <= 50, binwidth=1") + labs(y="Count", x="Insertion Lengths")


ppdf(print(plt1))



plt2 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(2, 50)) + ylim(0, 50010) + 
    labs(title="Chrom21, Deletions lengths <= 50, binwidth=1") + labs(y="Count", x="Deletion Lengths")  



ppdf(print(plt2))


########### now the rest

print(dim(insdt)[1])
## 159752  

print(dim(deldt)[1])  
## 229236


plt3 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 229236)) +
    labs(title="Chrom21, Insertion lengths > 50, binwidth=5") + labs(y="Count", x="Insertion Lengths")  


ppdf(print(plt3))


plt4 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 229236)) +
    labs(title="Chrom21, Deletion lengths > 50, binwidth=5") + labs(y="Count", x="Deletion Lengths")

ppdf(print(plt4))


print('end')

