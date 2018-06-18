

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


chrom14 = readVcf('chr14_only.vcf') 

fixed_chrom14 = fixed(chrom14)

fixed_chrom14$QUAL = NULL

fixed_chrom14$FILTER = NULL

fixed_chrom14$ALT = CharacterList(fixed_chrom14$ALT)

dt = as.data.table(fixed_chrom14)

print(dim(dt))

## CHR 14
## 654187      2 


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
## 740311 

print(length(total_deletions))
## 665281 



### remove the zeros, as we won't need them

total_insertions = total_insertions[total_insertions != 0]

total_deletions = total_deletions[total_deletions != 0]


print(length(total_insertions))
## 227076 

print(length(total_deletions))
## 298615



insdt = data.table(insertions = total_insertions)


deldt = data.table(deletions = total_deletions)


##  > head(table(deldt))
##  deldt
##  1      2      3      4      5      6
##  271587   5357   2125   1755   5824   1773  

## > head(table(insdt))
## insdt
## 2     3     4     5     6     7
## 50700  6189  5343 79015 38900 15917 


## first plot lengths from 2 to 50, and then the rest



print(max(insdt$insertions))
## 62030  

print(max(deldt$deletions))
## 201048 


library(skitools)


plt1 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(1, 50)) + ylim(0, 79500) + 
    labs(title="Chrom14, Insertion lengths <= 50, binwidth=1") + labs(y="Count", x="Insertion Lengths")


ppdf(print(plt1))



plt2 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(1, 50)) + ylim(0, 79500) + 
    labs(title="Chrom14, Deletions lengths <= 50, binwidth=1") + labs(y="Count", x="Deletion Lengths")  



ppdf(print(plt2))


########### now the rest

print(dim(insdt)[1])
## 227076

print(dim(deldt)[1])  
## 298615 


plt3 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 298615)) +
    ylim(0, 170) + 
    labs(title="Chrom14, Insertion lengths > 50, binwidth=5") + labs(y="Count", x="Insertion Lengths")  


ppdf(print(plt3))


plt4 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 298615)) +
    ylim(0, 170) + 
    labs(title="Chrom14, Deletion lengths > 50, binwidth=5") + labs(y="Count", x="Deletion Lengths")

ppdf(print(plt4))



fwrite(insdt, file = "chrom14_insertions.csv", sep = ",")


fwrite(deldt, file = "chrom14_deletions.csv", sep = ",") 


print('end')


