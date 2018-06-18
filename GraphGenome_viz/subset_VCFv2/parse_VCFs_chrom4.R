

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


chrom4 = readVcf('chr4_only.vcf') 

fixed_chrom4 = fixed(chrom4)

fixed_chrom4$QUAL = NULL

fixed_chrom4$FILTER = NULL

fixed_chrom4$ALT = CharacterList(fixed_chrom4$ALT)

dt = as.data.table(fixed_chrom4)

print(dim(dt))

## CHR 4
## 1453909       2 


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
## 1697927 

print(length(total_deletions))
## 1596935



### remove the zeros, as we won't need them

total_insertions = total_insertions[total_insertions != 0]

total_deletions = total_deletions[total_deletions != 0]


print(length(total_insertions))
## 559496

print(length(total_deletions))
## 715576 



insdt = data.table(insertions = total_insertions)


deldt = data.table(deletions = total_deletions)


##  > head(table(deldt))
##  deldt
##  1      2      3      4      5      6
##   535071  12855   4969   4141  13238   4097 

## > head(table(insdt))
## insdt
## 2     3     4     5     6     7
## 101781  12576  12628 195267  93146  38550 


## first plot lengths from 2 to 50, and then the rest



print(max(insdt$insertions))
## 147694

print(max(deldt$deletions))
## 147106

library(skitools)


plt1 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(1, 50)) + ylim(0, 200000) + 
    labs(title="Chrom4, Insertion lengths <= 50, binwidth=1") + labs(y="Count", x="Insertion Lengths")


ppdf(print(plt1))



plt2 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(1, 50)) + ylim(0, 200000) + 
    labs(title="Chrom4, Deletions lengths <= 50, binwidth=1") + labs(y="Count", x="Deletion Lengths")  



ppdf(print(plt2))


########### now the rest

print(dim(insdt)[1])
## 559496

print(dim(deldt)[1])  
## 715576  

 
plt3 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 715576)) + ylim(0, 430) + 
    labs(title="Chrom4, Insertion lengths > 50, binwidth=5") + labs(y="Count", x="Insertion Lengths")  


ppdf(print(plt3))


plt4 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 715576)) + ylim(0, 430) + 
    labs(title="Chrom4, Deletion lengths > 50, binwidth=5") + labs(y="Count", x="Deletion Lengths")

ppdf(print(plt4))



fwrite(insdt, file = "chrom4_insertions.csv", sep = ",")

fwrite(deldt, file = "chrom4_deletions.csv", sep = ",")



print('end')

