

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


chrom8 = readVcf('chr8_only.vcf') 

fixed_chrom8 = fixed(chrom8)

fixed_chrom8$QUAL = NULL

fixed_chrom8$FILTER = NULL

fixed_chrom8$ALT = CharacterList(fixed_chrom8$ALT)

dt = as.data.table(fixed_chrom8)

print(dim(dt))

## CHR 8
##  141904      2  


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
## 164831 

print(length(total_deletions))
## 145887 



### remove the zeros, as we won't need them

total_insertions = total_insertions[total_insertions != 0]

total_deletions = total_deletions[total_deletions != 0]


print(length(total_insertions))
## 53282 

print(length(total_deletions))
##  53213



insdt = data.table(insertions = total_insertions)


deldt = data.table(deletions = total_deletions)


##  > head(table(deldt))
##  deldt
##  1      2      3      4      5      6
## 45644  1128   455   388  1311   411  

## > head(table(insdt))
## insdt
## 2     3     4     5     6     7
## 7888   932  1267 20899 10041  4049 


## first plot lengths from 2 to 50, and then the rest



print(max(insdt$insertions))
## 75470        

print(max(deldt$deletions))
##  64626 



library(skitools)


plt1 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(2, 50)) + ylim(0, 21000) +  
    labs(title="Chrom8, Insertion lengths <= 50, binwidth=1") + labs(y="Count", x="Insertion Lengths")


ppdf(print(plt1))



plt2 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=1, col="#46ACC8") + scale_x_continuous(limits = c(2, 50)) + ylim(0, 21000) + 
    labs(title="Chrom8, Deletions lengths <= 50, binwidth=1") + labs(y="Count", x="Deletion Lengths")  



ppdf(print(plt2))


########### now the rest

print(dim(insdt)[1])
## 53282  

print(dim(deldt)[1])  
## 53213 


plt3 = ggplot(insdt, aes(x=insertions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 229236)) +
    ylim(0, 110) + 
    labs(title="Chrom8, Insertion lengths > 50, binwidth=5") + labs(y="Count", x="Insertion Lengths")  


ppdf(print(plt3))


plt4 = ggplot(deldt, aes(x=deletions)) + geom_histogram(binwidth=5, col="#46ACC8") + scale_x_continuous(limits = c(51, 229236)) +
    ylim(0, 110) + 
    labs(title="Chrom8, Deletion lengths > 50, binwidth=5") + labs(y="Count", x="Deletion Lengths")

ppdf(print(plt4))


fwrite(insdt, file = "chrom8_insertions.csv", sep = ",")

fwrite(deldt, file = "chrom8_deletions.csv", sep = ",")



print('end')

