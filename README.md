# Graph Genome of Multiple De Novo Assemblies

An algorithmically novel approach to construct a genome graph representation of multiple de novo sequence assemblies. We then provide a proof of principle by constructing a genome graph of seven ethnically-diverse human genomes. 

This project was initiated at an NCBI-style hackathon held before the 2016 Biological Data Science meeting at Cold Spring Harbor Laboratory in October, 2016.

## Genome Graph of Seven Human Assemblies

Relying upon a linear, one-dimensional character string as the monoploid reference is severely restrictive for genomic research, as such a reference cannot encompass the full breadth of genetic variation which exists. Within human genetics, relying upon a linear monoploid reference genome consequently both limits and biases our understanding into the full diversity of subpopulation variation. Motivated by the potential of genome graphs to address these shortcomings, we present a pipeline for constructing a graph genome from multiple de novo assemblies. We then focused directly on building a graph genome composed of seven human assemblies. We included the following assemblies within the genome graph:
* AK1, Korean
* HX1, Han Chinese
* CHM1, European
* CHM13, European
* HG003, Ashkenazim
* HG004, Ashkenazim
* NA19240, Yoruba


## Pipeline

* reference file GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa (GRCh38 without ALTs)
* contigs file input_sequences.fa

```
## Requires samtools version >= 1.7

bwa index GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa
bwa mem GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa input_sequences.fa  | samtools view -Sb - > input_sequences_unsorted.bam
samtools sort -o SevenGenomesPlusGRCh38Alts.bam input_sequences_unsorted.bam
samtools index SevenGenomesPlusGRCh38Alts.bam

## Check that there are no unmapped reads in the input BAM because this might lead to unknown behaviour:
samtools view -c -f 0x4 SevenGenomesPlusGRCh38Alts.bam

perl BAM2ALIGNMENT.pl --BAM SevenGenomesPlusGRCh38Alts.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta AllContigs.fa --outputFile /intermediate_files/AlignmentInput

## Expect output '../intermediate_files/AlignmentInput.sortedWithHeader'

perl checkBAM_SVs_and_INDELs.pl --BAM /data/projects/phillippy/projects/hackathon/shared/alignments/SevenGenomesPlusGRCh38Alts.bam --referenceFasta /data/projects/phillippy/projects/hackathon/shared/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta /data/projects/phillippy/projects/hackathon/shared/contigs/AllContigs.fa

perl FIND_GLOBAL_ALIGNMENTS.pl --alignmentsFile ../intermediate_files/AlignmentInput.sortedWithHeader --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --outputFile forMAFFT.bam --outputTruncatedReads ../intermediate_files/truncatedReads --outputReadLengths ../intermediate_files/postGlobalAlignment_readLengths

perl BAM2MAFFT.pl --BAM ../intermediate_files/forMAFFT.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --readsFasta AllContigs.fa --outputDirectory .../intermediate_files/forMAFFT --inputTruncatedReads .../intermediate_files/truncatedReads 

perl countExpectedGlobalAlignments.pl --BAM .../intermediate_files/forMAFFT.bam

## assumes you are using SGE to submit jobs
perl CALLMAFFT.pl --action kickOff --mafftDirectory .../intermediate_files/forMAFFT --qsub 1
perl CALLMAFFT.pl --action check --mafftDirectory .../intermediate_files/forMAFFT
perl CALLMAFFT.pl --action reprocess --mafftDirectory .../intermediate_files/forMAFFT

perl globalize_windowbams.pl --fastadir /intermediate_files/forMAFFT/ --msadir /intermediate_files/forMAFFT/ --contigs /intermediate_files/postGlobalAlignment_readLengths --output /intermediate_files/combined_2.sam

perl validate_BAM_MSA.pl --BAM /intermediate_files/combined_sorted.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa

perl BAM2VCF.pl --BAM /intermediate_files/combined_sorted.bam --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output VCF/uber_vcf.vcf

# get CRAM, via SAM sort
samtools view -h -t GRCh38.headerfile.txt /intermediate_files/combined_2.sam > /intermediate_files/combined_2_with_header.sam;\
samtools sort /intermediate_files/combined_2_with_header.sam -o /intermediate_files/combined_2_with_header_sorted.sam

cat /intermediate_files/combined_2_with_header_sorted.sam | samtools view -C -T GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa - > /intermediate_files/combined_2.cram

samtools index /intermediate_files/combined_2.cram

perl checkMAFFT_input_and_output.pl --MAFFTdir /intermediate_files/forMAFFT_2/ --contigLengths /intermediate_files/postGlobalAlignment_readLengths_2 --preMAFFTBAM /intermediate_files/forMAFFT_2.bam --finalOutputCRAM /intermediate_files/combined_2.cram &> output_checkMAFFT_v2

perl CRAM2VCF.pl --CRAM /intermediate_files/combined_2.cram --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output VCF/graph_v2.vcf --contigLengths /intermediate_files/postGlobalAlignment_readLengths_2

perl CRAM2VCF_checkVariantDistribution.pl --output VCF/graph_v2.vcf

perl launch_CRAM2VCF_C++.pl --output VCF/graph_v2.vcf

perl CRAM2VCF_createFinalVCF.pl --CRAM /intermediate_files/combined_2.cram --referenceFasta GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa --output VCF/graph_v2.vcf

```

## Instructions to Download and Process Input Assemblies

The following commands were used to download the assembly FASTAs used for this project:

```
## AK1, Korean:
for file in `echo LPVO02.1.fsa_nt.gz LPVO02.2.fsa_nt.gz LPVO02.3.fsa_nt.gz LPVO02.4.fsa_nt.gz LPVO02.5.fsa_nt.gz LPVO02.6.fsa_nt.gz`; do wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LP/VO/LPVO02/$file; done

## CHM1, European:
for file in `echo LJII01.1.fsa_nt.gz LJII01.10.fsa_nt.gz LJII01.11.fsa_nt.gz LJII01.12.fsa_nt.gz LJII01.13.fsa_nt.gz LJII01.14.fsa_nt.gz LJII01.15.fsa_nt.gz LJII01.2.fsa_nt.gz LJII01.3.fsa_nt.gz LJII01.4.fsa_nt.gz LJII01.5.fsa_nt.gz LJII01.6.fsa_nt.gz LJII01.7.fsa_nt.gz LJII01.8.fsa_nt.gz LJII01.9.fsa_nt.gz`; do echo $file; wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LJ/II/LJII01/$file; done
 
## CHM13, European:
for file in `echo LDOC03.1.fsa_nt.gz LDOC03.2.fsa_nt.gz LDOC03.3.fsa_nt.gz LDOC03.4.fsa_nt.gz LDOC03.5.fsa_nt.gz LDOC03.6.fsa_nt.gz LDOC03.7.fsa_nt.gz`; do echo $file; wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LD/OC/LDOC03/$file; done

## HX1, Han Chinese:
wget http://hx1.wglab.org/data/hx1f4.3rdfixedv2.fa.gz

## HG003, Ashkenazim father:
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/MtSinai_PacBio_Assembly_falcon_03282016/hg003_p_and_a_ctg.fa
 
## HG004, Ashkenazim mother:
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/MtSinai_PacBio_Assembly_falcon_03282016/hg004_p_and_a_ctg.fa

## NA19240, Yoruba:
for file in `echo LKPB01.1.fsa_nt.gz LKPB01.2.fsa_nt.gz LKPB01.3.fsa_nt.gz LKPB01.4.fsa_nt.gz LKPB01.5.fsa_nt.gz LKPB01.6.fsa_nt.gz`; do wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LK/PB/LKPB01/$file; done
```
Upon the successful download of these FASTAs, users should produce a sorted BAM of each assembly by both running BWA-MEM to align against the human reference and sorting the results via SAMtools.



## Contributors

Evan Biederstedt (NYGC, WCM)

Alexander Dilthey (NHGRI-NIH, HHU/UKD)

Nathan Dunn (LBNL)

Aarti Jajoo (Baylor)

Nancy Hansen (NIH)

Jeff Oliver (Arizona)

Andrew Olsen (CSHL)
