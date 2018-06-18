# Graph Genome of Human Assemblies

A genome graph representation of seven ethnically-diverse human genomes, constructed using an algorithmically novel approach.

This project was initiated at an NCBI-style hackathon held before the 2016 Biological Data Science meeting at Cold Spring Harbor Laboratory in October, 2016.


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

## Contributors

Evan Biederstedt (NYGC, WCM)

Alexander Dilthey (NHGRI-NIH, HHU/UKD)

Nathan Dunn (LBNL)

Aarti Jajoo (Baylor)

Nancy Hansen (NIH)

Jeff Oliver (Arizona)

Andrew Olsen (CSHL)
