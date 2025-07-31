#!/bin/bash

Combining results from Repeatmasker on hemMar

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/hemMar_complete_sorted_JBAT.FINAL.simple_mask.cat.gz \
02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL.tetrapoda_mask.cat.gz \
03_known_out/hemMar_complete_sorted_JBAT.FINAL.known_mask.cat.gz \
04_unknown_out/hemMar_complete_sorted_JBAT.FINAL.unknown_mask.cat.gz \
06_sparrow_TE/hemMar_complete_sorted_JBAT.FINAL.sparrow_mask.cat.gz \
> 05_full_out/hemMar_complete_sorted_JBAT.FINAL.full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/hemMar_complete_sorted_JBAT.FINAL.simple_mask.out \
<(cat 02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL.tetrapoda_mask.out | tail -n +4) \
<(cat 03_known_out/hemMar_complete_sorted_JBAT.FINAL.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/hemMar_complete_sorted_JBAT.FINAL.unknown_mask.out | tail -n +4) \
<(cat 06_sparrow_TE/hemMar_complete_sorted_JBAT.FINAL.sparrow_mask.out | tail -n +4) \
> 05_full_out/hemMar_complete_sorted_JBAT.FINAL.full_mask.out

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/hemMar_complete_sorted_JBAT.FINAL.simple_mask.out > 05_full_out/hemMar_complete_sorted_JBAT.FINAL.simple_mask.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL.tetrapoda_mask.out \
<(cat 03_known_out/hemMar_complete_sorted_JBAT.FINAL.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/hemMar_complete_sorted_JBAT.FINAL.unknown_mask.out | tail -n +4) \
<(cat 06_sparrow_TE/hemMar_complete_sorted_JBAT.FINAL.sparrow_mask.out | tail -n +4) \
> 05_full_out/hemMar_complete_sorted_JBAT.FINAL.complex_mask.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_out/hemMar_complete_sorted_JBAT.FINAL.simple_mask.align \
02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL.tetrapoda_mask.align \
03_known_out/hemMar_complete_sorted_JBAT.FINAL.known_mask.align \
04_unknown_out/hemMar_complete_sorted_JBAT.FINAL.unknown_mask.align \
06_sparrow_TE/hemMar_complete_sorted_JBAT.FINAL.sparrow_mask.align \
> 05_full_out/hemMar_complete_sorted_JBAT.FINAL.full_mask.align

# calculate the length of the genome sequence in the FASTA
allLen=`seqtk comp hemMar_complete_sorted_JBAT.FINAL.fa | datamash sum 2`;

# calculate the length of the N sequence in the FASTA
nLen=`seqtk comp hemMar_complete_sorted_JBAT.FINAL.fa | datamash sum 9`;

# tabulate repeats per subfamily with total bp and proportion of genome masked
cat 05_full_out/hemMar_complete_sorted_JBAT.FINAL.full_mask.out | tail -n +4 | awk -v OFS="\t" '{ print $6, $7, $11 }' |
awk -F '[\t/]' -v OFS="\t" '{ if (NF == 3) print $3, "NA", $2 - $1 +1; else print $3, $4, $2 - $1 +1 }' |
datamash -sg 1,2 sum 3 | grep -v "\?" |
awk -v OFS="\t" -v genomeLen="${allLen}" '{ print $0, $3 / genomeLen }' > 05_full_out/hemMar_complete_sorted_JBAT.FINAL.full_mask.tabulate
