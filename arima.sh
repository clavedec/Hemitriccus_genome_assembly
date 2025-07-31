#! /bin/bash
SRA='hemMar_HiC'
LABEL='hemMar_assembly'
BWA='/n/home01/cfcarvalho/bwa/bwa'
SAMTOOLS='/n/home01/cfcarvalho/samtools-1.18/samtools'
IN_DIR='./'
REF='./Hmar.HiFi.HMRG_6433.p_ctg.purged.fa'
FAIDX='$REF.fai'
RAW_DIR='raw_dir'
FILT_DIR='filt_dir'
FILTER='./arima/filter_five_end.pl'
COMBINER='./arima/two_read_bam_combiner.pl'
STATS='./arima/get_stats.pl'
PICARD='/n/home01/cfcarvalho/picard.jar'
TMP_DIR='tmp_dir'
PAIR_DIR='pair_dir'
REP_DIR='rep_dir'
REP_LABEL=${LABEL}_rep1
MAPQ_FILTER=10
CPU=12

echo "### Step 0: Check output directories’ existence & create them as
needed"
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

## Getting index
$SAMTOOLS faidx $REF
$BWA index $REF

## Lane 003
echo "### Step 1.A1: FASTQ to BAM (1st)"
$BWA mem -t $CPU $REF T_2774_HiC_S1_L003_R1_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb > HiC1_R1.bam
echo "### Step 1.B1: FASTQ to BAM (2nd)"
$BWA mem -t $CPU $REF T_2774_HiC_S1_L003_R2_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb > HiC1_R2.bam
## Lane 004
echo "### Step 1.A2: FASTQ to BAM (3nd)"
$BWA mem -t $CPU $REF T_2774_HiC_S1_L004_R1_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb > HiC2_R1.bam
echo "### Step 1.B2: FASTQ to BAM (4nd)"
$BWA mem -t $CPU $REF T_2774_HiC_S1_L004_R2_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb > HiC2_R2.bam

echo "### Step 2.A1: Filter 5' end (1st)"
$SAMTOOLS view -h HiC1_R1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/HiC1_filt_R1.bam
echo "### Step 2.B1: Filter 5' end (2nd)"
$SAMTOOLS view -h HiC1_R2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/HiC1_filt_R2.bam
echo "### Step 2.A2: Filter 5' end (3rd)"
$SAMTOOLS view -h HiC2_R1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/HiC2_filt_R1.bam
echo "### Step 2.B2: Filter 5' end (4th)"
$SAMTOOLS view -h HiC2_R2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/HiC2_filt_R2.bam

echo "### Step 3.A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/HiC1_filt_R1.bam $FILT_DIR/HiC1_filt_R2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/HiC1_merged.bam
perl $COMBINER $FILT_DIR/HiC2_filt_R1.bam $FILT_DIR/HiC2_filt_R2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/HiC2_merged.bam

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/HiC1_merged.bam OUTPUT=$PAIR_DIR/HiC1_merged.ReadGroups.bam ID=HiC1.merged.ReadGroups LB=HiC1.merged.ReadGroups SM=$LABEL PL=ILLUMINA PU=none
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/HiC2_merged.bam OUTPUT=$PAIR_DIR/HiC2_merged.ReadGroups.bam ID=HiC2.merged.ReadGroups LB=HiC2.merged.ReadGroups SM=$LABEL PL=ILLUMINA PU=none

echo "### STEP 3.C Removing technical replicates"
java -jar ~/picard.jar MergeSamFiles INPUT=$PAIR_DIR/HiC1_merged.ReadGroups.bam INPUT=$PAIR_DIR/HiC2_merged.ReadGroups.bam OUTPUT=$TMP_DIR/HiC.merged.runs.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -jar ~/picard.jar MarkDuplicates INPUT=$TMP_DIR/HiC.merged.runs.bam OUTPUT=$REP_DIR/HiC.merged.MarkDups.bam \
METRICS_FILE=$REP_DIR/metrics.HiC.merged.MarkDups.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
$SAMTOOLS index $REP_DIR/HiC.merged.MarkDups.bam
perl $STATS $REP_DIR/HiC.merged.MarkDups.bam > $REP_DIR/HiC.merged.MarkDups.bam.stats

java -jar ~/picard.jar CollectRawWgsMetrics \
        I=$REP_DIR/HiC.merged.MarkDups.bam \
        O=$REP_DIR/HiC.merged.MarkDups.dups.bam.raw_wgs_metrics.txt \
        R=$REF \
        INCLUDE_BQ_HISTOGRAM=true

echo "Finished Mapping Pipeline through Duplicate Removal"

