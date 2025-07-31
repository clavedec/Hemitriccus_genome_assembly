# Genome assembly of _Hemitriccus margaritaceiventer_

Scripts used to perform the assembly, scaffolding, masking repeats and annotation of _Hemitriccus margaritaceiventer_ (Tyrannidae, code: hemMar).

## Genome assembly

We used a Snakemake workflow from https://github.com/harvardinformatics/pacbio_hifi_assembly and the config file `config.assembly.yaml`.

```
snakemake --snakefile Snakefile --cores 20 --use-conda --rerun-incomplete
``` 

## Scaffolding

To scaffold the genome assembly, we used ArimaGenomics pipeline for Hi-C data processing, YAHS, and Juicebox Assembly tools.

### Purge duplicates
Hmar.HiFi.HMRG_6433.p_ctg.fa = assembled genome
Hmar.HiFi.HMRG_6433_combined_reads.fastq.gz = combined fastq reads after bam2fastx

```
minimap2 -x map-ont Hmar.HiFi.HMRG_6433.p_ctg.fa Hmar.HiFi.HMRG_6433_combined_reads.fastq.gz | gzip -c - > hemMar.paf.gz

./purge_dups/bin/pbcstat *.paf.gz
./purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log
./purge_dups/bin/split_fa Hmar.HiFi.HMRG_6433.p_ctg.fa > Hmar.HiFi.HMRG_6433.p_ctg.split

minimap2 -xasm5 -DP Hmar.HiFi.HMRG_6433.p_ctg.split Hmar.HiFi.HMRG_6433.p_ctg.split | gzip -c - > Hmar.HiFi.HMRG_6433.p_ctg.split.self.paf.gz

./purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov Hmar.HiFi.HMRG_6433.p_ctg.split.self.paf.gz > dups.bed 2> purge_dups.log
./purge_dups/bin/get_seqs -e dups.bed Hmar.HiFi.HMRG_6433.p_ctg.fa > haplotype.fasta
```
### ArimaGenomics 

We followed the pipeline from ArimaGenomics (https://github.com/ArimaGenomics/mapping_pipeline) to process HiC data, using with the script `arima.sh`.

### YaHS

We used YaHS (https://github.com/c-zhou/yahs) to perform the scaffolding of the hemMar genome assembly using HiC.
```
yahs Hmar.HiFi.HMRG_6433.p_ctg.purged.fa HiC.merged.MarkDups.bam
```
The output was then used to generate Hi-C contact maps in Juicebox Assembly Tools (https://github.com/aidenlab/Juicebox) for manual refinement. The scaffolded genome has the name `hemMar_complete_sorted_JBAT.FINAL.fa`. 

### Synteny mapping

We next compared the refined scaffolded genome with Lance-tailed manakin (_Chiroxiphia lanceolata)_ (ChiLan1.pri; GCF_009829145.1) and Zebra Finch (_Taeniopygia guttata)_ (TaeGut1.4.pri; GCF_003957565.2) assemblies, which were two of the closest relatives with chromosome-level assemblies. We identified putative chromosomes based on the synteny mapping results. 

```
~/mummer-4.0.0rc1/nucmer -c 100 -t 8 GCF_009829145.1_bChiLan1.pri_genomic.fna  hemMar_complete_sorted_JBAT.FINAL.fa -p Hmar.chiLan.alignment.fix.nucmer

~/mummer-4.0.0rc1/delta-filter -1 -l 500 Hmar.chiLan.alignment.fix.nucmer.delta > Hmar.chiLan.alignment.fix.nucmer.delta.filt

~/mummer-4.0.0rc1/show-coords -H -c -l -o -r -T Hmar.chiLan.alignment.fix.nucmer.delta.filt > Hmar.chiLan.alignment.fix.nucmer.delta.filt.coords
```
The synteny was later assessed using Circo plots.

## Annotating and masking repeats

We annotated and masked repeatable elements using RepeatModeler and RepeatMasker following the pipeline from https://darencard.net/blog/2022-07-09-genome-repeat-annotation/.

### Building a repeat database

Building a database using RepeatModeler and then labeling the annotated repeats as 'known', or 'unknown' if not annotated.
```
BuildDatabase -name hemMar -engine ncbi hemMar_complete_sorted_JBAT.FINAL.fa
RepeatModeler -threads 16 -engine ncbi -database hemMar 2>&1 | tee 00_repeatmodeler.log

cat hemMar-families.fa | seqkit fx2tab | awk '{ print "hemMar_"$0 }' | seqkit tab2fx > hemMar-families.prefix.fa
cat hemMar-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > hemMar-families.prefix.fa.known
cat hemMar-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > hemMar-families.prefix.fa.unknown
```
Using `repclassifier` from https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier to help build the database with known repeats.
```
repclassifier -t 3 -d Tetrapoda \
-u 01.build.database/hemMar-families.prefix.fa.unknown \
-k 01.build.database/hemMar-families.prefix.fa.known \
-a 01.build.database/hemMar-families.prefix.fa.known \
-o round-1_RepbaseTetrapoda-Self ## rounds from 1-5
```
### Annotating and masking repeats with RepeatMasker

```
# Round1: Simple repeats
~/repeat-annotation/RepeatMasker/RepeatMasker -pa 16 -a -e ncbi -dir 01_simple_out -noint -xsmall hemMar_complete_sorted_JBAT.FINAL.fa 2>&1 | tee logs/01_simplemask.log

# Round 2: Tetrapoda elements
~/repeat-annotation/RepeatMasker/RepeatMasker -pa 16 -a -e ncbi -dir 02_tetrapoda_out -nolow \
-species tetrapoda 01_simple_out/hemMar_complete_sorted_JBAT.FINAL.simple_mask.masked.fasta 2>&1 | tee logs/02_tetrapodamask.log

rename simple_mask.masked.fasta tetrapoda_mask 02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL*
rename .masked .masked.fasta 02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL*

# Round 3: Known elements from RepeatModeler database
~/repeat-annotation/RepeatMasker/RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow \
-lib ../02.repeat.classifier/round-4_Self/round-4_Self.known \
02_tetrapoda_out/hemMar_complete_sorted_JBAT.FINAL.tetrapoda_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log

rename tetrapoda_mask.masked.fasta known_mask 03_known_out/yahs.out_scaffolds_final*
rename .masked .masked.fasta 03_known_out/yahs.out_scaffolds_final*

# Round 4: Unknown elements from RepeatModeler
~/repeat-annotation/RepeatMasker/RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow \
-lib ../02.repeat.classifier/round-4_Self/round-4_Self.unknown \
03_known_out/hemMar_complete_sorted_JBAT.FINAL.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log

rename known_mask.masked.fasta unknown_mask 04_unknown_out/hemMar_complete_sorted_JBAT.FINAL*
rename .masked .masked.fasta 04_unknown_out/hemMar_complete_sorted_JBAT.FINAL*

# Round 5: Repeats found in sparrow (Benham et al. 2024, GBE)
~/repeat-annotation/RepeatMasker/RepeatMasker -pa 16 -a -e ncbi -dir 06_sparrow_TE -nolow \
-lib Passerellidae.final_TElibrary.fa \
04_unknown_out/hemMar_complete_sorted_JBAT.FINAL.unknown_mask.masked.fasta 2>&1 | tee logs/06_sparowmask.log
```
The results were combined using the script `combine.repeatmasker.output.sh`.

## Annotation

We annotated the genome using three different approaches: using (1) transcriptome-based annotation from the RNA sequencing data using StringTie and TransDecoder; (2) BRAKER 3; and (3) TOGA.

### RNA-seq
Trimming RNA fastq.gz paired reads with Trim_Galore.

```
trim_galore --paired -j 1 --retain_unpaired --phred33 --output_dir ./ --length 35 -q 0 --stringency 5 $FQ1 $FQ2
```
Mapping filtered RNA sequences from each individual to the hemMar draft genome using STAR.
```
# Indexing draft genome
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles $GENOME

# Mapping to hemMar draft genome
STAR --twopassMode Basic  --runThreadN 8 --genomeDir $GENOME \
--readFilesCommand "gunzip -c" --readFilesIn $FQ1 $FQ2 \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${FQ1%%merged.trimmed.R1.fq.gz}

# Sorting the BAM files
samtools sort -m 10G -o ${BAM%%.out.bam}.sorted.bam $BAM
samtools index ${BAM%%.out.bam}.sorted.bam

# Collecting statistics of the mapping
java -jar picard.jar CollectAlignmentSummaryMetrics \
--INPUT ${BAM%%.out.bam}.sorted.bam \
--REFERENCE_SEQUENCE $REFERENCE \
--OUTPUT $OUTPUT_METRICS_FILE
```
### Stringtie and TransDecoder
Our first approach was to use the mapped RNA reads to generate a GFF3 file using Stringtie (https://github.com/gpertea/stringtie) and TransDecoder (https://github.com/TransDecoder/TransDecoder).

```
stringtie "$BAM_FILE" -o "${BASENAME}.gtf" -p 4
stringtie --merge -o RNA.merged.gtf  gtf.list
gffcompare -o gffcompare_output RNA.merged.gtf

~/TransDecoder-TransDecoder-v5.7.0/util/gtf_genome_to_cdna_fasta.pl RNA.merged.gtf $GENOME > transcriptome.fasta
~/TransDecoder-TransDecoder-v5.7.0/util/gtf_to_alignment_gff3.pl RNA.merged.gtf > transcriptome.gff3  # GTF to GFF
~/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t transcriptome.fasta  # Generate best candidate ORF predictions
~/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t transcriptome.fasta  # Generate best candidate ORF predictions
~/TransDecoder-TransDecoder-v5.7.0/util/cdna_alignment_orf_to_genome_orf.pl \
    transcriptome.fasta.transdecoder.gff3 \
    transcriptome.gff3 \
    transcriptome.fasta > transcriptome_transdecoder.gff3   
   ```

### BRAKER 3
Our second approach was to use BRAKER3 (https://github.com/Gaius-Augustus/BRAKER) to annotate the hemMar draft genome using the RNA-seq data.

```
brakersif=braker3.sif
myspecies=Hemitriccus_margaritaceiventer_v1
genome=hemMar_complete_sorted_JBAT.FINAL.full.soft.mask.fasta
proteindbase=Vertebrata.fa
bamsdir=bamfiles.txt

bamfiles=$(cat $bamsdir | awk '$1=$1' RS= OFS=,)

singularity exec --cleanenv $brakersif cp -Rs /opt/Augustus/config/ augustus_config

singularity exec --no-home \
                 --home /opt/gm_key \
                 --cleanenv \
                 --env AUGUSTUS_CONFIG_PATH=${PWD}/augustus_config \
                 $brakersif braker.pl \
                 --prot_seq=${proteindbase} \
                 --bam=${bamfiles} \
                 --species=${myspecies} \
                 --genome=${genome} \
                 --threads=48 \
                 --softmasking \
                 --etpmode \
                 --gff3
```
Functional annotation was done using Interproscan.

```
./interproscan.sh -appl SUPERFAMILY,PANTHER,Pfam -i proteins.clean.fa -o interproscan_output.gff -f gff3 -goterms -iprlookup -pa -cpu 48 --verbose
```                 
### TOGA

Our last approach was to perform annotation based on ortholog genes of different species using TOGA (https://github.com/hillerlab/TOGA). This method uses pairwise genome alignment chains between an annotated reference genome (here chicken, _Gallus gallus,_ GRCg6a_;_ Zebra Finch, _Taenopygia guttata_, bTaeGut1.4, and _Chiroxiphia lanceolata_, bChiLan1) and the query species (here _H. margaritaceiventer_).

First, we performed an aligment using Cactus. Below, the example using _C. lanceolata_.
```
singularity exec --cleanenv cactus_v2.9.3.sif cactus-preprocess temp/0 ./hemMar.chiLan.txt cactus_prepare/hemMar.chiLan.txt --inputNames hemMar chiLan   --logFile cactus_prepare/logs/preprocess-hemMar.log
singularity exec --cleanenv cactus_v2.9.3.sif cactus-blast temp/1 cactus_prepare/hemMar.chiLan.txt cactus_prepare/Anc0.paf --root Anc0   --logFile cactus_prepare/logs/blast-Anc0.log
singularity exec --cleanenv cactus_v2.9.3.sif cactus-align temp/2 cactus_prepare/hemMar.chiLan.txt cactus_prepare/Anc0.paf cactus_prepare/hemMar.chiLan.hal --root Anc0   --logFile cactus_prepare/logs/align-Anc0.log
singularity exec --cleanenv cactus_v2.9.3.sif hal2fasta cactus_prepare/hemMar.chiLan.hal Anc0 --hdf5InMemory > cactus_prepare/Anc0.fa
```
We then performed an exon-aware liftover to prepare the files for TOGA following the steps described here: https://github.com/harvardinformatics/AnnotationTOGA.

Finally, we performed TOGA.
```
./toga.py hemMar.chiLan.chain \
CDS_chiLan.corrected.bed chiLan.2bit \
hemMar.2bit --kt --pn hemMar.chiLan.toga \
-i chiLan_CDS_isoforms.tsv \
--nc nextflow_config_files \
--cb 10,100,200 --cjn 750
```


