#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem=15Gb
#SBATCH -t 0-8:00:00
#SBATCH --output=final_project.%j.out

# =============================================================================
# FINAL PROJECT: Vhl loss in the developing mouse retina (bulk RNA-seq)
# Author:  Sreya Mitra (smitra7), BIOL4324/6324 Bioinformatics, University of Houston
#
# Question: when the Von Hippel-Lindau gene (Vhl) is lost in all retinal cells,
#           which retinal cell types are the main source of the gene-expression changes?
#
# Data:     GEO GSE109102 (Arreola et al. 2018, JCI Insight), whole mouse retina at P7,
#           Illumina HiSeq 2 x 150 bp paired-end reads.
#           4 samples, listed one per line in accessions.txt (SRA run IDs):
#             SRR6464389  WT_rep1       control
#             SRR6464390  WT_rep2       control
#             SRR6464393  Vhl_null_rep1 Vhl-null
#             SRR6464394  Vhl_null_rep2 Vhl-null
#           (sample names and groups are in metadata.tsv)
#
# What this script does (this is the part run on the Sabine cluster):
#   Step 1  download the first 2 million read pairs of each sample
#   Step 2  trim adapters and low-quality bases with fastp
#   Step 3  download the mouse genome (GRCm39) and gene annotation from NCBI
#   Step 4  build a reduced STAR index and align each sample to the genome
#   Step 5  count reads per gene with htseq-count
#   Step 6  merge the four count files into one table (gene.table.tsv)
#
# What happens afterwards (run in R on a laptop, see final_project_analysis.R):
#   DESeq2 differential expression, PCA, cell-type marker tests, pathway tests (fgsea),
#   and sensitivity checks (removing muscle/fat genes and rod genes).
#
# How to run:   sbatch final_project.sh
# Needs:        accessions.txt and metadata.tsv in the working folder; sratoolkit module;
#               fastp, STAR 2.7.11b, htseq-count and the NCBI "datasets" tool installed.
# Output:       gene.table.tsv (genes x samples), aligned/ (BAMs, counts, STAR logs),
#               trimmed/ and fastp_reports/ (trimmed reads and QC reports)
# Limits:       jobs on this cluster are limited to 5 CPUs and 15 Gb of memory.
# =============================================================================

# Working folder for the whole project (everything below is relative to it)
WORK=/project/stuckert/smitra7/final_project
cd $WORK

# -----------------------------------------------------------------------------
# Step 1: download the reads.
# We take only the first 2 million read pairs of each sample (a subset), because
# the shared disk space is limited. "-X" limits the number of reads, "--split-files"
# writes the two mates of each pair to separate files (_1 and _2), "--gzip" compresses them.
# If a sample is already downloaded, it is skipped, so the script can be rerun safely.
# -----------------------------------------------------------------------------
module load sratoolkit/3.4.1
NREADS=2000000
mkdir -p raw_data
while read acc
do
  if [ ! -f raw_data/${acc}_1.fastq.gz ]; then
    echo "Downloading the first $NREADS read pairs of $acc"
    fastq-dump -X $NREADS --split-files --gzip -O raw_data $acc
  else
    echo "$acc already downloaded"
  fi
done < accessions.txt

# -----------------------------------------------------------------------------
# Step 2: quality control and trimming with fastp (paired-end mode).
# fastp finds and removes adapter sequences ("--detect_adapter_for_pe") and low-quality
# bases, and writes an html and a json quality report for each sample into fastp_reports/.
# The raw reads are made read-only so they are never overwritten; trimmed copies go to trimmed/.
# -----------------------------------------------------------------------------
FASTP=$WORK/software/fastp
chmod a-w raw_data/*.fastq.gz
mkdir -p trimmed fastp_reports
while read acc
do
  if [ ! -f trimmed/${acc}_1.trim.fastq.gz ]; then
    echo "Trimming $acc"
    $FASTP -i raw_data/${acc}_1.fastq.gz -I raw_data/${acc}_2.fastq.gz \
      -o trimmed/${acc}_1.trim.fastq.gz -O trimmed/${acc}_2.trim.fastq.gz \
      --detect_adapter_for_pe -w 5 \
      -h fastp_reports/${acc}.html -j fastp_reports/${acc}.json
  else
    echo "$acc already trimmed"
  fi
done < accessions.txt

# -----------------------------------------------------------------------------
# Step 3: download the mouse reference genome and annotation from NCBI.
# Genome: GRCm39 (RefSeq assembly GCF_000001635.27). The GTF file lists where each gene and exon is.
# The "datasets" command-line tool downloads a zip, which is unzipped into reference/.
# -----------------------------------------------------------------------------
export PATH=$PATH:/project/stuckert/bioinformatics/software
mkdir -p reference
if [ ! -d reference/ncbi_dataset ]; then
  datasets download genome accession GCF_000001635.27 --include genome,gtf --filename reference/mouse.zip
  unzip -o reference/mouse.zip -d reference
  rm reference/mouse.zip
fi
GENOME=$(ls reference/ncbi_dataset/data/GCF_000001635.27/*.fna)
GTF=reference/ncbi_dataset/data/GCF_000001635.27/genomic.gtf
ls -lh $GENOME $GTF

# -----------------------------------------------------------------------------
# Step 4: align the reads to the genome with STAR 2.7.11b.
# A full mouse index needs more memory than the 15 Gb job limit, so the index is built in a
# reduced form: "--genomeSAsparseD 8" keeps only every 8th position of the suffix array, and
# "--genomeSAindexNbases 14" is the setting recommended for a genome this size. The temporary
# files are written to the compute node's local disk ($TMPDIR) instead of the shared disk.
# The index is built without the annotation; the GTF is given at the alignment step
# ("--sjdbGTFfile", "--sjdbOverhang 149" = read length 150 minus 1) so STAR knows the splice junctions.
# "--quantMode GeneCounts" makes STAR also count reads per gene, which was used to check that the
# library is unstranded. Output is an unsorted BAM file per sample in aligned/.
# The index is only built if at least one sample still needs aligning.
# -----------------------------------------------------------------------------
STAR=$HOME/STAR-2.7.11b/bin/Linux_x86_64/STAR
mkdir -p aligned
NEEDINDEX=no
while read acc
do
  if [ ! -s aligned/${acc}.Aligned.out.bam ]; then NEEDINDEX=yes; fi
done < accessions.txt

if [ $NEEDINDEX = yes ]; then
  IDX=$TMPDIR/star_index
  mkdir -p $IDX
  $STAR --runMode genomeGenerate --runThreadN 4 --genomeDir $IDX \
    --genomeFastaFiles $GENOME --genomeSAsparseD 8 --genomeSAindexNbases 14 \
    --limitGenomeGenerateRAM 14000000000 \
    --outTmpDir $TMPDIR/tmp_gen --outFileNamePrefix $TMPDIR/idx_
  cp $TMPDIR/idx_Log.out aligned/star_index_Log.out
  while read acc
  do
    if [ ! -s aligned/${acc}.Aligned.out.bam ]; then
      echo "Aligning $acc"
      rm -rf $TMPDIR/tmp_aln
      $STAR --runMode alignReads --genomeDir $IDX \
        --readFilesIn trimmed/${acc}_1.trim.fastq.gz trimmed/${acc}_2.trim.fastq.gz \
        --readFilesCommand zcat --outSAMtype BAM Unsorted \
        --sjdbGTFfile $GTF --sjdbOverhang 149 --quantMode GeneCounts \
        --outTmpDir $TMPDIR/tmp_aln --outFileNamePrefix aligned/${acc}. --runThreadN 5
    fi
  done < accessions.txt
fi

# -----------------------------------------------------------------------------
# Step 5: count how many reads fall in each gene with htseq-count.
# "-f bam -r name" reads the BAM file ordered by read name (paired reads are counted once),
# "-s no" means the library is unstranded, "-t exon" counts reads overlapping exons, and
# "-i gene_id" adds the exon counts up per gene. One count file per sample is written to aligned/.
# -----------------------------------------------------------------------------
export PATH=$PATH:$HOME/.local/bin
while read acc
do
  if [ ! -s aligned/${acc}.gene.counts ]; then
    echo "Counting $acc"
    htseq-count -f bam -r name -s no -t exon -i gene_id \
      aligned/${acc}.Aligned.out.bam $GTF > aligned/${acc}.gene.counts
  fi
done < accessions.txt

# -----------------------------------------------------------------------------
# Step 6: build one count table for the statistics in R.
# The gene names come from the first sample, then the count column of each sample is pasted
# next to it. The htseq summary lines that start with "__" (for example unassigned reads) are
# left out. The header row has "gene" followed by the sample names from metadata.tsv.
# Result: gene.table.tsv, one row per gene and one column per sample.
# -----------------------------------------------------------------------------
if [ ! -s gene.table.tsv ]; then
  first=$(head -1 accessions.txt)
  grep -v "^__" aligned/${first}.gene.counts | cut -f1 > gene.table.tmp
  while read acc
  do
    grep -v "^__" aligned/${acc}.gene.counts | cut -f2 > col.tmp
    paste gene.table.tmp col.tmp > gene.table.tmp2
    mv gene.table.tmp2 gene.table.tmp
  done < accessions.txt
  (printf "gene"; tail -n +2 metadata.tsv | cut -f2 | sed 's/^/\t/' | tr -d '\n'; printf "\n") > gene.table.tsv
  cat gene.table.tmp >> gene.table.tsv
  rm -f gene.table.tmp col.tmp
fi
echo "Number of genes in the count table:" $(( $(wc -l < gene.table.tsv) - 1 ))
head -3 gene.table.tsv

# -----------------------------------------------------------------------------
# Next: the count table (gene.table.tsv) is analysed in R with final_project_analysis.R
# (DESeq2, cell-type marker tests, pathway tests and sensitivity checks).
# -----------------------------------------------------------------------------
