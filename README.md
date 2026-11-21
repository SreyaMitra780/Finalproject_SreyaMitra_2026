# Which retinal cell types drive gene-expression changes after Vhl loss?

Final project for BIOL4324/6324 Bioinformatics, University of Houston. Author: Sreya Mitra.

**Question.** When the Von Hippel-Lindau gene (Vhl) is lost in all cells of the developing mouse retina, which retinal cell types are the main source of the gene-expression changes?

**Data.** Bulk RNA-seq from GEO [GSE109102](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109102) (Arreola et al., 2018, JCI Insight): whole mouse retina at postnatal day 7, UBC-CreER, HiSeq 2 x 150 bp. Two control samples (SRR6464389, SRR6464390) and two Vhl-null samples (SRR6464393, SRR6464394). Raw reads were reprocessed from FASTQ; the deposited alignments were not used.

## Pipeline

| Step | Tool | Notes |
| --- | --- | --- |
| Download | SRA Toolkit 3.4.1 (`fastq-dump`) | first 2,000,000 read pairs per sample (storage limit) |
| Trim | fastp 1.4.0 | paired-end, automatic adapter detection |
| Align | STAR 2.7.11b | GRCm39 (RefSeq GCF_000001635.27); reduced sparse index to fit 15 GB memory |
| Count | htseq-count 2.0.9 | `-s no -t exon -i gene_id` |
| Test | DESeq2 | design `~ condition`, Vhl_null vs control |
| Interpret | R: Wilcoxon marker-set test, fgsea (MSigDB Hallmark, GO:BP) | 17 retinal cell-type marker sets |

## Files

```
scripts/final_project_pipeline.sh   SLURM script: download, trim, align, count (run on the Sabine cluster)
scripts/final_project_analysis.R    R analysis: DESeq2, PCA, cell-type markers, pathways, sensitivity checks
data/gene.table.tsv                 gene-by-sample count table (2M read pairs per sample)
data/metadata.tsv                   sample table
results/                            DESeq2 results and test tables (.csv)
plots/                              PCA, cell-type, and pathway figures (.png)
```

## How to run

1. On the cluster: `sbatch scripts/final_project_pipeline.sh` (produces `gene.table.tsv`).
2. In R/RStudio: set the working directory to the folder with `gene.table.tsv` and `metadata.tsv`, then run `scripts/final_project_analysis.R`. Needs DESeq2, ggplot2, msigdbr and fgsea.

## Main findings

- 37 genes changed at adjusted p < 0.1 (12 up, 25 down).
- HIF targets were higher in Vhl-null retina (Egln1, Bnip3, Slc2a1, Vegfa), and mitochondrial respiration genes were lower.
- Rod photoreceptor genes (Rho, Pde6b, Rbp3) were about 2-fold higher, and the rod marker set was significant. Vessel and glial markers showed only weak trends.
- Skeletal-muscle and brown-fat genes were lower in both Vhl-null samples, which suggests a tissue-composition difference between groups. A sensitivity check that reruns the analysis without these genes is included in `scripts/final_project_analysis.R` (section 5).

## Limitations

Two samples per group; bulk data cannot separate changed expression from changed cell proportions; first 2M read pairs per sample only; reduced STAR index; marker lists curated by hand (to be checked against Macosko et al. 2015 and Clark et al. 2019). Results are consistent with rods as the main source of change but do not prove it.

## Reference

Arreola A, et al. (2018) JCI Insight 3(4):e92193. GEO GSE109102.
