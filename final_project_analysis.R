# =============================================================================
# Final project: which retinal cell types drive gene-expression changes
# after Vhl loss in the developing mouse retina?
# Data: GEO GSE109102 (Arreola et al. 2018), bulk RNA-seq, P7 mouse retina
# Input:  gene.table.tsv (htseq-count table from final_project.sh), metadata.tsv
# Run:    in RStudio, set the working directory to the folder with those files.
# =============================================================================

# setwd("~/Desktop/Second year PhD 2026-2027/final_project")   # <- set to your folder

library(DESeq2)
library(ggplot2)
library(msigdbr)
library(fgsea)

# ---- 1. Load counts and sample table -----------------------------------------
counts_matrix <- read.table("gene.table.tsv", header = TRUE, sep = "\t",
                            row.names = 1, check.names = FALSE)
meta <- read.table("metadata.tsv", header = TRUE, sep = "\t")
colData <- data.frame(condition = factor(meta$condition, levels = c("control", "Vhl_null")),
                      row.names = meta$sample)

# ---- 2. Differential expression (DESeq2) ---------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = colData, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Vhl_null", "control"))
summary(res)
write.csv(as.data.frame(res), "deseq2_results.csv")

# PCA of variance-stabilized counts
vsd <- vst(dds, blind = FALSE)
p_pca <- plotPCA(vsd, intgroup = "condition")
ggsave("pca.png", p_pca, width = 6, height = 5)

# ---- 3. Cell-type marker test --------------------------------------------------
# Marker lists were written from canonical markers; verify against published
# retina atlases (Macosko et al. 2015; Clark et al. 2019) before citing.
markers_full <- list(
  Rod_photoreceptor    = c("Rho", "Nrl", "Gnat1", "Pde6a", "Pde6b", "Sag", "Nr2e3", "Cnga1"),
  Cone_photoreceptor   = c("Opn1sw", "Opn1mw", "Arr3", "Gnat2", "Gngt2", "Pde6c", "Pde6h", "Thrb"),
  Rod_bipolar          = c("Prkca", "Pcp2", "Grm6", "Trpm1", "Cabp5"),
  Cone_bipolar         = c("Grik1", "Vsx1", "Syt2", "Lhx4"),
  Amacrine_GABA        = c("Gad1", "Gad2", "Slc32a1", "Tfap2a", "Tfap2b"),
  Amacrine_glycine     = c("Slc6a9", "Slc6a5"),
  Amacrine_starburst   = c("Chat", "Slc18a3", "Megf10"),
  Horizontal           = c("Onecut1", "Onecut2", "Lhx1", "Prox1", "Ptf1a", "Calb1"),
  Ganglion             = c("Pou4f1", "Pou4f2", "Rbpms", "Thy1", "Sncg", "Nefl", "Slc17a6", "Isl2"),
  Muller_glia          = c("Glul", "Rlbp1", "Slc1a3", "Aqp4", "Sox9", "Apoe", "Clu", "Dkk3", "Vim"),
  Retinal_progenitor   = c("Mki67", "Top2a", "Pcna", "Ccnd1", "Hes1", "Hes5", "Fgf15"),
  Astrocyte            = c("Gfap", "Pax2", "S100b", "Aldh1l1"),
  Microglia_macrophage = c("Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Csf1r", "Hexb", "C1qa", "C1qb"),
  Endothelial          = c("Pecam1", "Cldn5", "Kdr", "Flt1", "Cdh5", "Tie1", "Esam", "Emcn"),
  Pericyte_mural       = c("Pdgfrb", "Rgs5", "Cspg4", "Acta2", "Kcnj8", "Abcc9", "Des"),
  Fibroblast           = c("Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"),
  RPE                  = c("Rpe65", "Best1", "Tyr", "Pmel", "Lrat", "Rdh5", "Mitf")
)

# For each marker set: Wilcoxon test of the DESeq2 statistic of its genes vs all other genes
test_markers <- function(res, marker_list) {
  d <- as.data.frame(res)
  d <- d[!is.na(d$stat), ]
  out <- do.call(rbind, lapply(names(marker_list), function(ct) {
    g <- intersect(marker_list[[ct]], rownames(d))
    if (length(g) < 2) return(NULL)
    s_in  <- d[g, "stat"]
    s_out <- d[setdiff(rownames(d), g), "stat"]
    data.frame(cell_type = ct, markers_found = length(g),
               mean_log2FC = mean(d[g, "log2FoldChange"]),
               median_stat = median(s_in),
               p_value = wilcox.test(s_in, s_out)$p.value)
  }))
  out$adj_p <- p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}

ct_full <- test_markers(res, markers_full)
print(ct_full)
write.csv(ct_full, "cell_type_enrichment_full.csv", row.names = FALSE)

p_full <- ggplot(ct_full, aes(x = reorder(cell_type, median_stat), y = median_stat, fill = median_stat > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
  labs(x = "", y = "Median DESeq2 statistic of marker genes (Vhl-null vs control)",
       title = "Cell-type marker genes in Vhl-null retina") +
  theme_bw()
ggsave("cell_type_markers_full.png", p_full, width = 8, height = 6.5)

# ---- 4. Pathway tests (fgsea on MSigDB Hallmark and GO:BP) ----------------------
run_fgsea <- function(res, use_go = TRUE) {
  d <- as.data.frame(res)
  d <- d[!is.na(d$stat), ]
  ranks <- sort(setNames(d$stat, rownames(d)), decreasing = TRUE)
  gs <- msigdbr(species = "Mus musculus", collection = "H")
  if (use_go) gs <- rbind(gs, msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP"))
  pathways <- lapply(split(gs$gene_symbol, gs$gs_name), unique)
  set.seed(1)
  fg <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500)
  as.data.frame(fg[order(fg$pval), c("pathway", "size", "NES", "pval", "padj")])
}

themes <- "ANGIOGENESIS|HYPOXIA|GLYCOLYSIS|NOTCH|TGF_BETA|SMOOTH_MUSCLE|BLOOD_VESSEL|VASCULATURE|LIPID|FATTY_ACID|CHOLESTEROL|STEROL"

fg_hall <- run_fgsea(res, use_go = FALSE)          # the 50 Hallmark sets only
print(head(fg_hall, 15))
write.csv(fg_hall, "pathways_hallmark.csv", row.names = FALSE)

fg_all <- run_fgsea(res, use_go = TRUE)            # Hallmark + GO:BP, theme-selected
sel <- fg_all[grepl(themes, fg_all$pathway), ]
print(head(sel, 25))
write.csv(sel, "pathways_selected.csv", row.names = FALSE)

p_path <- ggplot(head(sel, 20), aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
  labs(x = "", y = "Normalized enrichment score (Vhl-null vs control)",
       title = "Selected pathways in Vhl-null retina") +
  theme_bw()
ggsave("pathways_selected.png", p_path, width = 10, height = 7)

# One enrichment plot and one top-genes plot per chosen pathway
chosen <- c("HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
            "HALLMARK_NOTCH_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING",
            "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS")
hall <- msigdbr(species = "Mus musculus", collection = "H")
hall_sets <- lapply(split(hall$gene_symbol, hall$gs_name), unique)
d_all <- as.data.frame(res)
d_all <- d_all[!is.na(d_all$stat), ]
ranks <- sort(setNames(d_all$stat, rownames(d_all)), decreasing = TRUE)
dir.create("pathway_plots", showWarnings = FALSE)
for (pw in chosen) {
  p1 <- plotEnrichment(hall_sets[[pw]], ranks) + labs(title = pw)
  ggsave(file.path("pathway_plots", paste0(pw, "_enrichment.png")), p1, width = 6, height = 4)
  gd <- d_all[intersect(hall_sets[[pw]], rownames(d_all)), ]
  gd$gene <- rownames(gd)
  gd <- gd[order(abs(gd$stat), decreasing = TRUE), ]
  gd <- gd[seq_len(min(20, nrow(gd))), ]
  p2 <- ggplot(gd, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
    labs(x = "", y = "log2 fold change (Vhl-null vs control)", title = paste(pw, "- top genes")) +
    theme_bw()
  ggsave(file.path("pathway_plots", paste0(pw, "_genes.png")), p2, width = 6, height = 5)
}

# ---- 5. Sensitivity check A: remove muscle and brown-fat genes ------------------
contam <- c("Myl1", "Acta1", "Tnni2", "Tnnc2", "Myh8", "Mb", "Ckm", "Tpm2",
            "Myl4", "Smpx", "Fabp3", "Fabp4", "Ucp1", "Cox8b")
print(round(counts(dds, normalized = TRUE)[intersect(contam, rownames(dds)), ], 1))

keep <- !(rownames(dds) %in% contam)
cd <- data.frame(condition = dds$condition, row.names = colnames(dds))
dds_clean <- DESeqDataSetFromMatrix(countData = counts(dds)[keep, ], colData = cd, design = ~ condition)
dds_clean <- DESeq(dds_clean)
res_clean <- results(dds_clean, contrast = c("condition", "Vhl_null", "control"))
summary(res_clean)
print(res_clean[c("Rho", "Pde6b", "Rbp3", "Vegfa", "Slc2a1", "Bnip3", "Egln1"),
                c("log2FoldChange", "pvalue", "padj")])
ct_clean <- test_markers(res_clean, markers_full)
print(ct_clean)
write.csv(as.data.frame(res_clean), "deseq2_results_no_muscle_fat.csv")
write.csv(ct_clean, "cell_type_enrichment_no_muscle_fat.csv", row.names = FALSE)

# ---- 6. Sensitivity check B: remove rod photoreceptor genes ---------------------
rod_genes <- unique(c(markers_full$Rod_photoreceptor, "Gngt1", "Pde6g", "Rcvrn", "Crx", "Rom1", "Prph2",
                      "Cngb1", "Guca1a", "Guca1b", "Rp1", "Abca4", "Pdc", "Rbp3", "Aipl1", "Rs1", "Tulp1"))
keep <- !(rownames(dds) %in% rod_genes)
cat("Rod genes removed:", sum(!keep), "\n")
dds_norod <- DESeqDataSetFromMatrix(countData = counts(dds)[keep, ], colData = cd, design = ~ condition)
dds_norod <- DESeq(dds_norod)
res_norod <- results(dds_norod, contrast = c("condition", "Vhl_null", "control"))
summary(res_norod)
write.csv(as.data.frame(res_norod), "deseq2_results_no_rod_genes.csv")
ct_norod <- test_markers(res_norod, markers_full[names(markers_full) != "Rod_photoreceptor"])
print(ct_norod)
write.csv(ct_norod, "cell_type_enrichment_no_rod_genes.csv", row.names = FALSE)

# ---- 7. Link results: changed genes -> pathways -> cell-type markers --------------
d <- as.data.frame(res)
d$gene <- rownames(d)
changed <- d[!is.na(d$pvalue) & d$pvalue < 0.01 & abs(d$log2FoldChange) > 0.5, ]
changed <- changed[order(changed$pvalue), ]
lookup <- function(g, sets) paste(names(sets)[sapply(sets, function(s) g %in% s)], collapse = "; ")
changed$hallmark_pathways <- sapply(changed$gene, lookup, sets = hall_sets)
changed$cell_type_marker  <- sapply(changed$gene, lookup, sets = markers_full)
linked <- changed[, c("gene", "baseMean", "log2FoldChange", "pvalue", "padj", "hallmark_pathways", "cell_type_marker")]
print(head(linked, 40))
write.csv(linked, "linked_genes_pathways_cell_types.csv", row.names = FALSE)
nrow(linked)

sessionInfo()
