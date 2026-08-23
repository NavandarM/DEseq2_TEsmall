suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(DESeq2)
})

# Subset/reorder a counts data frame's sample columns to match metadata$sample order.
reorder_columns <- function(counts_df, metadata_df, fixed_cols = c("fid", "ftype")) {
  fixed_cols <- intersect(fixed_cols, colnames(counts_df))
  matching_samples <- intersect(metadata_df$sample, colnames(counts_df))
  counts_df %>% select(all_of(c(fixed_cols, matching_samples)))
}

# Load counts + metadata and guarantee colnames(counts) == rownames(metadata) in the
# same order, dropping metadata rows for samples absent from the counts file.
# (DESeqDataSetFromMatrix does not verify this itself when colData has no rownames,
# so a mismatch here silently reassigns conditions to the wrong samples.)
load_counts_and_metadata <- function(counts_file, metadata_file, Type) {
  Undata <- read.delim(counts_file, check.names = FALSE)
  metadata <- read.delim(metadata_file)

  if ("smallRNA" %in% Type) {
    data <- reorder_columns(Undata, metadata)
  } else {
    matching_samples <- intersect(metadata$sample, colnames(Undata))
    data <- Undata[, c(colnames(Undata)[1], matching_samples), drop = FALSE]
  }
  colnames(data)[1] <- "Ids"

  if ("ftype" %in% names(data)) data$ftype <- NULL

  rownames(data) <- data$Ids
  data$Ids <- NULL

  sample_cols <- colnames(data)
  missing_from_counts <- setdiff(metadata$sample, sample_cols)
  if (length(missing_from_counts) > 0) {
    warning("Samples listed in metadata but absent from the counts file were dropped: ",
            paste(missing_from_counts, collapse = ", "))
  }

  metadata <- metadata[match(sample_cols, metadata$sample), , drop = FALSE]
  rownames(metadata) <- metadata$sample

  if (!identical(colnames(data), rownames(metadata))) {
    stop("Could not align counts columns with metadata samples.")
  }

  list(counts = data, metadata = metadata)
}

# Build and fit the DESeqDataSet shared by every comparison.
build_dds <- function(counts, metadata, group_size, ref_condition = "") {
  metadata$sample <- as.factor(metadata$sample)
  metadata$condition <- factor(metadata$condition)

  if (nzchar(ref_condition)) {
    if (!ref_condition %in% levels(metadata$condition)) {
      stop("RefCondition '", ref_condition, "' not found among conditions: ",
           paste(levels(metadata$condition), collapse = ", "))
    }
    metadata$condition <- relevel(metadata$condition, ref = ref_condition)
  }

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1 + condition)

  smallestGroupSize <- as.integer(group_size)
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep, ]

  DESeq(dds, fitType = "local")
}

PCA_tm <- function(dds, pca_file, groups = "condition", trans_func = varianceStabilizingTransformation) {
  color_var <- if (length(groups) > 0) sym(groups[1]) else NULL

  PCA <- dds %>% trans_func() %>% plotPCA(intgroup = groups, returnData = TRUE)
  percentVar <- round(100 * attr(PCA, "percentVar"), 1)

  n_groups <- length(unique(PCA[[groups[1]]]))
  # A 2-color manual palette silently breaks (recycled/missing colors) once there
  # are 3+ groups, so size the palette to however many groups are actually present.
  palette <- if (n_groups <= 2) {
    c("dodgerblue", "darkorange")
  } else {
    scales::hue_pal()(n_groups)
  }

  p <- ggplot(PCA, aes(x = PC1, y = PC2, color = {{color_var}})) +
    geom_point(size = 6, alpha = 0.5) +
    geom_text_repel(aes(label = name), size = 5, show.legend = FALSE) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    scale_color_manual(values = palette) +
    theme_bw() +
    theme(axis.text = element_text(size = 17),
          axis.title = element_text(size = 17))

  ggsave(filename = pca_file, plot = p, width = 12, height = 12, dpi = 300, units = "in")
  p
}

VolcanoPlot <- function(res, NormalizedCounts, Name, vol_file, deg_file, heatmap_file,
                        padj_cutoff = 0.05, label_length = 20) {

  Filtered <- res[res$padj < padj_cutoff, ]
  Up <- Filtered[Filtered$log2FoldChange > 0, ]
  Down <- Filtered[Filtered$log2FoldChange < 0, ]

  write.table(Filtered, deg_file, sep = "\t", quote = FALSE)

  # min(label_length, nrow(...)) with seq_len avoids indexing past the last row
  # (which previously padded top_genes with NA rows whenever fewer than
  # label_length genes were up/down).
  n_up <- min(label_length, nrow(Up))
  n_down <- min(label_length, nrow(Down))

  top_Up <- Up[order(-Up$log2FoldChange), ][seq_len(n_up), , drop = FALSE]
  top_Up$Names <- sub("^([^:]+:[^:]+):.*", "\\1", row.names(top_Up))

  top_Down <- Down[order(Down$log2FoldChange), ][seq_len(n_down), , drop = FALSE]
  top_Down$Names <- sub("^([^:]+:[^:]+):.*", "\\1", row.names(top_Down))

  top_genes <- rbind(top_Up, top_Down)

  legend_labels <- c("Non-Significant",
    paste("\nSignificant\nUp:", nrow(Up), "\nDown:", nrow(Down)))

  # Symmetric axis limits must be based on the largest-magnitude fold change in
  # EITHER direction; using min() twice (as before) clips whichever tail has the
  # larger absolute value.
  lfc_limit <- max(abs(res$log2FoldChange), na.rm = TRUE)

  p <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < padj_cutoff), size = 1.5) +
    scale_color_manual(name = "",
                       values = c('FALSE' = '#7f7f7f', 'TRUE' = '#f62728'),
                       labels = legend_labels) +
    xlim(-lfc_limit, lfc_limit) +
    geom_text_repel(data = top_genes, aes(label = Names), size = 4, max.overlaps = 40) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(1, 'cm')) +
    guides(color = guide_legend(byrow = TRUE))

  ggsave(filename = vol_file, plot = p, width = 12, height = 12, dpi = 300, units = "in")

  common_ids <- intersect(rownames(Filtered), rownames(NormalizedCounts))
  Filtered_sub <- Filtered[common_ids, ]
  Norm_sub <- NormalizedCounts[common_ids, ]
  Merged <- cbind(Filtered_sub, Norm_sub)

  if (nrow(Merged) >= 2) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      install.packages("pheatmap", repos = "https://cloud.r-project.org/")
    }
    library(pheatmap)

    my_colors <- colorRampPalette(c("brown", "white", "blue"))(40)
    # Count columns start right after Filtered_sub's own columns, rather than a
    # hardcoded "7" that silently breaks if res's column layout ever changes.
    count_cols <- (ncol(Filtered_sub) + 1):ncol(Merged)
    pheatmap(Merged[, count_cols], filename = heatmap_file, scale = "row",
             color = my_colors, fontsize_row = 8, cluster_cols = FALSE)
  } else {
    message("Fewer than 2 significant genes for ", Name, " -- skipping heatmap.")
  }

  invisible(p)
}
