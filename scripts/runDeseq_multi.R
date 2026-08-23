# Multi-condition DESeq2 comparison: fits one DESeqDataSet on the full dataset,
# then runs every pairwise contrast between condition levels. For a simple
# two-group comparison use runDeseq.R instead.

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  getwd()
}
source(file.path(get_script_dir(), "deseq2_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
metadata_file <- args[2]
vsd_file <- args[3]
nrc_file <- args[4]
pca_file <- args[5]
DE_file <- args[6]
GroupSize <- args[7]
Type <- args[8]
RefCondition <- if (length(args) >= 9) args[9] else ""

print(metadata_file)
print("Loading the files...")
loaded <- load_counts_and_metadata(counts_file, metadata_file, Type)
data <- loaded$counts
metadata <- loaded$metadata
print("Files loaded!")

n_conditions <- nlevels(factor(metadata$condition))
if (n_conditions < 3) {
  stop("runDeseq_multi.R expects 3+ conditions (found ", n_conditions,
       "). Use runDeseq.R for a single two-group comparison.")
}

dds <- build_dds(data, metadata, GroupSize, RefCondition)

NormalizedCounts <- counts(dds, normalized = TRUE)
write.table(NormalizedCounts, nrc_file, sep = '\t', quote = FALSE)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local")
write.table(assay(vsd), vsd_file, sep = '\t', quote = FALSE)

PCA_tm(dds, pca_file)

levels_condition <- levels(colData(dds)$condition)
pairs <- combn(levels_condition, 2, simplify = FALSE)

for (pair in pairs) {
  ref_level <- pair[1]
  comp_level <- pair[2]
  Name <- paste0("condition_", comp_level, "_vs_", ref_level)

  print(paste("Running contrast:", Name))

  res <- results(dds, contrast = c("condition", comp_level, ref_level))
  res <- res %>% data.frame() %>% drop_na()

  if (nrow(res) == 0) {
    message("No genes with finite results for ", Name, " -- skipping.")
    next
  }

  VolcanoPlot(res, NormalizedCounts, Name,
              vol_file = file.path(DE_file, paste0("volcanoplot_", Name, ".pdf")),
              deg_file = file.path(DE_file, paste0("DEG_", Name, ".txt")),
              heatmap_file = file.path(DE_file, paste0("Heatmap_", Name, ".pdf")),
              padj_cutoff = 0.1, label_length = 10)
}

dir.create("Robj", showWarnings = FALSE)
save.image(file = "Robj/session.RData")
