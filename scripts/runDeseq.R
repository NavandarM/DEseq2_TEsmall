# Two-condition DESeq2 comparison. For 3+ conditions use runDeseq_multi.R instead.

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
vol_file <- args[6]
DE_file <- args[7]
GroupSize <- args[8]
Type <- args[9]
RefCondition <- if (length(args) >= 10) args[10] else ""

print(metadata_file)
print("Loading the files...")
loaded <- load_counts_and_metadata(counts_file, metadata_file, Type)
data <- loaded$counts
metadata <- loaded$metadata
print("Files loaded!")

n_conditions <- nlevels(factor(metadata$condition))
if (n_conditions != 2) {
  stop("runDeseq.R only supports exactly two conditions (found ", n_conditions,
       "). Use runDeseq_multi.R for 3+ conditions.")
}

dds <- build_dds(data, metadata, GroupSize, RefCondition)

ref_level <- levels(colData(dds)$condition)[1]
comp_level <- levels(colData(dds)$condition)[2]
Name <- paste0("condition_", comp_level, "_vs_", ref_level)

res <- results(dds, contrast = c("condition", comp_level, ref_level))
res <- res %>% data.frame() %>% drop_na()

NormalizedCounts <- counts(dds, normalized = TRUE)
write.table(NormalizedCounts, nrc_file, sep = '\t', quote = FALSE)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local")
write.table(assay(vsd), vsd_file, sep = '\t', quote = FALSE)

VolcanoPlot(res, NormalizedCounts, Name,
            vol_file = vol_file,
            deg_file = file.path(DE_file, paste0("DEG_", Name, ".txt")),
            heatmap_file = file.path(DE_file, paste0("Heatmap_", Name, ".pdf")),
            padj_cutoff = 0.1, label_length = 10)

PCA_tm(dds, pca_file)

dir.create("Robj", showWarnings = FALSE)
save.image(file = "Robj/session.RData")
