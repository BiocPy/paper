library(scRNAseq)
library(scrapper)
library(celldex)
library(SingleR)

# 1. Fetch a scRNA-seq dataset
sce <- fetchDataset("zeisel-brain-2015", "2023-12-14", realize.assays = TRUE)
print(sce)

# 2. Run full scran analysis with defaults
res <- analyze(sce, num.threads = 2)

# 3. Fetch the reference from celldex
mouse_rnaseq_ref <- fetchReference("mouse_rnaseq", "2024-02-26", realize.assays = TRUE)

# 4. Run SingleR
cell_labels <- SingleR(
    test = assay(sce, "counts"),
    ref = mouse_rnaseq_ref,
    labels = mouse_rnaseq_ref$label.main,
    num.threads = 4
)
table(cell_labels$labels)

# 5. Attach cell type annotations on the original dataset
colData(sce)["label"] <- cell_labels$labels

