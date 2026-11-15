import celldex
import scranpy
import scrnaseq
import seaborn as sns
import singler

# 1. Fetch a scRNA-seq dataset
zeisel_brain = scrnaseq.fetch_dataset(
    "zeisel-brain-2015", "2023-12-14", realize_assays=True
)
print(zeisel_brain)

# 2. Run full scran analysis with defaults
results = scranpy.analyze(zeisel_brain, num_threads=2)

# 3. Fetch the reference from celldex
mouse_rnaseq_ref = celldex.fetch_reference(
    "mouse_rnaseq", "2024-02-26", realize_assays=True
)

# 4. Run SingleR
annotation_results = singler.annotate_single(
    test_data=zeisel_brain,
    ref_data=mouse_rnaseq_ref,
    ref_labels="label.main",
    num_threads=4,
)
print(annotation_results)

# 5. Attach cell type annotations on the original dataset
zeisel_brain.col_data.set_column(
    "label", annotation_results["best"], in_place=True
)

# 6. Plot t-SNE coordinates from scranpy analysis with singleR annotated clusters
to_keep = [i for i, x in enumerate(results.rna_qc_filter) if x is True]
filtered_matches = [annotation_results["best"][i] for i in to_keep]
sns.scatterplot(
    x=results.tsne[0],
    y=results.tsne[1],
    hue=filtered_matches,
    palette="Paired",
)

