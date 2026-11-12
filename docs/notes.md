# Methods notes

Background notes on methods used in this project. The notebooks explain what was done and why; this file goes one level deeper into how some of the methods work.

## Contents

1. [Why combine all timepoints before processing](#1-why-combine-all-timepoints-before-processing)
2. [Why QC and normalization do not remove batch effects](#2-why-qc-and-normalization-do-not-remove-batch-effects)
3. [How Harmony works](#3-how-harmony-works)
4. [When to use batch correction](#4-when-to-use-batch-correction)
5. [Caveats for this dataset](#5-caveats-for-this-dataset)
6. [Practical notes](#6-practical-notes)

---

## 1. Why combine all timepoints before processing

- Concatenating the timepoints does not change any expression values. It only puts all spots into one matrix.
- The point is consistency. All timepoints share one set of highly variable genes, one PCA space, and one set of clusters, so a label like "neoblast" means the same thing at every timepoint.
- If each timepoint were processed separately, each would get its own clusters with their own numbering and boundaries, and matching them up afterwards would be guesswork.
- The working assumption: the identity of a cell type (for example, the marker genes of neoblasts) is stable over time. What changes during regeneration is mainly how many of each type there are, where they are, and their state.
- Typical workflow: combine and integrate to define cell types once, then split by timepoint to compare proportions and expression.

## 2. Why QC and normalization do not remove batch effects

- **QC** removes poor-quality spots and rarely detected genes. It works one spot at a time.
- **Normalization** removes differences in sequencing depth. Each spot is scaled to the same total count, which is a single scaling factor per spot.
- **Batch effects** come from differences between experiments: tissue sections, reagent lots, permeabilization time, slides, sequencing runs. They do not scale all genes equally. Some genes go up in one batch and others go down.
- Neither QC nor normalization compares one batch to another as a whole, and a single scaling factor cannot fix gene-specific shifts. So after normalization, the same tissue type from two batches can still sit in different places in PCA/UMAP, and clusters split by batch.

## 3. How Harmony works

**Input.** The PCA coordinates of every spot (here 50 dimensions) and the batch label of every spot (here timepoint). Harmony does not touch the expression matrix.

**Goal.** In each group of similar spots, the mix of batches should match the mix in the whole dataset. A group made almost entirely of one batch is likely a batch artifact, not a real cell type.

Harmony repeats three steps until the coordinates stop changing (usually a few rounds):

1. **Soft clustering with a diversity penalty.**
   - Spots are grouped into many clusters (about 100 by default). Each spot gets a probability of belonging to each cluster rather than a single assignment (for example 70% A, 30% B).
   - A penalty discourages clusters dominated by one batch: if a batch is already over-represented in a cluster, spots from that batch are less likely to be assigned there.
   - `theta` sets the strength of this penalty. Higher theta forces more mixing.

2. **Estimate the batch shift within each cluster.**
   - In each cluster, Harmony estimates how far each batch's spots sit from the cluster center. That offset is the batch effect for that batch in that cluster.
   - This is a weighted linear regression (coordinates = cluster center + batch offset), weighted by each spot's membership probability.
   - A ridge penalty (`lambda`) keeps the offsets from being over-estimated.

3. **Remove the shift.**
   - Each spot subtracts its batch's offset. Because a spot belongs to several clusters with different probabilities, it subtracts a weighted average of those clusters' offsets.
   - The cluster centers are kept, so real differences between cell types are preserved.

**Why correct per cluster rather than globally.** Different cell types are affected by a batch in different ways. In one batch, T cells might be shifted left and B cells shifted up. A single global shift would correct neither well. Correcting per cluster is a local, cell-type-aware correction.

**Output.** Corrected coordinates with the same shape as the input PCA (`obsm['X_pca_harmony']`). The neighbor graph, UMAP, and clustering are computed from these.

**Where the distances live.** Every spot is a point in the 50-dimensional PC space; its coordinates are its PC scores. Because PCA was run on all timepoints together, all spots share the same axes, so cluster centers and batch offsets are directly comparable vectors in that space.

## 4. When to use batch correction

- **One sample:** not needed.
- **Several samples or batches combined:** very common, close to a standard step.
- **Check first:** make a UMAP without correction. If spots clearly separate by batch, correct; if they already mix, correction may not be necessary.
- **Common methods:** Harmony (fast and simple), BBKNN, scVI, Seurat integration.

## 5. Caveats for this dataset

- **Batch and biology are confounded.** Here the batch variable (timepoint) is also the variable of biological interest. Harmony cannot tell technical from biological differences, so it may also smooth out real changes over time. In particular, a cell state that exists at only one timepoint can be pulled toward other timepoints.
- **What is safe to use.** Cell-type proportions per timepoint are still meaningful, because they are counts of spots within clusters. Changes in gene expression over time should be tested on the normalized expression matrix (or raw counts), never on the Harmony coordinates.
- **Finer batch structure.** Each timepoint consists of four Visium capture areas (`library_id`), and each capture area holds several tissue sections (`Section_id`). This analysis corrects only at the timepoint level. Correcting by `library_id` is an option to explore, but it would also partly remove timepoint differences.
- **Spots are not cells.** Each Visium spot is about 55 µm across and usually covers several cells, so clusters describe tissue composition, not pure single-cell types.

## 6. Practical notes

- `sc.external.pp.harmony_integrate` mis-handles the output shape with newer harmonypy versions (PyTorch backend), so the notebook calls `harmonypy.run_harmony` directly, moves the result off GPU/MPS if needed, and transposes it to (spots × PCs) if necessary.
- pandas 3.0 removed `DataFrame.applymap`; use `DataFrame.map` instead.
- `sc.settings.set_figure_params` is deprecated; use `sc.set_figure_params`.
