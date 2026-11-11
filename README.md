# Spatial dynamics of neural regeneration in planarians

Re-analysis of public Visium spatial transcriptomics data from a planarian regeneration time course (0 hpa to 7 dpa; [Cui et al., 2023](https://www.nature.com/articles/s41467-023-39016-0)). The project follows one question chain: which cell types change after amputation, where neural tissue is lost and rebuilt, how it is rebuilt, and whether the rebuilt tissue is molecularly mature.

**Main finding:** neural tissue at the wound is spatially rebuilt by 7 dpa, but its gene expression has only partly returned to the uninjured state. Spatial completion comes before molecular completion.

## Questions

1. Which cell types change the most during regeneration?
2. Are neurons lost and rebuilt specifically at the wound?
3. How are neurons repopulated at the wound?
4. Is the regenerated tissue at 7 dpa molecularly mature?

## Key results

**Q1. Neurons are the most dynamic differentiated cell type.** All seven cell types change over time; neurons follow a V-shaped trajectory (high at 12 hpa, lowest at 3 dpa, recovered by 7 dpa).

![Q1](figures/q1_celltype_kinetics.png)

**Q2. Neuronal loss and recovery are wound-specific.** Neurons are depleted at the wound at 12 hpa and 1 dpa, then exceed wound-distal levels at 7 dpa (0.34 vs 0.23).

![Q2](figures/q2_wound_vs_distal.png)

**Q3. Spatial statistics best fit a local burst of neuron production followed by self-organization.** At 7 dpa, neurons are the most spatially clustered of the time course and most strongly separated from neoblasts.

![Q3](figures/q3_morans_i.png)

**Q4. The regenerated tissue is only partly mature.** About 35% of the baseline wound-region gene signature is recovered at 7 dpa; genes specific to regeneration are enriched for muscle and extracellular matrix, not for cell cycle or stem-cell programs.

![Q4](figures/q4_enrichment.png)

## Methods

- **Cell-type proportions:** cell2location estimates (provided with the dataset) divided by their per-spot sum.
- **Anatomical alignment:** pharynx center located from FOXA1-ortholog expression; wound-proximal spots defined by distance from the cut edges.
- **Statistics:** Kruskal-Wallis (Q1), two-sided Mann-Whitney U (Q2), Spearman correlation and univariate/bivariate Moran's I on a k-nearest-neighbor graph with permutation tests and spatial block bootstrap (Q3), Wilcoxon differential expression with Benjamini-Hochberg correction and Fisher's exact test for gene-set enrichment (Q4).
- **Gene annotation:** planarian (SMED) genes mapped to human gene symbols by BLAST best hit (`data/smed_blast.csv`, not included in this repository; see Data source).

## Limitations

- **No biological replicates.** One section per timepoint, so p-values reflect variation between spots within one section, not between animals.
- **Cell-type proportions are model estimates**, not direct measurements; biases in cell2location carry into every downstream result.
- **The wound mask changes meaning over time.** At 0 hpa it covers the intact head and tail tips, used here as an anatomical baseline.
- **Fixed thresholds** (wound distance, marker percentile) were not varied in a sensitivity analysis.
- **Ortholog mapping by BLAST best hit** is approximate, and planarian genes without a human match are excluded from enrichment analysis.

## How to run

```bash
uv sync
# put the .h5ad files in data/Scanpy_Visium/ (see data/README.md)
uv run jupyter lab notebooks/analysis.ipynb
```

Running the notebook also writes the figures to `figures/`.

## Repository layout

```
notebooks/analysis.ipynb   main analysis
data/                      inputs (large files not tracked)
figures/                   figures exported by the notebook
```

## Data source

Cui et al., 2023, *Nature Communications*. Processed Visium data: OMIX003867. The SMED-to-human BLAST table (`data/smed_blast.csv`) is not included in this repository yet; it is needed to run the notebook.
