# Data

- `smed_blast.csv`: SMED gene ID to human gene symbol lookup (BLAST). Not tracked in git.
- `Scanpy_Visium/`: Visium sections for six timepoints (0 hpa, 6 hpa, 12 hpa, 1 dpa, 3 dpa, 7 dpa). Not tracked in git.

Expected files in `Scanpy_Visium/`:

```
Visium_0hpa.h5ad
Visium_6hpa.h5ad
Visium_12hpa.h5ad
Visium_1dpa.h5ad
Visium_3dpa.h5ad
Visium_7dpa.h5ad
```

## Source

Cui et al., 2023, Nature Communications (https://www.nature.com/articles/s41467-023-39016-0).

- Raw sequencing data: GSA, accession CRA007941 (project PRJCA011425)
- H&E images and processed data: OMIX, accessions OMIX003867 and OMIX003889
- Reference transcriptome (smed_20140614): GEO, accession GSE72389

## Download

The `.h5ad` files come directly from OMIX003867, file `OMIX003867-03` (`Scanpy_Visium.zip`, about 760 MB), with no further conversion:
https://ngdc.cncb.ac.cn/omix/release/OMIX003867

Unzip it so the six `.h5ad` files sit in `data/Scanpy_Visium/`.
