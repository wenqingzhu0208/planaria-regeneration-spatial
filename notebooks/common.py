"""Shared loading helpers for notebooks 03-05.

Loads the largest Visium section per timepoint (the same sections used in
02_wound_region_neoblasts.ipynb) and attaches the wound labels saved there.
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc

DATA_DIR = '../data/Scanpy_Visium'
PROCESSED_DIR = '../data/processed'

TP_INFO = {
    '0hpa':  ('Visium_0hpa.h5ad',  0),
    '6hpa':  ('Visium_6hpa.h5ad',  6),
    '12hpa': ('Visium_12hpa.h5ad', 12),
    '1dpa':  ('Visium_1dpa.h5ad',  24),
    '3dpa':  ('Visium_3dpa.h5ad',  72),
    '7dpa':  ('Visium_7dpa.h5ad',  168),
}
TP_ORDER = list(TP_INFO)
HOURS = [h for _, h in TP_INFO.values()]
TICK_LBL = ['0h', '6h', '12h', '1d', '3d', '7d']

CELL_TYPES = ['Epidermal', 'Gut', 'Muscle', 'Neoblast',
              'Neuronal', 'Parenchymal', 'Secretory']
CT_COLS = [f'cell2location_{ct}' for ct in CELL_TYPES]


def load_max_section(fname):
    """Load the section with the most spots; attach scaled spatial coords."""
    # Open in backed mode so only the selected section is read into memory
    a = sc.read_h5ad(os.path.join(DATA_DIR, fname), backed='r')
    sec_id = a.obs['Section_id'].value_counts().idxmax()
    sec = a[(a.obs['Section_id'] == sec_id).values].to_memory()
    lib_id = sec.obs['library_id'].iloc[0]
    sf = a.uns['spatial'][lib_id]['scalefactors']['tissue_lowres_scalef']
    sec.obsm['spatial_scaled'] = sec.obsm['spatial'] * sf
    sec.uns.pop('spatial', None)  # H&E images (large) are not needed
    a.file.close()
    return sec


def load_sections():
    """Return {timepoint: AnnData} with a boolean obs['wound_proximal'] column.

    Wound labels come from 02_wound_region_neoblasts.ipynb
    (data/processed/02_wound_regions.csv); run that notebook first.
    """
    wound = pd.read_csv(os.path.join(PROCESSED_DIR, '02_wound_regions.csv'))
    sections = {}
    for tp, (fname, _) in TP_INFO.items():
        sec = load_max_section(fname)
        labels = wound[wound['timepoint'] == tp].set_index('spot')['wound_proximal']
        sec.obs['wound_proximal'] = labels.reindex(sec.obs_names).values
        assert sec.obs['wound_proximal'].notna().all(), f'missing wound labels for {tp}'
        sec.obs['wound_proximal'] = sec.obs['wound_proximal'].astype(bool)
        sections[tp] = sec
    return sections


def cell_props(sec):
    """Per-spot cell-type proportions from the cell2location estimates."""
    est = sec.obs[CT_COLS]
    props = est.div(est.sum(axis=1), axis=0)
    props.columns = CELL_TYPES
    return props


def load_blast():
    """SMED ID <-> human gene symbol lookups (BLAST best hit)."""
    blast = pd.read_csv('../data/smed_blast.csv').dropna(subset=['gene_symbol'])
    blast = blast.drop_duplicates('smed_id')
    id2sym = dict(zip(blast['smed_id'], blast['gene_symbol']))
    sym2ids = blast.groupby('gene_symbol')['smed_id'].apply(list).to_dict()
    return id2sym, sym2ids
