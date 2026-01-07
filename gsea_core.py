"""
Core GSEA logic for the GSEA Helper App.

This module is responsible for:
- Preparing a ranked gene list from a DE results DataFrame
- Running GSEA (preranked) using gseapy
- Returning a tidy pandas DataFrame with useful summary columns
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
# NOTE: we import gseapy lazily inside run_gsea to avoid test-time dependency


# Simple mapping from (organism, collection) -> gseapy gene_sets name.
# You can adjust these if you prefer other libraries.
GENESET_MAP: Dict[Tuple[str, str], str] = {
    ("Homo sapiens", "HALLMARK"): "MSigDB_Hallmark_2020",
    ("Homo sapiens", "KEGG"): "KEGG_2016",
    ("Mus musculus", "HALLMARK"): "MSigDB_Hallmark_2020",  # via orthologs; OK for demo
    ("Mus musculus", "KEGG"): "KEGG_2019_Mouse",
}


@dataclass
class GSEAConfig:
    """Configuration for a single GSEA run."""
    organism: str
    gene_set_collection: str
    min_size: int = 15
    max_size: int = 500
    permutation_num: int = 100
    seed: int = 42

    @property
    def gene_sets(self) -> str:
        key = (self.organism, self.gene_set_collection)
        if key not in GENESET_MAP:
            available = ", ".join(
                f"{org} / {coll}" for (org, coll) in sorted(GENESET_MAP.keys())
            )
            raise ValueError(
                f"No gene set mapping for organism={self.organism!r}, "
                f"collection={self.gene_set_collection!r}. "
                f"Available combinations: {available}"
            )
        return GENESET_MAP[key]


def prepare_ranked_list(
    df: pd.DataFrame,
    gene_col: str,
    score_col: str,
    ascending: bool = False,
) -> pd.DataFrame:
    """
    Prepare a ranked list for GSEA (prerank).

    Parameters
    ----------
    df : DataFrame
        DE results table.
    gene_col : str
        Column name with gene identifiers (e.g. gene symbols).
    score_col : str
        Column name with ranking scores (e.g. log2FC, stat).
    ascending : bool
        Sort scores in ascending order. Default False (highest first).

    Returns
    -------
    rnk : DataFrame
        DataFrame indexed by gene, single column 'score', sorted.
    """
    if gene_col not in df.columns:
        raise ValueError(f"Gene column {gene_col!r} not found in DataFrame.")
    if score_col not in df.columns:
        raise ValueError(f"Score column {score_col!r} not found in DataFrame.")

    # Keep only needed columns, drop missing scores
    rnk = (
        df[[gene_col, score_col]]
        .dropna(subset=[score_col])
        .copy()
    )

    # Convert to string and uppercase (recommended for Enrichr / MSigDB libs)
    rnk[gene_col] = rnk[gene_col].astype(str).str.upper()

    # Drop duplicate genes, keep the first occurrence
    rnk = rnk.drop_duplicates(subset=[gene_col], keep="first")

    # Rename and sort
    rnk = rnk.rename(columns={gene_col: "gene", score_col: "score"})
    rnk = rnk.sort_values("score", ascending=ascending)

    # GSEApy expects gene as index or first column
    rnk = rnk.set_index("gene")

    return rnk


def _find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    """Return the first column in candidates that exists in df.columns."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _parse_tag_fraction(tag_value: str) -> Tuple[float, float, float]:
    """
    Parse the 'Tag %' column from gseapy (like '16/52').

    Returns
    -------
    gene_ratio : float
        hits / set_size
    hits : float
    set_size : float
    """
    try:
        # Some versions give '16/52', others may have spaces. Be generous.
        parts = str(tag_value).split("/")
        if len(parts) != 2:
            return np.nan, np.nan, np.nan
        hits = float(parts[0])
        size = float(parts[1])
        if size == 0:
            return np.nan, hits, size
        return hits / size, hits, size
    except Exception:
        return np.nan, np.nan, np.nan


def tidy_gsea_results(pre_res: gp.prerank) -> pd.DataFrame:
    """
    Convert gseapy prerank result into a tidy DataFrame with extra columns.

    Added columns:
    - pathway      : pathway / term name
    - NES          : normalized enrichment score
    - fdr          : FDR q-value
    - direction    : 'upregulated' / 'downregulated' based on NES
    - gene_ratio   : leading-edge hits / gene set size (from 'Tag %')
    - num_genes    : number of leading-edge genes
    - geneset_size : total number of genes in the gene set
    - lead_genes   : list string of leading-edge genes
    """
    # gseapy stores a nice summary table in res2d
    df = pre_res.res2d.reset_index(drop=True).copy()

    # Harmonize column names across gseapy versions
    term_col = _find_col(df, ["Term", "term"])
    nes_col = _find_col(df, ["NES", "nes"])
    fdr_col = _find_col(df, ["FDR q-val", "FDR_q-val", "fdr", "FDR"])
    tag_col = _find_col(df, ["Tag %", "tag %", "tag_percent", "Tag_percent"])
    lead_col = _find_col(df, ["Lead_genes", "lead_genes"])
    name_col = _find_col(df, ["Name", "name"])

    if term_col is None or nes_col is None:
        raise RuntimeError(
            f"Unexpected GSEApy result columns: {list(df.columns)}. "
            "Cannot find 'Term' / 'NES' columns."
        )

    # Rename to standard names we will use everywhere
    rename_map = {
        term_col: "pathway",
        nes_col: "NES",
    }
    if fdr_col:
        rename_map[fdr_col] = "fdr"
    if tag_col:
        rename_map[tag_col] = "Tag_fraction"
    if lead_col:
        rename_map[lead_col] = "lead_genes"
    if name_col:
        rename_map[name_col] = "comparison_name"

    df = df.rename(columns=rename_map)

    # Direction based on NES
    df["direction"] = np.where(df["NES"] > 0, "upregulated", "downregulated")

    # Gene ratio and sizes from Tag %
    if "Tag_fraction" in df.columns:
        parsed = df["Tag_fraction"].apply(_parse_tag_fraction)
        df["gene_ratio"], df["num_genes"], df["geneset_size"] = zip(*parsed)
    else:
        df["gene_ratio"] = np.nan
        df["num_genes"] = np.nan
        df["geneset_size"] = np.nan

    # If we have lead_genes, also compute num_genes from that as a fallback
    if "lead_genes" in df.columns:
        # gseapy uses ';' as separator in many versions
        df["num_genes"] = df["lead_genes"].astype(str).str.split("[,;]").str.len()

    # Fallback if fdr is missing: use nominal p-val if present
    if "fdr" not in df.columns:
        pval_col = _find_col(df, ["NOM p-val", "pval", "Pval", "P-value"])
        if pval_col:
            df = df.rename(columns={pval_col: "fdr"})
        else:
            df["fdr"] = np.nan

    # Sort by NES by default (can be overridden later)
    df = df.sort_values("NES", ascending=False)

    # Keep original columns too, but ensure our main ones are present
    # For convenience, reorder to put the important columns first
    preferred_order = [
        "pathway",
        "NES",
        "fdr",
        "direction",
        "gene_ratio",
        "num_genes",
        "geneset_size",
        "lead_genes",
    ]
    cols = [c for c in preferred_order if c in df.columns] + [
        c for c in df.columns if c not in preferred_order
    ]
    df = df[cols]

    return df

def run_gsea(
    df: pd.DataFrame,
    gene_col: str,
    score_col: str,
    config: GSEAConfig,
) -> pd.DataFrame:
    """
    High-level helper: prepare ranked list, run prerank GSEA, return tidy DataFrame.

    Parameters
    ----------
    df : DataFrame
        DE results.
    gene_col : str
        Column with gene identifiers.
    score_col : str
        Column with ranking metric.
    config : GSEAConfig
        Configuration (organism, collection, permutations, etc.).

    Returns
    -------
    results : DataFrame
        Tidy table of enriched pathways.
    """
    # Import gseapy lazily so tests that only use prepare_ranked_list
    # do not require gseapy to be installed.
    try:
        import gseapy as gp
    except ImportError as e:
        raise ImportError(
            "gseapy is required to run GSEA but is not installed. "
            "Install it inside your virtualenv with: 'pip install gseapy'."
        ) from e

    # 1) Prepare ranked list
    rnk = prepare_ranked_list(df, gene_col=gene_col, score_col=score_col)

    if rnk.empty:
        raise ValueError(
            "Ranked list is empty after filtering. "
            "Check that the score column is numeric and not all NaN."
        )

    # 2) Run GSEA prerank
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=config.gene_sets,
        outdir=None,  # don't write to disk, we keep everything in memory
        min_size=config.min_size,
        max_size=config.max_size,
        permutation_num=config.permutation_num,
        seed=config.seed,
        verbose=False,
    )

    # 3) Convert to tidy DataFrame
    results = tidy_gsea_results(pre_res)
    return results

