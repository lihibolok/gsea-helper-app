"""
Plotting utilities for GSEA Helper App.

Uses matplotlib + seaborn to build:
- barplot of top enriched pathways
- dotplot of top enriched pathways
"""

from __future__ import annotations

from typing import Literal

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def barplot_enrichment(
    results: pd.DataFrame,
    top_n: int = 20,
    sort_by: Literal["NES", "fdr"] = "NES",
    ascending: bool = False,
) -> plt.Figure:
    """
    Create a barplot of top enriched pathways.

    Parameters
    ----------
    results : DataFrame
        Tidy GSEA results (from gsea_core.run_gsea).
    top_n : int
        Number of pathways to show.
    sort_by : {"NES", "fdr"}
        Metric used to select top pathways.
    ascending : bool
        Sort order (for NES usually False, for fdr usually True).

    Returns
    -------
    fig : matplotlib Figure
    """
    if results.empty:
        raise ValueError("Results DataFrame is empty. Nothing to plot.")

    if sort_by not in results.columns:
        raise ValueError(f"Column {sort_by!r} not found in results.")

    # Select top_n by chosen metric
    df = results.sort_values(sort_by, ascending=ascending).head(top_n).copy()

    # For plotting, we usually want small values at bottom and large at top
    df = df.sort_values(sort_by, ascending=True)

    # Auto-adjust figure height based on number of pathways
    fig_height = max(4, 0.4 * len(df))
    fig, ax = plt.subplots(figsize=(8, fig_height))

    sns.barplot(
        data=df,
        x=sort_by,
        y="pathway",
        hue="direction" if "direction" in df.columns else None,
        ax=ax,
    )

    ax.set_xlabel(sort_by)
    ax.set_ylabel("Pathway")
    if "direction" in df.columns:
        ax.legend(title="Direction", bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        ax.legend_.remove()

    fig.tight_layout()
    return fig


def dotplot_enrichment(
    results: pd.DataFrame,
    top_n: int = 20,
    sort_by: Literal["NES", "fdr"] = "NES",
    ascending: bool = False,
    size_col: str = "num_genes",
    color_col: str = "fdr",
) -> plt.Figure:
    """
    Create a dotplot of top enriched pathways.

    Parameters
    ----------
    results : DataFrame
        Tidy GSEA results (from gsea_core.run_gsea).
    top_n : int
        Number of pathways to show.
    sort_by : {"NES", "fdr"}
        Metric used to select top pathways.
    ascending : bool
        Sort order for selecting top_n.
    size_col : str
        Column controlling dot size.
    color_col : str
        Column controlling dot color.

    Returns
    -------
    fig : matplotlib Figure
    """
    if results.empty:
        raise ValueError("Results DataFrame is empty. Nothing to plot.")

    for col in [sort_by, size_col, color_col]:
        if col not in results.columns:
            raise ValueError(f"Column {col!r} not found in results.")

    df = results.sort_values(sort_by, ascending=ascending).head(top_n).copy()
    df = df.sort_values(sort_by, ascending=True)

    fig_height = max(4, 0.4 * len(df))
    fig, ax = plt.subplots(figsize=(8, fig_height))

    # Normalize size a bit for plotting: scale to something reasonable
    sizes = df[size_col].fillna(1).astype(float)
    size_scale = 100 / sizes.max() if sizes.max() > 0 else 1.0

    scatter = ax.scatter(
        x=df[sort_by],
        y=df["pathway"],
        s=sizes * size_scale * 50,  # size factor
        c=df[color_col],
        cmap="viridis",
        alpha=0.8,
    )

    ax.set_xlabel(sort_by)
    ax.set_ylabel("Pathway")

    # Colorbar for fdr or other color_col
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(color_col)

    fig.tight_layout()
    return fig
