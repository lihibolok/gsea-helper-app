"""
Streamlit interface for the GSEA Helper App.

Steps:
1. Upload DE results (CSV/TSV)
2. Choose gene and score columns
3. Select organism and gene set collection
4. Run GSEA
5. View/download results and plots
"""

from __future__ import annotations

import io

import pandas as pd
import streamlit as st

from gsea_core import GSEAConfig, run_gsea
from plotting import barplot_enrichment, dotplot_enrichment


def _detect_sep(filename: str) -> str:
    """Very simple separator heuristic based on file extension."""
    if filename.lower().endswith(".tsv") or filename.lower().endswith(".txt"):
        return "\t"
    return ","


def main() -> None:
    st.set_page_config(
        page_title="GSEA Helper App",
        layout="wide",
    )

    st.title("GSEA Helper App")
    st.markdown(
        """
This app helps you run **Gene Set Enrichment Analysis (GSEA)** on differential
expression (DE) results without writing code.

**Workflow:**

1. Upload a CSV or TSV file with your DE results.
2. Select the gene column and the score column (e.g. log2FoldChange or stat).
3. Choose the organism and gene set collection.
4. Click **Run GSEA** to see enriched pathways, tables, and plots.
        """
    )

    # --- Sidebar configuration ---
    st.sidebar.header("GSEA Configuration")

    organism = st.sidebar.selectbox(
        "Organism",
        options=["Homo sapiens", "Mus musculus"],
        index=0,
    )

    gene_set_collection = st.sidebar.selectbox(
        "Gene set collection",
        options=["HALLMARK", "KEGG"],
        index=0,
    )

    min_size = st.sidebar.number_input(
        "Minimum gene set size",
        min_value=5,
        max_value=1000,
        value=15,
        step=1,
    )

    max_size = st.sidebar.number_input(
        "Maximum gene set size",
        min_value=min_size,
        max_value=5000,
        value=500,
        step=10,
    )

    permutation_num = st.sidebar.number_input(
        "Number of permutations",
        min_value=10,
        max_value=2000,
        value=100,
        step=10,
        help="Higher values give more stable p-values but take longer to run.",
    )

    fdr_cutoff = st.sidebar.slider(
        "FDR q-value cutoff",
        min_value=0.0,
        max_value=0.25,
        value=0.05,
        step=0.01,
    )

    top_n = st.sidebar.slider(
        "Number of top pathways in plots",
        min_value=5,
        max_value=50,
        value=20,
        step=1,
    )

    # --- File upload ---
    st.subheader("1. Upload differential expression (DE) results")

    uploaded_file = st.file_uploader(
        "Upload a CSV or TSV file with at least a gene column and a score column.",
        type=["csv", "tsv", "txt"],
    )

    if uploaded_file is None:
        st.info("Please upload a DE results file to continue.")
        return

    # Read the file
    sep = _detect_sep(uploaded_file.name)
    try:
        df = pd.read_csv(uploaded_file, sep=sep)
    except Exception as e:
        st.error(f"Could not read file: {e}")
        return

    if df.empty:
        st.error("The uploaded file is empty.")
        return

    st.markdown("**Preview of uploaded data:**")
    st.dataframe(df.head())

    # --- Column selection ---
    st.subheader("2. Select gene and score columns")

    all_columns = list(df.columns)

    if not all_columns:
        st.error("No columns found in the uploaded file.")
        return

    gene_col = st.selectbox(
        "Gene identifier column (e.g. gene symbol)",
        options=all_columns,
        index=0,
    )

    # Try to guess a score column
    default_score_idx = 0
    for candidate in ["log2FoldChange", "log2FC", "stat", "score"]:
        if candidate in all_columns:
            default_score_idx = all_columns.index(candidate)
            break

    score_col = st.selectbox(
        "Score column used to rank genes (e.g. log2FoldChange or stat)",
        options=all_columns,
        index=default_score_idx,
    )

    st.markdown(
        """
The score column should be numeric and reflect the **direction and strength**
of differential expression (e.g. a test statistic or log2 fold change).
        """
    )

    # --- Run GSEA button ---
    st.subheader("3. Run GSEA")

    # Initialize session state for GSEA results
    if "gsea_results" not in st.session_state:
        st.session_state["gsea_results"] = None

    run_clicked = st.button("Run GSEA", type="primary")

    # If the button was clicked, run GSEA and store results
    if run_clicked:
        with st.spinner("Running GSEA... this may take a moment."):
            try:
                config = GSEAConfig(
                    organism=organism,
                    gene_set_collection=gene_set_collection,
                    min_size=int(min_size),
                    max_size=int(max_size),
                    permutation_num=int(permutation_num),
                    seed=42,
                )
                results = run_gsea(
                    df=df,
                    gene_col=gene_col,
                    score_col=score_col,
                    config=config,
                )
            except Exception as e:
                st.error(f"GSEA failed: {e}")
                st.session_state["gsea_results"] = None
            else:
                if results.empty:
                    st.warning(
                        "GSEA completed, but no pathways were returned. "
                        "Try relaxing the parameters or checking the input file."
                    )
                    st.session_state["gsea_results"] = None
                else:
                    # Store full results; filtering by FDR is done below
                    st.session_state["gsea_results"] = results

    # --- Show results if we have them in session_state ---
    results = st.session_state.get("gsea_results")
    if results is None:
        # Nothing to show yet (either not run, or last run failed)
        return

    # Filter by FDR using the current slider value
    if "fdr" in results.columns:
        sig_results = results[results["fdr"] <= fdr_cutoff].copy()
    else:
        sig_results = results.copy()

    if sig_results.empty:
        st.warning(
            f"GSEA completed, but no pathways passed the FDR cutoff ({fdr_cutoff}). "
            "Showing full results below."
        )
        sig_results = results.copy()

    # --- Show table ---
    st.subheader("4. Enriched pathways (table)")

    st.markdown(
        f"Showing pathways with FDR ≤ {fdr_cutoff} "
        f"(or all pathways if none passed the threshold)."
    )
    st.dataframe(sig_results)

    # Download as CSV
    csv_buffer = io.StringIO()
    sig_results.to_csv(csv_buffer, index=False)
    st.download_button(
        label="Download table as CSV",
        data=csv_buffer.getvalue(),
        file_name="gsea_results_filtered.csv",
        mime="text/csv",
    )

    # --- Plots ---
    st.subheader("5. Plots")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Barplot**")
        try:
            fig_bar = barplot_enrichment(
                sig_results,
                top_n=top_n,
                sort_by="NES",
                ascending=False,
            )
            st.pyplot(fig_bar)

            # Save barplot to PNG in memory
            bar_buf = io.BytesIO()
            fig_bar.savefig(bar_buf, format="png", bbox_inches="tight", dpi=300)
            bar_buf.seek(0)

            st.download_button(
                label="Download barplot as PNG",
                data=bar_buf,
                file_name="gsea_barplot.png",
                mime="image/png",
            )
        except Exception as e:
            st.error(f"Could not create barplot: {e}")

    with col2:
        st.markdown("**Dotplot**")
        try:
            fig_dot = dotplot_enrichment(
                sig_results,
                top_n=top_n,
                sort_by="NES",
                ascending=False,
                size_col="num_genes",
                color_col="fdr",
            )
            st.pyplot(fig_dot)

            # Save dotplot to PNG in memory
            dot_buf = io.BytesIO()
            fig_dot.savefig(dot_buf, format="png", bbox_inches="tight", dpi=300)
            dot_buf.seek(0)

            st.download_button(
                label="Download dotplot as PNG",
                data=dot_buf,
                file_name="gsea_dotplot.png",
                mime="image/png",
            )
        except Exception as e:
            st.error(f"Could not create dotplot: {e}")

    st.success("GSEA completed successfully.")


if __name__ == "__main__":
    main()
