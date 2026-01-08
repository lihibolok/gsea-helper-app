# GSEA Helper App 

## Overview

This project aims to develop a Python-based web application that simplifies running **Gene Set Enrichment Analysis (GSEA)** on differential expression (DE) results (e.g., from bulk RNA-seq or proteomics).  
The tool provides an easy-to-use interface where the user uploads a results table, selects basic parameters, and receives a ranked list of enriched pathways together with publication-ready visualizations.

The goal is to turn the manual process of formatting input files, calling GSEA functions, and plotting the results into a single, streamlined workflow that can be reused in future research projects.

---

## Background: What is GSEA and why is it used?

**Gene Set Enrichment Analysis (GSEA)** is a method for testing whether predefined sets of genes (e.g., pathways, Hallmark signatures, KEGG pathways) are **systematically enriched** at the top or bottom of a ranked list of genes.

Typical use case:

1. You run a differential expression (DE) analysis between two conditions (e.g. KO vs WT, treated vs control) and obtain a table with one row per gene and statistics such as log2 fold change and p-value.
2. Instead of looking at single genes in isolation, GSEA asks whether **groups of functionally related genes** tend to be collectively upregulated or downregulated.
3. This pathway-level view is often more stable and biologically interpretable than focusing only on individual genes, especially in noisy high-dimensional data (bulk RNA-seq, proteomics, etc.).

In this project, I use **pre-ranked GSEA**: genes are ranked by a score (e.g. log2 fold change or a test statistic), and enrichment is tested for each gene set against this ranked list.

### How the app simplifies the standard workflow

A typical GSEA workflow for a student or researcher might look like:

1. Run DESeq2 / edgeR / another DE tool to get a DE results table.
2. Manually clean and reformat the table into a ranked list of genes (choose gene column, choose ranking metric, uppercase gene symbols, remove duplicates, etc.).
3. Download or select a gene set collection (Hallmark, KEGG, etc.) and connect it to the GSEA tool.
4. Call a GSEA function or command-line tool (e.g. gseapy, the original GSEA Java tool), adjusting parameters.
5. Manually inspect the output files and write custom R/Python code to create barplots/dotplots and export them.

The **GSEA Helper App** automates and unifies these steps in one place:

- The user uploads a DE results table directly (CSV/TSV).
- The app takes care of ranking genes (preparing the list for GSEA), choosing the gene set library, and running GSEA under the hood (using `gseapy`).
- Results are returned as a tidy table and visualized immediately (barplot/dotplot), with buttons to **download** the table and plots.
- The whole workflow can be repeated on new datasets just by uploading a different DE file and changing a few options in the UI.

---

## Project Goals

1. **User-friendly GSEA interface**  
   Provide a simple graphical interface (built with Streamlit) for running GSEA without writing code.

2. **Standardized analysis pipeline**  
   Implement a reproducible GSEA pipeline in Python using widely used libraries (e.g., gseapy, pandas).

3. **Clear visual output**  
   Generate clear tables and plots (dotplot / barplot) that summarize which pathways are enriched, their direction (up/down), and statistical significance.

4. **Reusable tool for future projects**  
   Make the repository easy to clone and re-run on new datasets, so it can be used beyond this course.

---

## Functional Specifications

### 1. Input Data

The app accepts a **differential expression (DE) results table**, typically exported from bulk RNA-seq or proteomics pipelines.

**Expected input format:**

- File type: `.csv` or `.tsv`
- Required columns (names are chosen by the user inside the app):
  - A **gene identifier column** (e.g. gene symbol such as *Gnai3* or *CD24a*).  
    The app expects gene symbols for the gene set libraries used (MSigDB/Enrichr-style), and it automatically converts them to uppercase.
  - A **score column** used for ranking, e.g. `log2FoldChange` or a test statistic (e.g. `stat`).
- Optional columns (used for display or external filtering):
  - `pvalue`, `padj` (adjusted p-value / FDR)
  - any other metadata the user wants to keep

The user selects which column should be treated as the gene identifier and which as the ranking score. The app then constructs the ranked gene list internally.

### 2. GSEA Configuration

In the Streamlit interface, the user can:

- Upload a DE results file.
- Select:
  - The **gene column**.
  - The **score column** used to rank genes.
  - The **organism** (initially: *Homo sapiens* and *Mus musculus*).
  - The **gene set collection**, e.g.:
    - Hallmark gene sets (MSigDB Hallmark collection)
    - KEGG pathways

Gene sets are provided via predefined gene set library names supported by `gseapy` (Enrichr/MSigDB collections).  
The app internally constructs a ranked gene list and passes it to the GSEA (pre-rank) function.

### 3. GSEA Analysis

The core analysis is implemented in Python using `gseapy`:

- Run **pre-ranked GSEA** on the ranked gene list against the selected gene set collection.
- Extract, for each pathway:
  - NES (Normalized Enrichment Score)
  - FDR / adjusted p-value
  - Gene set size
  - Leading‐edge or contributing genes
- Add derived columns:
  - `direction` = *upregulated* if NES > 0, *downregulated* if NES < 0.
  - `gene_ratio` = (number of leading-edge genes) / (gene set size).
  - `num_genes` = number of leading-edge genes.
- Return the results as a tidy `pandas.DataFrame` ready for plotting and export.

If no pathways pass the chosen threshold or no gene sets survive the size filters, the app shows a clear message instead of failing.

### 4. Output and Visualization

The app presents the results in two main ways:

1. **Interactive table**
   - Displays the enriched pathways (sorted by NES or adjusted p-value).
   - Shows columns such as: pathway name, NES, FDR (`fdr`), `gene_ratio`, `direction`, and `num_genes`.
   - Allows the user to download the full (filtered) results as a `.csv` file.

2. **Plots**
   - A **barplot** summarizing enriched terms, e.g.:
     - x-axis: NES
     - y-axis: pathway names
     - color: direction (up/down).
   - A **dotplot**, for example:
     - x-axis: NES
     - y-axis: pathway names
     - point size: number of genes (`num_genes`)
     - color: FDR (`fdr`)
   - Plots are generated with `matplotlib` / `seaborn`, shown directly in the app, and can be downloaded as PNG images.

---

## Example Data Used in This Project

The `examples/` directory currently contains two DE result tables:

1. **Toy example (`example_de_results.csv`)**  
   A small, manually constructed DE table with a handful of human gene symbols and made-up log2 fold changes and p-values.  
   This file is included for quick testing of the app and for demonstration purposes (it is not based on real biological data).

2. **Real bulk RNA-seq example (mouse heart, Tv1 KO vs DKO)**  
   A DE results table derived from a **bulk RNA-seq experiment of mouse heart tissue**.  
   The experiment compares two genotypes:
   - **Tv1 KO** – transcription variant 1 of *Cd24* knocked out.
   - **DKO** – double knockout where both transcript variants of *Cd24* are knocked out.  

   Reads were processed and differential expression was calculated using a standard bulk RNA-seq DE pipeline outside of this project.  
   The exported DE table (one row per gene, with gene symbols, log2 fold changes, p-values and adjusted p-values) is used here as input to the GSEA Helper App to demonstrate how the tool can be applied to real biological data.

In both cases, the app does **not** re-do the differential expression analysis; it assumes the DE table is already computed and focuses on the downstream **GSEA** step and visualization.

---

## Quick start with the example data

To test the app quickly with the toy example:

1. Activate the virtualenv and run:

   ```bash
   streamlit run app.py
   ```

2. In the browser, upload:

   `examples/example_de_results.csv`

3. In step 2 of the app, choose:
   - **Gene column:** `gene`
   - **Score column:** `log2FoldChange`

4. In the left sidebar, use these settings:
   - **Organism:** `Homo sapiens`
   - **Gene set collection:** `HALLMARK`
   - **Minimum gene set size:** `5` (a small value is required for the toy example)
   - **Maximum gene set size:** `500` (default is fine)
   - **Number of permutations:** `100` (for a quick run)
   - **FDR q-value cutoff:** e.g. `0.25` or `0.1`

Because this toy table contains only a small number of genes, the minimum gene set size must be low.  
If `min_size` is too high, `gseapy` may return the error:

> GSEA failed: No gene sets passed through filtering condition !!!  
> Hint 1: Try to lower min_size or increase max_size !  
> Hint 2: Check gene symbols are identifiable to your gmt input.  
> Hint 3: Gene symbols curated in Enrichr web services are all upcases.

For the real bulk RNA-seq example (`examples/RNAsc_de_results.csv`, mouse heart Tv1 KO vs DKO), you can use more typical parameters, such as:
- **Organism:** `Mus musculus`
- **Minimum gene set size:** `10–15`
- **Maximum gene set size:** `500–1000`
- **FDR q-value cutoff:** `0.05–0.1`

---

## Technical Details

- **Language:** Python (3.10+)
- **Main libraries (planned and used):**
  - `pandas` (data loading and manipulation)
  - `numpy` (numerical operations)
  - `gseapy` (running GSEA and accessing gene sets)
  - `matplotlib` and `seaborn` (plotting)
  - `streamlit` (web interface / GUI)
  - `pytest` (basic tests)

A `requirements.txt` file lists all dependencies.

---

## Repository Structure (Planned / Implemented)

```text
gsea-helper-app/
├─ app.py                 # Streamlit entry point (UI)
├─ gsea_core.py           # Core GSEA logic: preparing ranked list, running analysis
├─ plotting.py            # Plot creation functions (barplot, dotplot)
├─ requirements.txt
├─ examples/
│  ├─ example_de_results.csv           # Toy DE table
│  └─ [real bulk RNA-seq DE table]     # Mouse heart Tv1 KO vs DKO
├─ tests/
│  ├─ __init__.py
│  └─ test_gsea_core.py    # Basic test for ranked list preparation
└─ README.md               # This document
```

---

## How to Install and Run (Planned)

These instructions describe the expected usage once the project is implemented.

1. **Clone the repository**

   ```bash
   git clone https://github.com/lihibolok/gsea-helper-app.git
   cd gsea-helper-app
   ```

2. **Create and activate a virtual environment (optional but recommended)**

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install dependencies**

   ```bash
   pip install -r requirements.txt
   ```

4. **Start the app**

   ```bash
   streamlit run app.py
   ```

Then open the local URL that Streamlit prints in your browser.  
From there you can upload your DE results file, configure the analysis, and view/download the GSEA results and plots.

---

## Notes

- The exact features and interface may evolve while implementing the project, but the core idea will remain:  
  a simple, reusable tool to run GSEA on differential expression results with minimal manual work.
-   After completion, this README will be updated with final usage instructions.
