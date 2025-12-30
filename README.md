# GSEA Helper App – Project Proposal

## Overview

This project aims to develop a Python-based web application that simplifies running Gene Set Enrichment Analysis (GSEA) on differential expression results (e.g., from bulk RNA-seq or proteomics).  
The tool will provide an easy-to-use interface where the user uploads a results table, selects basic parameters, and receives a ranked list of enriched pathways together with publication-ready visualizations.

The goal is to turn the manual process of formatting input files, calling GSEA functions, and plotting the results into a single, streamlined workflow that can be reused in future research projects.

---

## Project Goals

1. **User-friendly GSEA interface**  
   Provide a simple graphical interface (built with Streamlit) for running GSEA without writing code.

2. **Standardized analysis pipeline**  
   Implement a reproducible GSEA pipeline in Python using widely used libraries (e.g., `gseapy`, `pandas`).

3. **Clear visual output**  
   Generate clear tables and plots (dotplot / barplot) that summarize which pathways are enriched, their direction (up/down), and statistical significance.

4. **Reusable tool for future projects**  
   Make the repository easy to clone and re-run on new datasets, so it can be used beyond this course.

---

## Functional Specifications

### 1. Input Data

The app will accept a differential expression (DE) results table, typically exported from RNA-seq or proteomics pipelines.

**Expected input format:**

- File type: `.csv` or `.tsv`
- Required columns (names can be chosen by the user in the app):
  - A **gene identifier** column (e.g. gene symbol such as `Cd24a`).
  - A **score column** used for ranking, e.g.  
    - `log2FoldChange`, or  
    - a test statistic (e.g. `stat`).
- Optional columns:
  - `pvalue`, `padj` (for filtering or display only).

An example input file will be included in an `examples/` directory.

### 2. GSEA Configuration

In the user interface, the user will be able to:

- Upload a DE results file.
- Select:
  - The **gene column**.
  - The **score column** used to rank genes.
  - The **organism** (initially: `Homo sapiens` and `Mus musculus`).
  - The **gene set collection**, e.g.:
    - HALLMARK gene sets.
    - KEGG pathways.  

Gene sets will be provided via pre-downloaded GMT files or via `gseapy` utilities.

The app will internally construct a ranked gene list (mapping each gene to its score) and pass it to the GSEA function.

### 3. GSEA Analysis

The core analysis will be implemented in Python using `gseapy` (or a similar library):

- Run GSEA on the ranked gene list against the selected gene set collection.
- Extract, for each pathway:
  - Normalized Enrichment Score (`NES`)
  - Adjusted p-value / FDR (`p.adjust`)
  - Gene set size
  - Leading-edge or contributing genes
- Add derived columns:
  - `direction` = `"upregulated"` if NES > 0, `"downregulated"` if NES < 0.
  - `gene_ratio` = (number of leading-edge genes) / (gene set size).
  - `num_genes` = number of leading-edge genes.

Results will be returned as a tidy `pandas.DataFrame` ready for plotting and export.  
If no pathways pass the chosen threshold, the app will display a clear message instead of failing.

### 4. Output and Visualization

The app will present the results in two main ways:

1. **Interactive table**
   - Display the enriched pathways (sorted by NES or adjusted p-value).
   - Show columns such as: pathway name, NES, p.adjust, gene_ratio, direction, and number of genes.
   - Allow the user to download the full results as a `.csv` file.

2. **Plots**
   - A **dotplot** or **barplot** summarizing enriched terms. For example:
     - x-axis: enrichment metric (e.g. `gene_ratio` or `NES`)
     - y-axis: pathway names
     - Point size: number of genes
     - Color: direction and/or adjusted p-value
   - Plots will be generated with `matplotlib` / `seaborn` and shown directly in the app.  
     Optionally, users will be able to download the plot as a PNG image.

---

## Technical Details

### Language and Libraries

- **Language:** Python (3.10+)
- **Main libraries (planned):**
  - `pandas` – data loading and manipulation
  - `numpy` – numerical operations
  - `gseapy` – running GSEA and accessing gene sets
  - `matplotlib` and/or `seaborn` – plotting
  - `streamlit` – web interface / GUI
  - `pytest` (or built-in `unittest`) – basic tests for the analysis functions

A `requirements.txt` file will list all dependencies.

---

## Repository Structure (Planned)

~~~text
gsea-helper-app/
├─ app.py                 # Streamlit entry point
├─ gsea_core.py           # Core GSEA logic (loading data, running analysis)
├─ plotting.py            # Plot creation functions
├─ requirements.txt
├─ examples/
│  └─ example_de_results.csv
├─ tests/
│  └─ test_gsea_core.py
└─ README.md              # This proposal and later full documentation
~~~

---

## How to Install and Run (Planned)

These instructions describe the expected usage once the project is implemented.

1. **Clone the repository**

   ~~~bash
   git clone https://github.com/lihibolok/gsea-helper-app.git
   cd gsea-helper-app
   ~~~

2. **Create and activate a virtual environment (optional but recommended)**

   ~~~bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ~~~

3. **Install dependencies**

   ~~~bash
   pip install -r requirements.txt
   ~~~

4. **Start the app**

   ~~~bash
   streamlit run app.py
   ~~~

Then open the local URL that Streamlit prints in your browser.  
From there you can upload your DE results file, configure the analysis, and view/download the GSEA results and plots.

---

## Notes

- The exact features and interface may evolve while implementing the project, but the core idea will remain:  
  a simple, reusable tool to run GSEA on differential expression results with minimal manual work.
-   After completion, this README will be updated with final usage instructions.
