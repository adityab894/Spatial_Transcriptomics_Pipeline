# Spatial_Transcriptomics_Pipeline

An automated Python pipeline for analyzing spatial gene expression in cancer tissue samples. This tool processes spatial transcriptomics data to deliver key biological insights, including unsupervised tissue clustering, cancer marker analysis, and detailed characterization of the tumor microenvironment (TME).


## ✨ Features

* **Automated Data Processing:** Handles loading, alignment, and quality control of spatial transcriptomics data.
* **Unsupervised Clustering:** Uses the Leiden algorithm to identify distinct biological regions within the tissue based on gene expression profiles.
* **Differential Expression:** Identifies the unique gene signatures that define each cluster.
* **Cancer Marker Analysis:** Pinpoints "hotspots" of cancer activity by analyzing the expression of significant, user-defined cancer genes.
* **Tumor Microenvironment (TME) Analysis:** Scores each spot for the enrichment of various immune and stromal cell types (T-Cells, Macrophages, Fibroblasts, etc.) using a comprehensive database from PanglaoDB.
* **Rich Visualization:** Generates a suite of publication-quality plots, including spatial heatmaps of gene expression, TME cell locations, and clustering results.
* **Comprehensive Reporting:** Outputs all quantitative results to CSV files and generates a summary report for each sample.


## 📂 Project Structure

```
.
├── data/
│   └── sample1/
│       ├── sample1_stdata.csv
│       ├── sample1_coordinates.csv
│       └── sample1_tissue_image.jpg
├── results/
├── .gitignore
├── convert_10x_visium.py
├── convert_data_format.py
├── cancer_markers_comprehensive.json
├── config.json
├── PanglaoDB_markers_27_Mar_2020.tsv
├── README.md
├── requirements.txt
├── run_analysis.py
└── Spatial_Transcriptomics_Pipeline.py
```


## 🛠️ Installation and Setup

### Prerequisites
* Python 3.9+
* The required packages listed in `requirements.txt`.
* The PanglaoDB marker file (`PanglaoDB_markers_27_Mar_2020.tsv`).

  

### Steps
1.  **Clone the repository:**
    ```bash
    git clone (https://github.com/adityab894/Spatial_Transcriptomics_Pipeline.git)
    cd Spatial_Transcriptomics_Pipeline
    ```

2.  **Create a virtual environment (recommended):**
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

3.  **Install the required packages:**
    ```bash
    pip install -r requirements.txt
    ```

    

## 🚀 How to Run the Pipeline

1.  **Prepare Your Data:**
    * Place your sample folders inside the `data/` directory.
    * Each sample folder must contain three files: the expression data (`_stdata.csv`), the coordinate mapping file (`_coordinates.csv`), and the tissue image (`_tissue_image.jpg`).
    * convert_10x_visium.py and convert_data_format.py were used to prepare the data in the required data formats before running the pipeline.

2.  **Configure the Pipeline:**
    * Copy the configuration template: `cp config.json`.
    * Edit `config.json` to adjust parameters like file paths, QC thresholds, or clustering resolution as needed.

3.  **Run the Analysis:**
    * **To run the entire batch** on all samples in the `data/` directory (recommended):
        ```bash
        python run_analysis.py
        ```
    * **To run on a single sample:**
        ```bash
        python Spatial_pipeline.py --sample ./data/sample_name
        ```


## 🔬 Dataset

The sample data required to run this pipeline is publicly available from the NCBI Gene Expression Omnibus (GEO) under accession number **GSE144239**. One can use their own data to run the pipeline.

You can download the dataset here: [(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi)] and provide the accession number.

**Preparation:**
After downloading, please use the provided helper scripts (`convert_10x_visium.py`, `convert_data_format.py`) as needed to format the data into the structure expected by the pipeline (described under "Input Data Format").


## 📊 Input Data Format

For each sample, the pipeline expects the following files within its folder in `data/`:

* **`*_stdata.csv`:** The gene expression matrix. Rows should be spots (e.g., `10x12`) and columns should be gene symbols.
* **`*_coordinates.csv`:** A file mapping spot IDs to their pixel locations on the full-resolution tissue image. Must contain columns `spot`, `pixel_x`, and `pixel_y`.
* **`*_tissue_image.jpg` (or other format):** The high-resolution histology image for the tissue sample.

## 📈 Output Explanation

The pipeline will create a new folder for each sample inside the `results/` directory, containing:

* **`*_clustering.png`:** Plots showing the results of unsupervised clustering (PCA, UMAP, and spatial distribution).
* **`*_de_results.csv`:** A table of differentially expressed genes for each identified cluster.
* **`*_spatial_expression.png`:** A spatial plot showing the expression of significant cancer markers on the tissue image.
* **`*_tme_scores.csv`:** The calculated enrichment scores for each TME cell type at every spot.
* **`*_tme_spatial_scores.png`:** A multi-panel plot visualizing the location of TME cell type hotspots.
* **`*_analysis_report.md`:** A summary report of the analysis run.

## 📄 License

This project is licensed under the MIT License. See the `LICENSE` file for details.
