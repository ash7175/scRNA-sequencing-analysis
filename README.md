# scRNA-seq Analysis: PBMC Cell Type Annotation

Seurat (R) pipeline for single-cell RNA-seq analysis of PBMC data — quality control through cell type annotation. This is the R/Seurat counterpart to the Scanpy/Python pipeline in [SLE-Foundation-models](https://github.com/ash7175/SLE-Foundation-models), covering the same preprocessing stages in a different ecosystem.

## Data

Source: [GSE266852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266852) (10x Genomics PBMC scRNA-seq). The GEO series metadata references RA and SLE samples, but only 3 healthy control samples (barcodes/features/matrix) were available for download at the time of this analysis. This pipeline therefore covers **3 healthy donor samples only** — it is a preprocessing and cell type annotation pipeline, not a disease-vs-healthy comparison.

Raw data is not committed to this repo (see `.gitignore`); download GSE266852 from GEO and place files under `data/raw/` to reproduce.


## Pipeline

### 1. QC & Preprocessing (`01_QC_preprocessing.ipynb`)

- Loaded 3 samples (HC_01, HC_02, HC_03) into Seurat objects (`min.cells = 3`, `min.features = 200`)
- Calculated per-cell mitochondrial and ribosomal percentage
- Filtered: `200 < nFeature_RNA < 5000`, `nCount_RNA > 500`, `percent.mt < 12`

| Sample | Cells before | Cells after | Retained |
|---|---|---|---|
| HC_01 | 8,865 | 6,238 | 70.4% |
| HC_02 | 7,684 | 6,477 | 84.3% |
| HC_03 | 7,667 | 7,057 | 92.0% |
| **Total** | **24,216** | **19,772** | **81.6%** |

- Normalized (LogNormalize) and found variable features (vst, 2000 features) per sample
- Merged all 3 samples, joined layers, re-scaled on merged HVGs, ran PCA (50 PCs)
- Checked for batch effects by sample in PCA space — samples mixed well, no integration needed

### 2. Clustering & Cell Type Annotation (`02_clustering_celltype_annotation.ipynb`)

- `FindNeighbors`/`FindClusters` on PCs 1-15 (per elbow plot), resolution 0.5 → 16 clusters
- `RunUMAP` on the same PCs
- `FindAllMarkers` (positive markers, min.pct 0.25, logFC threshold 0.25) per cluster
- Annotated by cross-referencing cluster markers against canonical PBMC signatures

| Cell type | Cells |
|---|---|
| CD8 T cells | 7,678 |
| CD4 T cells | 4,671 |
| NK cells | 4,290 |
| B cells | 1,450 |
| CD14 Monocytes | 903 |
| Proliferating T/NK cells | 459 |
| FCGR3A Monocytes | 100 |
| Memory B cells | 69 |
| Platelets | 53 |
| Unknown (likely doublets/artifact) | 99 |
| **Total** | **19,772** |

One cluster (99 cells) showed a marker signature inconsistent with any PBMC lineage, alongside the highest UMI/feature counts of all clusters — consistent with a doublet or technical artifact rather than a real cell type. It was retained and labeled `Unknown` rather than silently dropped.

## Requirements

Run `renv::restore()` to install exact package versions from `renv.lock`.