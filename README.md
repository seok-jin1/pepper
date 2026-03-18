# ISS Chile Pepper Rhizosphere Microbiome

**Can rhizosphere microbiomes support plants in space? Microbiomic evidence for functional network collapse in ISS-grown *Capsicum annuum***

Seok-Jin Kang, Hongchul Shin
*Plant Physiology Letters* (2026, revised submission)
Department of Biotechnology, Korea University, Seoul, Republic of Korea

---

## Overview

This repository contains the complete bioinformatics pipeline used in the study. We perform a multi-level comparative analysis of the *Capsicum annuum* rhizosphere microbiome grown aboard the International Space Station (ISS) versus terrestrial soil controls, integrating:

- **16S rRNA amplicon sequencing** (QIIME2 + DADA2, V3–V4)
- **FAPROTAX functional trait annotation** (curated, v1.2.12; Louca et al. 2016)
- **Predicted functional pathways** (PICRUSt2 / MetaCyc)
- **Compositional differential abundance** (ANCOM-BC; Lin & Peddada 2020)
- **Community assembly null models** (βNTI/RCbray; picante)
- **Genus-level co-occurrence network analysis** (Spearman |r| > 0.4, p < 0.05)
- **Network robustness & threshold sensitivity simulations**
- **Functional redundancy and N-cycle completeness quantification**

**Key finding:** Spaceflight restructures the rhizosphere toward a stochastically assembled, generalist-dominated community in which nitrogen-cycling capacity has structurally collapsed — nitrification abundance is >1000× lower in Space (FAPROTAX Log₂FC = −10.68, q < 0.001), nitrification functional redundancy drops from 4.80 to 0.20 contributing genera, and the keystone hub *Rhodanobacter* is absent from the Space co-occurrence network at all tested correlation thresholds.

---

## Datasets

| Dataset | Description | n (analysis) | Access |
|---------|-------------|--------------|--------|
| **OSD-772** | ISS-grown *C. annuum* rhizosphere (NASA VEGGIE) | 106 (Space Flight) | [NASA OSDR](https://osdr.nasa.gov/bio/repo/data/studies/OSD-772) |
| **PRJNA1145089** | Terrestrial soil-grown *C. annuum* rhizosphere | 20 (Terrestrial Soil) | [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1145089) |

> **Note on OSD-772 sample heterogeneity:** The 106 Space Flight samples span diverse sample types (rhizosphere, arcillite, wick, stem, fruit, leaf, seed, foam, swab, water). Script 24 (`24_sample_filter_sensitivity.py`) validates that all filter levels (ALL n=106, ROOT_ZONE n=36, STRICT rhizosphere n=4) yield the same direction of difference vs. Terrestrial.

> **Note on sample groups:** The integrated metadata contains three groups — `Space_Flight` (n=106), `Ground_Seed` (ISS experiment ground controls within OSD-772, n=3), and `Terrestrial_Soil` (n=20). The manuscript compares **Space_Flight vs. Terrestrial_Soil** only. `Ground_Seed` is excluded from all main analyses.

> **Raw data not included.** Download instructions are in [Step 0](#step-0-download-terrestrial-reference-data-prjna1145089) and [Step 1](#step-1-prepare-iss-data-osd-772) below.

---

## Repository Structure

```
.
├── README.md
├── CLAUDE.md                         # Claude Code project instructions
├── RESEARCH_NOTES.md                 # Post-rejection analysis plan and notes
├── config.py                         # Centralized paths and analysis parameters
├── .gitignore
│
├── envs/
│   ├── qiime2-amplicon.yml           # Core environment (QIIME2, Python, R scripts)
│   ├── picrust2.yml                  # Functional prediction (PICRUSt2)
│   └── pepper-network.yml            # R-based SpiecEasi network (optional)
│
├── osd772_metadata.tsv               # OSD-772 sample metadata (NASA OSDR, bundled)
│
├── 00_download_data.py               # Download PRJNA1145089 from NCBI ENA
├── 01_create_metadata.py             # Merge OSD-772 and terrestrial metadata
├── 02_generate_manifest.py           # Generate paired-end QIIME2 manifest
├── 03_sync_ids.py                    # Synchronize sample IDs across datasets
├── 04_fix_ids.py                     # Fix sample ID formatting issues
├── 05_qiime2_pipeline.sh             # Full QIIME2 amplicon processing pipeline
├── 06_plot_depth.py                  # Sequencing depth → rarefaction threshold
├── 07_supp_figures.py                # → Main_Fig1B_Taxa_Barplot.png/.pdf
├── 08_picrust2_analysis.py           # → Main_Fig1C_Functional_Pathways.png/.pdf
├── 09_taxon_function_corr.py         # → Main_Fig2D_Taxon_Function_Correlation.png/.pdf
├── 10_network_analysis.py            # → Main_Fig2E_Network.png/.pdf
├── 13_supp_faith_pd.py               # → Supp_S1_FaithPD.png/.pdf  [Supp Fig S1]
├── 14_supp_temporal.py               # → Supp_S3_Temporal_Q1_Q4.png/.pdf  [Supp Fig S3]
│
│   ── Post-rejection mechanistic analyses (Scripts 15–27) ──────────────────
│
├── 15_faprotax_analysis.py           # → Main_Fig2a_FAPROTAX_Log2FC.png/.pdf  [Fig 2a]
├── 16_bnti_rcbray_analysis.R         # → Supp_S4_BNTI_RCbray.png/.pdf  [Supp Fig S4]
├── 17_faprotax_network.py            # FAPROTAX × network hub cross-reference (supporting)
├── 18_ancombc_plot.py                # ANCOM-BC differential abundance (supporting)
├── 19_indval.R                       # IndVal indicator species (supporting)
├── 20_network_robustness.py          # Network robustness simulation (supporting)
├── 21_beta_dispersion.R              # PERMDISP2 beta-dispersion (supporting)
├── 22_functional_redundancy.py       # → Main_Fig2c_FunctionalRedundancy.png/.pdf  [Fig 2c]
├── 23_pgp_index.py                   # PGP bacterial index (supporting)
├── 24_sample_filter_sensitivity.py   # Sample filter sensitivity (supporting)
├── 25_permanova.R                    # → Supp_S2_PERMANOVA.png/.pdf  [Supp Fig S2]
├── 26_ncycle_pgp_specificity.py      # → Main_Fig2b_NCycle_PGP.png/.pdf  [Fig 2b]
├── 27_network_threshold_sensitivity.py # → Supp_S5_NetworkSensitivity.png/.pdf  [Supp Fig S5]
│
└── version-2_integrated/             # All intermediate TSVs, CSVs, and figures
```

---

## Environment Setup

```bash
# Core analysis (QIIME2, Python scripts 06–27, R scripts via Rscript)
conda env create -f envs/qiime2-amplicon.yml
conda activate qiime2-amplicon

# Functional prediction (PICRUSt2 v2.5+)
conda env create -f envs/picrust2.yml

# R-based SpiecEasi network (optional, not used in main figures)
conda env create -f envs/pepper-network.yml
```

All Python scripts (06–27) and R scripts (16, 19, 21, 25) run in the `qiime2-amplicon` environment unless noted. Required R packages: `vegan`, `picante`, `ggplot2`, `dplyr`, `patchwork`, `indicspecies`.

---

## Full Pipeline — End-to-End Reproduction

> **All scripts must be run from the repository root** (the directory containing `config.py`):
> ```bash
> cd /path/to/pepper/
> ```
> Scripts use paths relative to this root (e.g., `version-2_integrated/`). Running from a different directory will cause file-not-found errors.

### Step 0: Download Terrestrial Reference Data (PRJNA1145089)

```bash
conda activate qiime2-amplicon
python 00_download_data.py
```

Downloads 16S amplicon FASTQ files from NCBI ENA to `external_data/PRJNA1145089/`. To obtain paired-end R1/R2 files required for QIIME2:

```bash
while read acc; do
  fasterq-dump --split-files --outdir external_data/PRJNA1145089 "$acc"
  gzip external_data/PRJNA1145089/${acc}_1.fastq external_data/PRJNA1145089/${acc}_2.fastq
done < external_list.txt
```

---

### Step 1: Prepare ISS Data (OSD-772)

Download raw paired-end FASTQ files from [NASA OSDR OSD-772](https://osdr.nasa.gov/bio/repo/data/studies/OSD-772) into `external_data/OSD-772/`:

```
external_data/OSD-772/GLDS-675_GAmplicon_<sample-id>_R1_raw.fastq.gz  (109 files)
external_data/OSD-772/GLDS-675_GAmplicon_<sample-id>_R2_raw.fastq.gz  (109 files)
```

`osd_r1_list.txt` and `osd_r2_list.txt` are pre-configured for this path and included in the repository.

---

### Step 2: Integrate Metadata

```bash
python 01_create_metadata.py
python 02_generate_manifest.py
python 03_sync_ids.py
python 04_fix_ids.py
```

Output: `version-2_integrated/integrated_metadata.tsv`
→ Three groups: `Space_Flight` (n=106), `Ground_Seed` (n=3), `Terrestrial_Soil` (n=20)

---

### Step 3: QIIME2 Amplicon Processing

```bash
conda activate qiime2-amplicon

# Download SILVA 138 V3-V4 reference (one-time):
mkdir -p classifiers
wget "https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza" \
  -O classifiers/silva-138-99-seqs.qza
wget "https://data.qiime2.org/2024.5/common/silva-138-99-tax.qza" \
  -O classifiers/silva-138-99-tax.qza

# Train Naive Bayes classifier (one-time, ~1 hr, ~16 GB RAM):
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  classifiers/silva-138-99-seqs.qza \
  --i-reference-taxonomy classifiers/silva-138-99-tax.qza \
  --o-classifier classifiers/silva-138-99-classifier.qza

# Run full pipeline:
bash 05_qiime2_pipeline.sh
```

Pipeline steps (in order):
1. Import paired-end reads (PairedEndFastqManifestPhred33V2)
2. DADA2 denoising (trunc-len-f 240, trunc-len-r 200)
3. Export raw table → `exported_table_raw/`
4. Taxonomic classification (SILVA 138)
5. Filter mitochondria & chloroplast
6. Phylogenetic tree (MAFFT + FastTree)
7. Core diversity metrics (rarefaction depth: 1,000 reads)
8. Export tables as TSV (combined, Space-only, Terrestrial-only)
9. PICRUSt2 functional prediction (run separately in `picrust2` env)

> **ANCOM-BC prerequisite:** Script 18 requires QIIME2 ANCOM-BC outputs. Run before `18_ancombc_plot.py`:
> ```bash
> # See ancombc_out/ for pre-generated outputs from this study
> ```

---

### Steps 4–9: Original Figure Scripts

```bash
python 06_plot_depth.py          # Supp: sequencing depth distribution
python 07_supp_figures.py        # Fig 1b — phylum barplot
python 08_picrust2_analysis.py   # Fig 1c — MetaCyc pathway Log₂FC
python 09_taxon_function_corr.py # Fig 2e — taxon–function heatmap
python 10_network_analysis.py    # Fig 2f — co-occurrence network
python 13_supp_faith_pd.py       # Supp Fig S1 — Faith's PD
python 14_supp_temporal.py       # Supp Fig S3 — Q1–Q4 temporal
```

---

### Steps 10–18: Post-Rejection Mechanistic Analyses (Scripts 15–27)

Run in numerical order. All output to `version-2_integrated/`.

```bash
# Step 10 — Fig 2a: FAPROTAX functional trait annotation
python 15_faprotax_analysis.py
# Requires FAPROTAX 1.2.12 at /home/laugh/pepper/FAPROTAX_1.2.12/
# Edit FAPROTAX_DIR in script if installed elsewhere

# Step 11 — Fig 2d: βNTI/RCbray community assembly null model (~10–30 min)
# Requires rooted tree exported from QIIME2:
qiime tools export \
  --input-path version-2_integrated/rooted-tree.qza \
  --output-path version-2_integrated/exported_tree
Rscript 16_bnti_rcbray_analysis.R

# Step 12 — Fig 2c: Functional redundancy
python 22_functional_redundancy.py

# Step 13 — Supp Fig S2: PERMANOVA + subsampling validation (1000 iterations, ~5–10 min)
Rscript 25_permanova.R

# Step 14 — Fig 2b: N-cycle completeness + PGP specificity
python 26_ncycle_pgp_specificity.py

# Step 15 — Supp Fig S4: Network threshold sensitivity (|r| = 0.3 / 0.4 / 0.5)
python 27_network_threshold_sensitivity.py

# Step 16 — Supporting: FAPROTAX × network hub cross-reference
python 17_faprotax_network.py

# Step 17 — Supporting: ANCOM-BC differential abundance (genus level)
python 18_ancombc_plot.py
# Reads from ancombc_out/ancombc-genus-exported/

# Step 18 — Supporting: IndVal, network robustness, beta-dispersion, PGP, filter sensitivity
Rscript 19_indval.R
python 20_network_robustness.py
Rscript 21_beta_dispersion.R
python 23_pgp_index.py
python 24_sample_filter_sensitivity.py
```

---

## Output Files

All files in `version-2_integrated/` unless otherwise noted.

### Main Figures

| Figure | Script | Figure output | Data output |
|--------|--------|---------------|-------------|
| **Fig 1b** | `07_supp_figures.py` | `Main_Fig1B_Taxa_Barplot.png/.pdf` | — |
| **Fig 1c** | `08_picrust2_analysis.py` | `Main_Fig1C_Functional_Pathways.png/.pdf` | — |
| **Fig 2a** | `15_faprotax_analysis.py` | `Main_Fig2a_FAPROTAX_Log2FC.png/.pdf` | `FAPROTAX_functional_table.tsv`<br>`FAPROTAX_group_comparison.csv`<br>`FAPROTAX_key_functions.csv`<br>`FAPROTAX_report.txt` |
| **Fig 2b** | `26_ncycle_pgp_specificity.py` | `Main_Fig2b_NCycle_PGP.png/.pdf` | `NCycle_Completeness_Results.csv`<br>`NCycle_Completeness_Summary.csv`<br>`PGP_Specificity_Results.csv` |
| **Fig 2c** | `22_functional_redundancy.py` | `Main_Fig2c_FunctionalRedundancy.png/.pdf` | `Functional_Redundancy_Results.csv`<br>`Functional_Redundancy_Summary.csv` |
| **Fig 2d** | `09_taxon_function_corr.py` | `Main_Fig2D_Taxon_Function_Correlation.png/.pdf` | — |
| **Fig 2e** | `10_network_analysis.py` | `Main_Fig2E_Network.png/.pdf` | `Network_Keystone_Genera.csv` |

### Supplementary Figures

| Supp Fig | Script | Figure output | Data output |
|----------|--------|---------------|-------------|
| **S1** | `13_supp_faith_pd.py` | `Supp_S1_FaithPD.png/.pdf` | — |
| **S2** | `25_permanova.R` | `Supp_S2_PERMANOVA.png/.pdf` | `PERMANOVA_results.csv`<br>`PERMANOVA_subsampling.csv` |
| **S3** | `14_supp_temporal.py` | `Supp_S3_Temporal_Q1_Q4.png/.pdf` | — |
| **S4** | `16_bnti_rcbray_analysis.R` | `Supp_S4_BNTI_RCbray.png/.pdf` | `BNTI_RCbray_results.csv`<br>`BNTI_RCbray_summary.csv` |
| **S5** | `27_network_threshold_sensitivity.py` | `Supp_S5_NetworkSensitivity.png/.pdf` | `NetworkSensitivity_Results.csv` |

---

## Figure–Script Mapping

### Main Figures

| Figure | Panel | Script | Key result |
|--------|-------|--------|------------|
| **Fig 1a** | Overview schematic | — | Conceptual diagram |
| **Fig 1b** | Phylum barplot | `07_supp_figures.py` | Proteobacteria shift |
| **Fig 1c** | PICRUSt2 Log₂FC | `08_picrust2_analysis.py` | Nitrifier denitrification −11.5 |
| **Fig 2a** | FAPROTAX Log₂FC barplot | `15_faprotax_analysis.py` | Nitrification Log₂FC = −10.68 |
| **Fig 2b** | N-cycle step presence heatmap | `26_ncycle_pgp_specificity.py` | Space 73% vs. Terrestrial 100% |
| **Fig 2c** | Functional redundancy barplot | `22_functional_redundancy.py` | Nitrification 0.20 vs. 4.80 genera |
| **Fig 2d** | Taxon–function heatmap | `09_taxon_function_corr.py` | Rhodanobacter–nitrification correlation |
| **Fig 2e** | Co-occurrence network | `10_network_analysis.py` | Rhodanobacter hub collapse |

### Supplementary Figures

| Supp Fig | Script | Description |
|----------|--------|-------------|
| **S1** | `13_supp_faith_pd.py` | Faith's PD: Space 10.05 vs. Terrestrial 24.98 |
| **S2** | `25_permanova.R` | PERMANOVA bootstrap validation (1000 iterations) |
| **S3** | `14_supp_temporal.py` | Q1–Q4 temporal stability (Faith's PD, Shannon) |
| **S4** | `16_bnti_rcbray_analysis.R` | βNTI/RCbray distribution: Terrestrial 93.2% deterministic vs. Space <0.1% |
| **S5** | `27_network_threshold_sensitivity.py` | Network metrics across |r| = 0.3/0.4/0.5 |

---

## Key Results Summary

| Analysis | Space Flight | Terrestrial Soil | Statistic |
|----------|-------------|-----------------|-----------|
| PERMANOVA R² (ROOT_ZONE) | — | — | F=70.15, R²=0.565, p=0.001 |
| Shannon diversity | 2.66 ± — | 6.23 ± — | p < 0.001 |
| Beta-dispersion (Bray-Curtis) | 0.574 ± 0.118 | 0.398 ± 0.082 | F=40.45, p=0.001 |
| Nitrification (FAPROTAX Log₂FC) | — | — | −10.68, q < 0.001 |
| N-cycle completeness | 0.730 ± 0.193 | 1.000 ± 0.000 | p < 0.001 |
| Nitrification redundancy (genera) | 0.20 | 4.80 | q < 0.001 |
| Deterministic assembly (βNTI) | <0.1% | 93.2% | — |
| Rhodanobacter hub degree | 0 (all thresholds) | 6–8 | — |

---

## Reproducibility

- Core analysis parameters and paths for scripts 08–10 are centralized in `config.py` (`RANDOM_SEED = 42`); scripts 00–07, 11–14 use hardcoded relative paths
- Validate all required input files before running:
  ```bash
  python config.py
  ```
- MetaCyc pathway descriptions: `version-2_integrated/metacyc_pathways_info.txt.gz`
- ANCOM-BC pre-generated outputs: `version-2_integrated/ancombc_out/`

---

## Funding

This work was supported by the Google Cloud Research Credits program (award GCP477932969).

---

## Contact

Hongchul Shin (corresponding author)
Department of Biotechnology, College of Life Sciences and Biotechnology, Korea University, Seoul 02841, Republic of Korea
saekomi5@korea.ac.kr
