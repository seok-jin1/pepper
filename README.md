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
├── 07_supp_figures.py                # Fig 1B (taxa barplot) + QC figures
├── 08_picrust2_analysis.py           # Fig 1C — MetaCyc pathway Log₂FC
├── 09_taxon_function_corr.py         # Fig 1D — Taxon–function Spearman heatmap
├── 10_network_analysis.py            # Fig 1E — Genus co-occurrence network + keystone
├── 11_pathogen_screening.py          # Exploratory: plant pathogen screening
├── 12_identify_asvs.py               # Identify dominant Space ASVs
├── 13_supp_faith_pd.py               # Supp Fig S1 — Faith's PD
├── 14_supp_temporal.py               # Supp Fig S3 — Q1–Q4 temporal (old S2)
│
│   ── Post-rejection mechanistic analyses (Scripts 15–28) ──────────────────
│
├── 15_faprotax_analysis.py           # Fig 2a — FAPROTAX Log₂FC barplot
├── 16_bnti_rcbray_analysis.R         # Supp Fig S4 — βNTI/RCbray assembly null model
├── 17_faprotax_network.py            # FAPROTAX × network hub cross-reference
├── 18_ancombc_plot.py                # ANCOM-BC differential abundance plot
├── 19_indval.R                       # IndVal indicator species analysis
├── 20_network_robustness.py          # Network robustness simulation (targeted vs. random)
├── 21_beta_dispersion.R              # PERMDISP2 beta-dispersion
├── 22_functional_redundancy.py       # Functional redundancy per FAPROTAX function
├── 23_pgp_index.py                   # Plant growth-promoting (PGP) bacterial index
├── 24_sample_filter_sensitivity.py   # Sample filter sensitivity (ALL/ROOT_ZONE/STRICT)
├── 25_permanova.R                    # PERMANOVA + 1000-iteration subsampling validation
├── 26_ncycle_pgp_specificity.py      # Fig 2b — N-cycle completeness + PGP specificity
├── 27_network_threshold_sensitivity.py # Network robustness across |r| = 0.3/0.4/0.5
├── 28_temporal_stability.py          # Temporal stability Q1–Q4 (FAPROTAX metrics)
│
└── version-2_integrated/             # All intermediate TSVs, CSVs, and figures
```

---

## Environment Setup

```bash
# Core analysis (QIIME2, Python scripts 06–28, R scripts via Rscript)
conda env create -f envs/qiime2-amplicon.yml
conda activate qiime2-amplicon

# Functional prediction (PICRUSt2 v2.5+)
conda env create -f envs/picrust2.yml

# R-based SpiecEasi network (optional, not used in main figures)
conda env create -f envs/pepper-network.yml
```

All Python scripts (06–28) and R scripts (16, 19, 21, 25) run in the `qiime2-amplicon` environment unless noted. Required R packages: `vegan`, `picante`, `ggplot2`, `dplyr`, `patchwork`, `indicspecies`.

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
python 07_supp_figures.py        # Fig 1B (taxa barplot) + QC figures
python 08_picrust2_analysis.py   # Fig 1C — MetaCyc pathway Log₂FC
python 09_taxon_function_corr.py # Fig 1D — taxon–function heatmap
python 10_network_analysis.py    # Fig 1E — co-occurrence network
python 13_supp_faith_pd.py       # Supp Fig S1 — Faith's PD
python 14_supp_temporal.py       # Supp Fig S3 — Q1–Q4 temporal
```

---

### Steps 10–23: Post-Rejection Mechanistic Analyses (Scripts 15–28)

Run in numerical order. All output to `version-2_integrated/`.

```bash
# Step 10 — FAPROTAX functional trait annotation
python 15_faprotax_analysis.py
# Requires FAPROTAX 1.2.12 at /home/laugh/pepper/FAPROTAX_1.2.12/
# Edit FAPROTAX_DIR in script if installed elsewhere

# Step 11 — βNTI/RCbray community assembly null model (~10–30 min)
# Requires rooted tree exported from QIIME2:
qiime tools export \
  --input-path version-2_integrated/rooted-tree.qza \
  --output-path version-2_integrated/exported_tree
Rscript 16_bnti_rcbray_analysis.R

# Step 12 — FAPROTAX × network hub cross-reference
python 17_faprotax_network.py

# Step 13 — ANCOM-BC differential abundance (genus level)
python 18_ancombc_plot.py
# Reads from ancombc_out/ancombc-genus-exported/

# Step 14 — IndVal indicator species
Rscript 19_indval.R

# Step 15 — Network robustness simulation
python 20_network_robustness.py

# Step 16 — Beta-dispersion (PERMDISP2)
Rscript 21_beta_dispersion.R

# Step 17 — Functional redundancy
python 22_functional_redundancy.py

# Step 18 — PGP bacterial index
python 23_pgp_index.py

# Step 19 — Sample filter sensitivity (ALL / ROOT_ZONE / STRICT)
python 24_sample_filter_sensitivity.py

# Step 20 — PERMANOVA + subsampling validation (1000 iterations, ~5–10 min)
Rscript 25_permanova.R

# Step 21 — N-cycle completeness + PGP specificity
python 26_ncycle_pgp_specificity.py

# Step 22 — Network threshold sensitivity (|r| = 0.3 / 0.4 / 0.5)
python 27_network_threshold_sensitivity.py

# Step 23 — Temporal stability Q1–Q4 (FAPROTAX metrics)
python 28_temporal_stability.py
```

---

## Output Files: Scripts 15–28

All files in `version-2_integrated/` unless otherwise noted.

| Script | Data output | Figure output |
|--------|-------------|---------------|
| **15** `15_faprotax_analysis.py` | `FAPROTAX_functional_table.tsv`<br>`FAPROTAX_group_comparison.csv`<br>`FAPROTAX_key_functions.csv`<br>`FAPROTAX_report.txt` | `Main_Fig_FAPROTAX.png/.pdf`<br>`Main_Fig_FAPROTAX_Log2FC.png/.pdf`<br>`Main_Fig_FAPROTAX_Ncycle.png/.pdf` |
| **16** `16_bnti_rcbray_analysis.R` | `BNTI_RCbray_results.csv`<br>`BNTI_RCbray_summary.csv` | `BNTI_RCbray_plot.png/.pdf` |
| **17** `17_faprotax_network.py` | `FAPROTAX_Network_CrossRef.csv` | `Main_Fig_FAPROTAX_Network.png/.pdf` |
| **18** `18_ancombc_plot.py` | `ancombc_out/ancombc-genus-exported/` | `Main_Fig_ANCOMBC.png/.pdf` |
| **19** `19_indval.R` | `IndVal_results.csv`<br>`IndVal_significant.csv` | `Main_Fig_IndVal.png/.pdf` |
| **20** `20_network_robustness.py` | `Network_Robustness_Results.csv`<br>`Network_Robustness_Summary.csv` | `Main_Fig_NetworkRobustness.png/.pdf` |
| **21** `21_beta_dispersion.R` | `BetaDispersion_results.csv`<br>`BetaDispersion_summary.csv` | `Main_Fig_BetaDispersion.png/.pdf` |
| **22** `22_functional_redundancy.py` | `Functional_Redundancy_Results.csv`<br>`Functional_Redundancy_Summary.csv` | `Main_Fig_FunctionalRedundancy.png/.pdf` |
| **23** `23_pgp_index.py` | `PGP_Index_Results.csv`<br>`PGP_Index_Summary.csv`<br>`PGP_Specificity_Results.csv` | `Main_Fig_PGP_Index.png/.pdf` |
| **24** `24_sample_filter_sensitivity.py` | — | `SampleType_Composition.png/.pdf` |
| **25** `25_permanova.R` | `PERMANOVA_results.csv`<br>`PERMANOVA_subsampling.csv` | `Main_Fig_PERMANOVA.png/.pdf` |
| **26** `26_ncycle_pgp_specificity.py` | `NCycle_Completeness_Results.csv`<br>`NCycle_Completeness_Summary.csv`<br>`PGP_Specificity_Results.csv` | `Main_Fig_NCycle_PGP.png/.pdf` |
| **27** `27_network_threshold_sensitivity.py` | `NetworkSensitivity_Results.csv` | `Main_Fig_NetworkSensitivity.png/.pdf` |
| **28** `28_temporal_stability.py` | `TemporalStability_Results.csv` | `Main_Fig_TemporalStability.png/.pdf` |

---

## Figure–Script Mapping

### Main Figures

| Figure | Panel | Script | Key result |
|--------|-------|--------|------------|
| **Fig 1A** | Overview schematic | — | Conceptual diagram |
| **Fig 1B** | Phylum barplot | `07_supp_figures.py` | Proteobacteria shift |
| **Fig 1C** | PICRUSt2 Log₂FC | `08_picrust2_analysis.py` | Nitrifier denitrification −11.5 |
| **Fig 1D** | Taxon–function heatmap | `09_taxon_function_corr.py` | Rhodanobacter–nitrification link |
| **Fig 1E** | Co-occurrence network | `10_network_analysis.py` | Hub collapse |
| **Fig 2A** | FAPROTAX Log₂FC | `15_faprotax_analysis.py` | Nitrification Log₂FC = −10.68 |
| **Fig 2B** | N-cycle completeness | `26_ncycle_pgp_specificity.py` | Space 73% vs. Terrestrial 100% |
| **Fig 2C** | Functional redundancy | `22_functional_redundancy.py` | Nitrification 0.20 vs. 4.80 genera |
| **Fig 2D** | FAPROTAX × network | `17_faprotax_network.py` | Nitrifiers as Terrestrial hubs |

### Supplementary Figures

| Supp Fig | Script | Description |
|----------|--------|-------------|
| **S1** | `13_supp_faith_pd.py` | Faith's PD: Space 10.05 vs. Terrestrial 24.98 |
| **S2** | `25_permanova.R` | PERMANOVA bootstrap validation (1000 iterations) |
| **S3** | `14_supp_temporal.py` | Q1–Q4 temporal stability (Faith's PD, Shannon) |
| **S4** | `16_bnti_rcbray_analysis.R` | βNTI/RCbray distribution: Terrestrial 93.2% deterministic vs. Space <0.1% |
| **S5** | `21_beta_dispersion.R` | Beta-dispersion: Space 0.574 vs. Terrestrial 0.398 |
| **S6** | `27_network_threshold_sensitivity.py` | Network metrics across |r| = 0.3/0.4/0.5 |
| **S7** | `28_temporal_stability.py` | Temporal stability of FAPROTAX metrics (KW p = 0.16–0.30) |

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
