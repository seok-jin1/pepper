# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Context

This is the bioinformatics pipeline for the paper *"Can rhizosphere microbiomes support plants in space? Microbiomic evidence for functional network collapse in ISS-grown Capsicum annuum"* (Kang & Shin, Plant Physiology Letters 2026). The pipeline compares ISS-grown *C. annuum* rhizosphere microbiomes (OSD-772, n=106 Space Flight) against terrestrial soil controls (PRJNA1145089, n=20).

**Current status:** Revised submission to Plant Physiology Letters (2026-03). Scripts were renumbered 06–17 after cleanup of supporting/exploratory scripts. See `RESEARCH_NOTES.md` for analysis notes.

---

## Environment Setup

```bash
# Core analysis (QIIME2, Python scripts 06–17, R scripts via Rscript)
conda env create -f envs/qiime2-amplicon.yml
conda activate qiime2-amplicon

# Functional prediction (PICRUSt2)
conda env create -f envs/picrust2.yml

# Optional: R-based SpiecEasi network
conda env create -f envs/pepper-network.yml
```

All Python scripts (06–17) and R scripts (13, 15) run in the `qiime2-amplicon` environment unless noted.

---

## Running Scripts

**All scripts must be run from the repository root** (the directory containing `config.py`). Scripts use paths relative to this root.

```bash
# Validate that all required input files are present
python config.py

# Quick-start: reproduce all figures from pre-processed data (no QIIME2 needed)
python download_processed_data.py   # downloads version-2_integrated/ from Zenodo
python 06_supp_figures.py           # Fig 1b
python 07_picrust2_analysis.py      # Fig 1c
python 08_taxon_function_corr.py    # Fig 2d
python 09_network_analysis.py       # Fig 2e
python 10_supp_faith_pd.py          # Supp Fig S1
python 11_supp_temporal.py          # Supp Fig S3
python 12_faprotax_analysis.py      # Fig 2a
Rscript 13_bnti_rcbray_analysis.R   # Supp Fig S4
python 14_functional_redundancy.py  # Fig 2c
Rscript 15_permanova.R              # Supp Fig S2
python 16_ncycle_pgp_specificity.py # Fig 2b
python 17_network_threshold_sensitivity.py  # Supp Fig S5

# Full pipeline (requires raw FASTQ downloads first)
bash 05_qiime2_pipeline.sh          # run after conda activate qiime2-amplicon
```

There are no unit tests or a test runner in this repository.

---

## Architecture

### Data Flow

```
Raw FASTQs (OSD-772 + PRJNA1145089)
  → 00–04: metadata prep, manifest, ID sync
  → 05: QIIME2 (DADA2 denoising, SILVA taxonomy, diversity metrics)
  → version-2_integrated/   ← all intermediate TSVs live here
  → 06–17: figure-generating Python/R scripts
  → version-2_integrated/ ← output PNGs/PDFs/CSVs
```

### Central Configuration: `config.py`

All file paths and analysis parameters for scripts 08–10 are centralized in `config.py`. Scripts import it as:

```python
import config
table = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t')
```

Key parameters in `config.py`:
- `RANDOM_SEED = 42` — must not change (peer review reproducibility)
- `NETWORK_PARAMS` — Spearman |r| > 0.4, raw p < 0.05
- `DIVERSITY_PARAMS` — rarefaction depth 1,000 reads (~86% sample retention)
- `PICRUST2_PARAMS` — pathway filtering thresholds

Scripts 00–06 and 10–11 use hardcoded relative paths rather than `config.py`. Scripts 07–09 and 12–17 import `config` for centralized paths.

### Sample Groups

The integrated dataset has three groups, but **only `Space_Flight` vs. `Terrestrial_Soil` are compared in the paper**. `Ground_Seed` (n=3, ISS ground controls) is excluded from all main analyses. Scripts using the combined table (`06_supp_figures.py`, `08_taxon_function_corr.py`) filter out `Ground_Seed` explicitly.

### Key Data Directories Under `version-2_integrated/`

| Directory | Contents |
|-----------|----------|
| `exported_table_clean/` | Combined ASV feature table (TSV) |
| `exported_table_space/` | Space Flight samples only |
| `exported_table_terrestrial/` | Terrestrial Soil samples only |
| `exported_taxonomy/` | SILVA 138 taxonomy assignments |
| `exported_diversity/` | Shannon, Faith's PD, ordination TSVs |
| `picrust2_out/` | PICRUSt2 pathway/EC/KO predictions |

### Figure–Script Mapping

| Figure | Script | Output file |
|--------|--------|-------------|
| Fig 1b (taxa barplot) | `06_supp_figures.py` | `Fig1B_Taxa_Barplot.png/.pdf` |
| Fig 1c (PICRUSt2 Log₂FC) | `07_picrust2_analysis.py` | `Fig1C_Functional_Pathways.png` |
| Fig 2a (FAPROTAX Log₂FC) | `12_faprotax_analysis.py` | `Fig2A_FAPROTAX_Log2FC.png/.pdf` |
| Fig 2b (N-cycle completeness) | `16_ncycle_pgp_specificity.py` | `Fig2B_NCycle_PGP.png/.pdf` |
| Fig 2c (functional redundancy) | `14_functional_redundancy.py` | `Fig2C_FunctionalRedundancy.png/.pdf` |
| Fig 2d (taxon–function heatmap) | `08_taxon_function_corr.py` | `Fig2D_Taxon_Function_Correlation.png/.pdf` |
| Fig 2e (co-occurrence network) | `09_network_analysis.py` | `Fig2E_Network.png/.pdf` |
| Supp Fig S1 (Faith's PD) | `10_supp_faith_pd.py` | `results/SuppS1_FaithPD.png/.pdf` |
| Supp Fig S2 (PERMANOVA bootstrap) | `15_permanova.R` | `SuppS2_PERMANOVA.png/.pdf` |
| Supp Fig S3 (temporal Q1–Q4) | `11_supp_temporal.py` | `results/SuppS3_Temporal_Q1_Q4.png/.pdf` |
| Supp Fig S4 (βNTI/RCbray) | `13_bnti_rcbray_analysis.R` | `SuppS4_BNTI_RCbray.png/.pdf` |
| Supp Fig S5 (network sensitivity) | `17_network_threshold_sensitivity.py` | `SuppS5_NetworkSensitivity.png/.pdf` |

### External Dependencies

- **FAPROTAX 1.2.12** — hardcoded path in `12_faprotax_analysis.py`: `/home/laugh/pepper/FAPROTAX_1.2.12`. Must exist locally; not bundled in repo.
- **SILVA 138 classifier** — must be trained once and placed at `classifiers/silva-138-99-classifier-515-806.qza` (see README Step 3).
- **Zenodo deposit** — processed data DOI is `https://doi.org/10.5281/zenodo.XXXXXXX` (placeholder; update before publication). Verify integrity with `md5sum -c processed_data.md5`.
