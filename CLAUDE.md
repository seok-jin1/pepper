# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Context

This is the bioinformatics pipeline for the paper *"Can rhizosphere microbiomes support plants in space? Microbiomic evidence for functional network collapse in ISS-grown Capsicum annuum"* (Kang & Shin, Plant Physiology 2026). The pipeline compares ISS-grown *C. annuum* rhizosphere microbiomes (OSD-772, n=106 Space Flight) against terrestrial soil controls (PRJNA1145089, n=20).

**Current status:** Desk-rejected by Plant Physiology (2026-03) for lacking mechanistic links. See `RESEARCH_NOTES.md` for the post-rejection analysis plan, FAPROTAX results, and resubmission targets.

---

## Environment Setup

```bash
# Core analysis (QIIME2, Python scripts 06–15)
conda env create -f envs/qiime2-amplicon.yml
conda activate qiime2-amplicon

# Functional prediction (PICRUSt2)
conda env create -f envs/picrust2.yml

# Optional: R-based SpiecEasi network
conda env create -f envs/pepper-network.yml
```

All Python analysis scripts (06–15) run in the `qiime2-amplicon` environment unless noted.

---

## Running Scripts

**All scripts must be run from the repository root** (the directory containing `config.py`). Scripts use paths relative to this root.

```bash
# Validate that all required input files are present
python config.py

# Quick-start: reproduce all figures from pre-processed data (no QIIME2 needed)
python download_processed_data.py   # downloads version-2_integrated/ from Zenodo
python 06_plot_depth.py             # Supp Fig 5
python 07_supp_figures.py           # Fig 1B + QC figures
python 08_picrust2_analysis.py      # Fig 1C
python 09_taxon_function_corr.py    # Fig 1D
python 10_network_analysis.py       # Fig 1E + keystone CSV
python 13_supp_faith_pd.py          # Supp Fig S1
python 14_supp_temporal.py          # Supp Fig S2
python 15_faprotax_analysis.py      # FAPROTAX functional trait analysis

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
  → 06–15: figure-generating Python scripts
  → results/ and version-2_integrated/ ← output PNGs/PDFs/CSVs
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

Scripts 00–07 and 11–15 use hardcoded relative paths rather than `config.py`.

### Sample Groups

The integrated dataset has three groups, but **only `Space_Flight` vs. `Terrestrial_Soil` are compared in the paper**. `Ground_Seed` (n=3, ISS ground controls) is excluded from all main analyses. Scripts using the combined table (`07_supp_figures.py`, `09_taxon_function_corr.py`) filter out `Ground_Seed` explicitly.

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

| Figure | Script |
|--------|--------|
| Fig 1B (taxa barplot) | `07_supp_figures.py` |
| Fig 1C (PICRUSt2 Log₂FC) | `08_picrust2_analysis.py` |
| Fig 1D (taxon–function heatmap) | `09_taxon_function_corr.py` |
| Fig 1E (co-occurrence network) | `10_network_analysis.py` |
| Supp Fig S1 (Faith's PD) | `13_supp_faith_pd.py` |
| Supp Fig S2 (temporal Q1–Q4) | `14_supp_temporal.py` |
| FAPROTAX (new, post-rejection) | `15_faprotax_analysis.py` |

### External Dependencies

- **FAPROTAX 1.2.12** — hardcoded path in `15_faprotax_analysis.py`: `/home/laugh/pepper/FAPROTAX_1.2.12`. Must exist locally; not bundled in repo.
- **SILVA 138 classifier** — must be trained once and placed at `classifiers/silva-138-99-classifier-515-806.qza` (see README Step 3).
- **Zenodo deposit** — processed data DOI is `https://doi.org/10.5281/zenodo.XXXXXXX` (placeholder; update before publication). Verify integrity with `md5sum -c processed_data.md5`.
