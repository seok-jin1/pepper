# Research Notes — Post-Rejection Analysis Plan

**Status:** Plant Physiology pre-proof submitted → desk rejection (2026-03)
**Rejection reason (verbatim):** *"While the paper presented some interesting correlations, the mechanistic links are missing. Another concern is that the work is difficult for others to reproduce, because of the specific conditions required."*

---

## 1. Reproducibility (addressed)

The rejection cited specific experimental conditions (microgravity) as a reproducibility concern — not the code. This is inherent to spaceflight research and cannot be resolved by software changes alone. However, the GitHub repository has been updated to:

- Include all processed intermediate files (`version-2_integrated/`) for direct analysis reproduction without re-running QIIME2
- Fix the SILVA classifier setup in `05_qiime2_pipeline.sh` (trained classifier is now correctly generated via `fit-classifier-naive-bayes`)
- Add `download_processed_data.py` for Zenodo-based data retrieval (Zenodo DOI pending upload)

---

## 2. Mechanistic Links — Root Cause

The paper shows *correlations* (co-occurrence networks, PICRUSt2 predicted pathways) but lacks *mechanistic* evidence for how the microbiome change affects plant physiology. Plant Physiology expects a clear chain: **microbiome change → functional mechanism → plant outcome**.

---

## 3. Plant Phenotype Data — Search Results

### OSD-772 (GLDS-675)
- **Contains:** 16S rRNA amplicon sequencing only (430 FASTQ files + metadata)
- **Does NOT contain:** leaf area, fresh weight, biomass, harvest yield, or any growth metrics
- Focus was food safety microbiology, not plant phenotyping

### LSDS Accession Search
- **Result: No separate LSDS-prefixed dataset for PH-04 phenotype data exists publicly as of March 2026**
- The OSDR API lists OSD-772 as the only dataset tagged with project identifier PH-04
- No companion LSDS dataset has been deposited

### What Was Measured but Not Deposited
From NTRS ID 20240015344 (Monje et al., ASGSRB 2024, San Juan — conference presentation):

| Measurement | Source |
|-------------|--------|
| Canopy CO₂ drawdown / gross photosynthesis rate | APH onboard sensors |
| Canopy Quantum Yield (CQY) | APH onboard sensors |
| Leaf/canopy temperature (IR) | APH IR transducer |
| Stomatal conductance (inferred) | Monje et al. 2024 |
| Germination delay (~14 days vs ground) | NASA press + Monje et al. 2024 |
| Fruit count (26 fruits, Day 109 + Day 137) | Khodadad et al. 2026 |
| Microbial CFU/gram fresh weight | Khodadad et al. 2026 (published) |

Destructive measurements (fresh weight, dry weight, leaf area, root biomass) were not deposited. Half the fruit was consumed by crew; the other half was used for microbial analysis only.

### Contact for Unpublished APH Sensor Data
**Oscar Monje** (oscar.a.monje@nasa.gov, Aetos Systems / Kennedy Space Center)
PI of the plant physiology component of PH-04. The CO₂/CQY/temperature data from the 2024 NTRS presentation have not been formally deposited. A data-sharing request may be appropriate as a collaboration.

---

## 4. FAPROTAX Analysis — Completed (script 15)

`15_faprotax_analysis.py` maps ASVs to experimentally validated functional traits using FAPROTAX 1.2.12 (Louca et al. 2016, *Science*). Unlike PICRUSt2, FAPROTAX is based on curated literature — not genome-scale inference.

### Key Results

| Function | Space Flight (%) | Terrestrial (%) | Log2FC | Interpretation |
|----------|-----------------|-----------------|--------|----------------|
| **nitrification** | 0.004 | 7.34 | **-10.7** | 99.9% collapse in space |
| **sulfate_respiration** | 0.000 | 0.091 | -16.5 | complete loss |
| **dark_sulfide_oxidation** | 0.000 | 0.054 | -15.7 | complete loss |
| **denitrification** | 0.148 | 0.010 | +3.9 | enriched in space |
| **nitrogen_fixation** | 4.74 | 0.66 | +2.9 | enriched in space |
| aerobic_chemoheterotrophy | 40.6 | 30.0 | +0.44 | slightly enriched |
| aromatic_compound_degradation | 1.53 | 1.92 | -0.32 | slightly reduced |

### Mechanistic Interpretation

The nitrogen cycle is disrupted in a specific, directional way:

1. **Nitrification collapses** (NH₄⁺ → NO₂⁻ → NO₃⁻ pathway, -10.7 log2FC) — the taxa responsible for converting ammonium to plant-available nitrate are nearly absent in spaceflight
2. **Denitrification increases** (+3.9 log2FC) — remaining fixed nitrogen is converted to N₂ gas and lost
3. **Net effect:** The rhizosphere loses its capacity to supply NO₃⁻ to the plant, while simultaneously accelerating nitrogen loss from the system

This provides a **literature-validated mechanistic link** (not PICRUSt2 prediction):
*Spaceflight → loss of nitrifying taxa → nitrification collapse → disrupted N supply → potential N limitation of plant growth*

If APH sensor data (canopy CO₂, CQY) shows reduced photosynthesis during the same period, the full chain would be:
**nitrification collapse → N deficiency → reduced photosynthesis → impaired growth**

---

## 5. Additional Analyses to Consider

### Priority 1 — Strengthen mechanistic narrative (feasible with existing data)

| Analysis | Tool | What it adds |
|----------|------|--------------|
| **βNTI / RCbray null model** | `picante` (R) | Quantifies deterministic vs stochastic assembly; "microgravity acts as strong environmental filter" → mechanistic |
| **Neutral Community Model (NCM)** | `sloan_fit` (Python/R) | Fraction of community explained by neutral drift; low R² in space = strong selection |
| **LEfSe or ANCOM-BC differential abundance** | `qiime2` plugin | Statistically rigorous taxa enrichment, replaces visual comparison |
| **Indicator species analysis** (IndVal) | `indicspecies` (R) | Identifies taxa specifically diagnostic for spaceflight condition |

### Priority 2 — Plant physiological connection

| Analysis | Requirement | What it adds |
|----------|------------|--------------|
| Correlate nitrification capacity with APH CO₂/CQY | APH sensor data from Monje et al. | Direct microbiome → plant physiology link |
| FAPROTAX results × keystone genus co-occurrence | Current data | "Hub taxa collapse coincides with functional collapse" narrative |

### Priority 3 — Resubmission target assessment

| Journal | Fit | Required additions |
|---------|-----|-------------------|
| **Microbiome** (BioMed Central) | High | Current + FAPROTAX sufficient |
| **ISME Journal** | High | Current + null model assembly analysis |
| **mSystems** (ASM) | High | Current sufficient |
| **Plant Physiology** (resubmit) | Conditional | Need APH sensor data (Monje contact) + βNTI |
| **Frontiers in Microbiology** (Space section) | High | Current sufficient |

---

## 6. Key References

- Louca et al. (2016) Decoupling function and taxonomy in the global ocean microbiome. *Science* 353:1272–1277. [FAPROTAX original paper]
- Khodadad et al. (2026) Evaluating microbial community profiles of Chile peppers grown on the ISS. *Scientific Reports* DOI: 10.1038/s41598-025-20440-9
- Monje et al. (2024) Plant Habitat-04: Characterizing pepper responses to environmental changes during spaceflight. ASGSRB 2024, NTRS ID: 20240015344
- OSD-772: https://osdr.nasa.gov/bio/repo/data/studies/OSD-772
