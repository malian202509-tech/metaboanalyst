# Placental Oxylipin Metabolomics Analysis

This repository contains Python scripts and documentation for a placental oxylipin metabolomics analysis workflow.

## Scope

Included in this repository:

- `scripts/`: analysis scripts for cohort alignment, feature filtering, imputation, transformation, and batch-effect diagnosis.
- `docs/`: workflow notes, SOP-style documentation, directory description, and script index.
- `results/tables/`: small audit/result tables needed to document current analysis progress.

Excluded by default:

- `data/`: raw and processed data files.
- generated figures and other large outputs under `results/`.

## Current Script Flow

```text
scripts/00b_align_omics_to_cohort.py
  -> scripts/01_feature_filtering.py
  -> scripts/02_imputation.py
  -> scripts/03_log_pareto.py
  -> scripts/04_batch_effect_diagnosis.py
```

`scripts/00_baseline_table1.py` generates the clinical baseline Table 1 and is maintained alongside the omics preprocessing workflow.

For detailed usage, see:

- `docs/脚本索引.md`
- `docs/数据处理与分析进展.md`
- `docs/目录结构说明.md`
