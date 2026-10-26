# Reproducing the analytical stages

## Workspace and environment

The recorded code used two directory prefixes: a study root containing `biology_rebuild_20261026`, and a legacy source-project root containing the original GEO downloads and `analysis_extensions`. The prefixes are configurable through `tools/stage_project.py`. Directory suffixes retain their historical names so existing file dependencies stay explicit. Historical directory names do not indicate additional analysed cohorts or experiments.

The stage tool copies files and rewrites these directory prefixes in Python string literals. Original repository scripts remain unchanged. It also installs the fixed receptor panel and programme definitions into their expected local locations. It does not download data, install dependencies, regenerate analysis protocols or assert that every input is present. Do not execute scripts until the required input files exist. Several plotting modules read tables at import time.

Use a dedicated Python environment. Compare installed versions with `config/software_versions_recorded.json`; the recorded environment was assembled on Windows. Memory and storage requirements depend on the downloaded matrices. The included `stream_rds_matrix.py` is a source-specific helper, not a general RDS reader.

## Dependency order

1. **Source resources.** The staged legacy `scripts` directory contains download and preparation helpers. Read each helper before use; some expect previously downloaded publication supplements or source file lists. Retain source-specific names and directory structures. These helpers are not a complete unattended acquisition pipeline.
2. **Receptor and protein inputs.** In the study `scripts` directory, use `extract_targets.py`, `analyze_targets.py` and `protein_validation.py` after their original matrices and annotations are present. `evidence_core.py` supplies patient-level operations. The target-panel smoke checks can run without downloading the datasets.
3. **Whole transcriptomes.** In `biology_rebuild_20261026/scripts`, run `rebuild_biology.py extract`, `process`, `de` and `summarise` in that order. `process` and `de` optionally take cohort names. Follow with `biology_sensitivity.py`, the source-specific export repair where required, and the plotting/verification scripts. Inspect the script entry points and existing outputs before deciding whether a repair step is needed.
4. **Expanded source tables.** In `reference_expansion_20261026/scripts`, run `expansion_analysis.py`, then the corresponding graph/figure scripts after the base tables exist. The functional-antigen comparison depends on the source-project antigen preparation and labels.
5. **Pathology extension.** In `clinical_validation_20261026/scripts`, use `clinical_download.py`, eligibility checks and `clinical_analysis.py`; then clinical plotting and verification. `config/clinical_analysis_protocol_frozen.json` records the earlier fixed programme analysis. The download helper requires the original supplementary clinical workbook and source annotations; acquisition is not inferred from expression.
6. **Whole-gene pathology and final plates.** In `multipanel_rebuild_20261026/scripts`, run `largepanel_clinical_genome.py`, followed by the separately documented `largepanel_clinical_sex_sensitivity.py`. The primary genome script also runs its verification; `--verify` repeats verification using existing results. Then draw `largepanel_main_1_4.py`, `largepanel_main_5_8.py`, `largepanel_supp_1_15.py`, `largepanel_supp_16_25.py` and `largepanel_supp_26.py`. Run the supplied figure-source verification scripts after output generation.

## Expected checks

Primary treatment cohorts contain 12,364 BCC, 36,172 lung and 8,767 oral CD8 cells, with 11,8 and6 patient pairs. The independent protein context contains 5,061 cells from 5 paired tissue donors. There are 8 main and 26 supplementary parent figures with 263 total panels in the completed local study.

The external clinical count extraction contains 15 patients and 33,514 recorded genes. The primary analysis excludes the two functional-label donors and retains 13 patients and 14,129 tested genes. The fixed source-prefix-adjusted model gives 31 gene-level BH discoveries; the unadjusted and post hoc Gender-adjusted sensitivities give 18 and 27. These are checkpoints for the recorded inputs and software, not numbers to enforce by changing filters or tuning models. Complete inferential families and raw-count conservation checks take precedence over a visual match.

This release has been checked for Python syntax, byte-preserving source copies, absence of disallowed file types/known credential patterns, and target-panel helper smoke checks. The upload preparation does not rerun the full biological analysis or constitute independent validation of the scientific findings. Full result and figure-source archives are maintained separately from this code-only repository.
