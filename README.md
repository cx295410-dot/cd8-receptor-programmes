# CD8 receptor programmes during cancer immunotherapy

Analysis and figure code for **Conserved CD8 receptor programmes accompany divergent transcriptional remodelling during cancer immunotherapy**.

The study compares receptor-independent transcriptional programmes with receptor RNA, paired surface-protein measurements, observed TCR sharing and contemporaneous pathology. The primary treatment analysis contains 57,303 CD8 cells from 25 paired patients across basal cell carcinoma, lung cancer and oral cancer. Independent breast RNA/ADT measurements contain 5,061 CD8 cells from 5 paired tissue donors. A separate resection resource contains 151,601 CD8 cells from 15 patients; the primary pathology analysis uses 13 patients and 125,887 cells.

## Repository contents

| Directory | Purpose |
| --- | --- |
| `pipelines/00_source_data` | Source acquisition and preparation helpers inherited from the source project |
| `pipelines/01_target_panel` | Fixed receptor extraction, donor-level operations and protein comparisons |
| `pipelines/02_transcriptome` | Full-transcriptome graphs, receptor-independent programmes and paired pseudobulk inference |
| `pipelines/03_expansion` | Additional marker, clone and measured functional-label comparisons |
| `pipelines/04_clinical` | Primary-resection eligibility, download, fixed programme analysis and checks |
| `pipelines/05_multipanel` | Full-gene pathology follow-up, post hoc sensitivity and 8 main / 26 supplementary composite figures |
| `config` | Target panel, programme definitions, analysis protocols and recorded software versions |
| `tools` | Safe staging of scripts into a configurable local study layout and code-package checks |
| `provenance` | Source-script hashes and preparation checks |

## Public sources

Original data remain in GEO: [GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813), [GSE179994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994), [GSE200996](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200996), [GSE267552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267552), [GSE176021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176021), [GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220), [GSE91061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061), and [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222). Source publications and input roles are listed in [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md).

This repository contains code and analytical definitions. The separate [aggregate data release v1.0.1](https://github.com/cx295410-dot/cd8-receptor-programmes/releases/tag/data-v1.0.1) contains group-, gene- and pathway-level estimates, eligible summary matrices and the complete clinical endpoint results. Download both Supplementary Data ZIP files and the README; manifests and SHA256 checksums are provided. The published version is immutable and corresponds to analysis commit `5ecda74e61bb768b6c52d344b028aa84bb108db9`. It has a version-specific public URL and no assigned DOI.

The public data release excludes patient-, donor-, cell- and clonotype-level derived tables, individual outcome records and individual expression/count matrices. It does not redistribute original article PDFs, author contact forms, authentication material or third-party packages. The aggregate deposit alone does not reproduce every individual-level analysis or plot; use the cited original public sources and the code documentation for the required local intermediate tables.

## Start here

1. Create a Python environment and install the needed scientific packages. `requirements-recorded.txt` records versions used in the study; `requirements.txt` lists unpinned package names for environment setup. Compatibility with other versions needs validation.
2. Read [docs/REPRODUCING.md](docs/REPRODUCING.md) before running analysis. The original scripts use a historical study directory layout. `tools/stage_project.py` copies them into a selected workspace and substitutes the two source-directory prefixes without changing analytical logic. It refuses to overwrite an existing file with different content.
3. Acquire the public input files and original annotations/supplements described in the documentation. Acquisition helpers are provided, but this is not an automated accession-to-figure workflow and some inputs require manual retrieval.
4. Run the documented stages in dependency order. Publication figures are produced by the scripts in `05_multipanel` after their input tables and models exist.

```powershell
python tools/check_repository.py
python tools/stage_project.py --study-root D:/CD8_Reanalysis/Virtual_Biotech_Study --legacy-root D:/CD8_Reanalysis/source_project --dry-run
```

The staging command performs no analysis and downloads no data. Remove `--dry-run` to create the script/configuration layout after selecting your directories.

## Interpretation

Patients define biological replication. Memory, Dysfunction and the other programme names are operational transcript-feature labels. Clone sharing is observed repertoire overlap, not a differentiation trajectory. Experimental MANA/viral labels originate from the source study and are not inferred for untested clones. The primary pathology endpoint is residual viable tumour percentage at surgery. Separate exploratory baseline bulk analyses examine objective response in GSE78220 and GSE91061, overall survival in GSE78220 and progression-free survival in GSE135222. None of the 24 programme–endpoint tests, or the separate eight-test Memory/Dysfunction family, met BH q < 0.05; these are association analyses and do not establish a validated clinical prediction model. The Gender-adjusted genome-wide model is a separately specified post hoc sensitivity. Complete testing families and nonsignificant results should be retained when reproducing the study.

## Citation and licensing

Use the repository URL and the exact commit identifier when referring to this code; `CITATION.cff` supplies author and title metadata. The associated manuscript has no assigned publication DOI in this release. A project reuse license has not yet been selected. Original datasets and third-party software retain their own terms.

AI assistance supported programming, numerical checks and plotting code. Numerical estimates and plots are computed from recorded source data. Human authors are responsible for the scientific interpretation and submitted article.
