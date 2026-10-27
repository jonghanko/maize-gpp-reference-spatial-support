# Maize GPP reference and spatial support

Software version 1.0.2 accompanies **Training–Test Representation and Reference-Path Sensitivity in a Satellite-Based Maize Gross Primary Production Benchmark**. The companion data record supplies processed observations, original model objects and analysis results.

## Citable version and reproducibility pin

- Software v1.0.2: [doi:10.5281/zenodo.23689895](https://doi.org/10.5281/zenodo.23689895).
- Processed data, original model objects and results v1.0.2: [doi:10.5281/zenodo.23689996](https://doi.org/10.5281/zenodo.23689996).
- [GitHub v1.0.2 release](https://github.com/jonghanko/maize-gpp-reference-spatial-support/releases/tag/v1.0.2), pinned source commit `f89593c17501d18e810d00d5629022d277dec1c5`.

For exact reproduction, use the version-specific records and tagged source. The attached `maize-gpp-code-v1.0.2.zip` remains unchanged: SHA-256 `b1ed403e1c71062e2d442c5ecf4338a75cc9176795573abeb21bd301dcae85ca`. This default-branch metadata update adds the published identifiers; the scientific runtime and original release tag/archive are preserved. Cite the software/data version records together with the original source datasets.

The principal comparison contains 642 site–date observations (489 distinct calendar dates); the separate Mead hindcast contains 437 site–date test observations (291 distinct calendar dates). The 480 later-year training observations reuse principal cases. Archived `n_days` fields count site–date observations within the reported cohort/role, not distinct dates across towers. See `OBSERVATION_COUNT_CLARIFICATION.md`. No schema or numerical values were renamed or recomputed.

## Assemble and reproduce

Use Linux CPython 3.11.13 with the recorded requirements for scientific replay/refit. Install dependencies in a separate environment; the analysis commands do not install software or acquire data. Download the companion data record and retain its directory structure:

```console
python tools/assemble_runtime.py --data ../maize-gpp-reference-spatial-support-data --out ../maize-gpp-runtime
python ../maize-gpp-runtime/crossprediction/portable_replay.py --mode validate
python ../maize-gpp-runtime/crossprediction/portable_replay.py --mode replay --start-method spawn --out ../maize-cross-results
python ../maize-gpp-runtime/baseline/run_reproduce.py --mode replay --out ../maize-baseline-results
python ../maize-gpp-runtime/baseline/run_reproduce.py --mode refit --workers 4 --out ../maize-baseline-refit
python tools/render_publication_figures.py --runtime ../maize-gpp-runtime --out ../maize-publication-figures
```

Use fresh output directories outside immutable input components. Assembly verifies both full release manifests. Scientific inference validates component manifests, blocks network access and confines input reads. Only load the supplied trusted models after their recorded SHA checks pass. The optional 220 baseline tasks reproduce existing fits; they are not new scientific experiments.

## Current and historical artwork

`runtime/publication_figures12/` reproduces the current four main and seven supplementary figures at 7.16in width with all plot text at least 12pt. Panel descriptions are outside the artwork in its current `CAPTIONS_EN.md`. The module includes only its minimal 15 numerical/geometry inputs, source and provider/software/font notices. Display dependency environments are recorded separately from model-execution versions. A single CPU style environment reproduced all 11 plots; numerical values and geometry are the target, not future bitwise graphics identity.

The companion data record's `runtime/publication_artwork12/` contains the 11 current reference PDFs and their SHA map. The original `runtime/figures/` module and earlier data-record reference artwork are retained as historical replay targets, not labelled as the latest publication layout. To reproduce that earlier layout, use `python ../maize-gpp-runtime/figures/run_figures.py --out ../maize-historical-figures`.

## Scientific boundaries

The package supports processed-input score replay, 120-model B/C crossprediction, the original 200 principal plus 20 hindcast reproduction tasks, 37 publication CSVs and current figure regeneration. Future platform/library bit identity is not promised. The 2,431-case context is an archived aggregate comparison rather than a completely reconstructed fit pipeline. Native EC partitioning, every satellite/weather acquisition/QA step, all 84 native footprint rasters, new FFP and independent geographic/forecast validation are outside the portable replay. Three original native raster subsets support the spatial illustration. Monthly conditional cores do not establish complete daily flux mass or crop purity. Ro1 crop dates and Ne1 registration remain unresolved. The hindcast uses later-year labels at already represented sites; cross-input comparisons remain post-review exposed-domain diagnostics.

## Provenance, credit and rights

Original model bytes, features, targets, dates, splits, metric and aggregation roles remain unchanged. Public logical IDs replace private source locators; historical original protocol SHA values embedded in estimators remain separate from sanitized public hashes. `runtime/PUBLIC_EXPORT_RELATION.json` is the historical v1.0.0 export receipt, not the current version declaration. `PUBLICATION_ALIGNMENT.json` describes this artifact version and publication-artwork additions. Certain archived original-engine/protocol descriptions retain the earlier Rice software lineage; they remain unchanged for source identity. Active release entrypoints and case bindings describe the Maize study. Historical source labels do not add Rice observations or Rice validation to this article.

Project-owned source is MIT. Project-created data/model/figure artifacts and documentation retain CC BY4.0 scopes; source datasets, software and fonts retain their original rights and attribution. Consult `FILE_LEVEL_LICENSES.json` and `runtime/licenses/`. Provider citations are not replaced by a project DOI. No dependency binaries, cached fonts, author photograph, manuscript Word files, private correspondence or credentials are included. Creator/ORCID, 14 confirmed CRediT roles and the two confirmed grants are in `CITATION.cff`, `CREDIT.md` and `FUNDING.md`.

## Graphical Abstract

`graphical_abstract/` contains the matching compact PNG, vector PDF, larger PNG and
English caption. It is one summary diagram in three renditions, distinct from the
11 publication figures. The caption scopes the principal random-forest percentage
and the conditional training contrasts, and distinguishes the two cohorts. The
module is outside `runtime/`, so it adds no scientific replay input or model change.
Project-created diagram/caption rights are CC BY4.0; font exceptions remain explicit.
