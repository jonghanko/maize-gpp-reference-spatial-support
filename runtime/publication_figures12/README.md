# Publication figures at12pt

This source-only addendum reproduces the four main and seven supplementary figures in the article-associated12pt layout. It belongs to the companion publication-artifact bundle. The original v1.0.0 numerical records and source archive remain unchanged. 

## Contents and scope

- Three layout programs: main Fig.1/2, Fig.3/4/S7, and Fig.S1-S6; one unchanged native-geometry helper; one optional independent contrast verifier; a portable runner.
- Exactly15 numerical/source inputs: geometry JSON plus3 native retained-core GeoTIFFs,2 contrast CSVs, and9 supplementary JSON/CSV inputs. A filtered manifest lists the9 SI inputs; their bytes and original hashes are unchanged.
- The authoritative current captions are in `CAPTIONS_EN.md`; older caption-migration notes are a historical display-change ledger. Editable caption-migration notes, per-file SHA256 and license records, seven existing source citations/provider notices, and software/font license texts.
- Recorded reference-artwork checks and CPU portability checks. Generated PDFs/PNGs/SVGs, fitted models, predictions, dependency binaries, cached fonts, photos, Word files and official templates are excluded.

These programs only redraw already-frozen numerical results and geometry. They do not fit models, run model predictions, acquire data, create footprints or select new scientific cases. All figure canvases are exactly7.16in wide, all actual PDF text states are at least12pt, PDF/SVG artwork is vector and normal PNG exports use600dpi. Keep this width when inserting figures; do not shrink the figures to force a caption onto a page.

## Portable execution

Use a separate Python environment. The single CPU reference environment was smoke-tested for all11 figures; its display dependencies are pinned in `requirements-cpu-style.txt`. Installing dependencies is a separate user action.

```sh
python -m pip install -r requirements-cpu-style.txt
python run_figures.py --out output
```

Choose a fresh output directory. Rendering individual groups is also possible:

```sh
python run_figures.py --groups main --out output_main
python run_figures.py --groups contrasts --out output_contrasts
python run_figures.py --groups supp --out output_supp
```

Optional contrast verification needs the separate dependencies in `requirements-optional-checks.txt` and the complete generated contrast outputs:

```sh
python run_figures.py --groups contrasts --out output_checked --verify-contrasts
```

The SI renderer may create a Matplotlib cache inside its output directory. No cache is shipped in this archive. Relative defaults and explicit `--source`/`--out` options replace task-specific path bindings. The native geometry helper and scientific plot bodies are preserved; all adapter changes are listed in `docs/PORTABLE_PATCH_LEDGER.json`.

## What was checked

The unpacked-style CPU replay rendered all11 figures, retained source-input hashes, exactly reproduced3 native geometry records and all96 displayed conditional/aggregate values, and reproduced every serialized SI plot value with maximum difference0.0. The original SI reviewer independently checked216 source-table values. All11 CPU output PDFs retain exact7.16in width and actual font states of at least12pt; SI text-bounding boxes report zero intersections or outside-canvas text. The original final previews were separately visually reviewed. Checks are in `checks/`.

The reference PDF hashes are identifiers for files supplied with the manuscript, not files bundled here. Cross-version byte-identical graphics are not promised: timestamps, font serialization and Matplotlib versions may affect bytes. Numerical equivalence, geometry, dimensions and font checks are the reproducibility target.

## Versions and scientific records

`STYLE_ENVIRONMENTS.json` distinguishes the final display environments: F1/F2 used CPU Python3.11.13, Matplotlib3.11.1, pandas2.3.3 and NumPy2.4.6; F3/F4/S7 used Python3.11.15, Matplotlib3.11.1, pandas3.0.5 and NumPy2.4.6; S1-S6 used bundled Python3.12.14, Matplotlib3.11.2, pandas3.0.1 and NumPy2.3.5. Optional reference requirement files retain those plotting pins. The source-only replay used the CPU display environment for all11 figures. These records describe typography/layout work and do not replace or change the historical baseline/crossprediction model-execution environment. The existing upstream software notices remain historical notices of their recorded distribution, supplemented by the Matplotlib3.11.2 license/font texts. No packages or font binaries are bundled.

## Attribution and licensing

Original project source follows the MIT license in `LICENSE-CODE-MIT.txt`; project-created tables, documentation and figures follow the project's CC-BY4.0 scope in `LICENSE-DATA.md`. Third-party sources keep their existing rights. The3 native footprint raster subsets retain the Chu etal. dataset's CC-BY4.0 attribution. The seven source records, source DOI/version strings, four delivered AmeriFlux policy notices, required acknowledgments and font/software notices remain in `licenses/upstream/`. US-IAC's citation concerns the separately reported broader context and does not add a tower to the principal comparison. These notices do not imply source endorsement, a new release accession or a new scientific approval.

`FILE_LEVEL_LICENSES.json` maps this addendum's files. The historical license inventories inside `licenses/upstream/` describe the existing broader candidate; they do not inventory this new archive. Preserve both current and upstream notices. No third-party material is relicensed by this addendum.
