# Maize GPP figure reproduction

This component reproduces the four main and seven supplementary figures from source-bound aggregate data and the archived native geometry. It performs no model fitting, prediction or scientific rescoring. The original plotting bodies retain model/reference priorities, adverse cases, axes, colors and the fixed July geometry choice. Paths, entrypoints, public figure numbering and additional PDF export were adapted.

Use the recorded Linux CPython3.11.13 environment and `requirements_figures.txt`. From this directory:

```console
python run_figures.py --out ../maize_figures_replay
```

The output must be a new directory outside this immutable component. The program verifies all43 input files, writes11 figures in PNG/SVG/PDF and records a replay receipt. It reads no original server tree and needs no network credentials. Relative paths are resolved from the script location, not the shell's working directory.

`data/` contains the37 published aggregate CSVs, the20-row own-path conditional table, the exact-rational omission summary, source geometry and three original footprint raster subsets. The path-only privacy edits in TableS23 replace84 source-path cells with logical source IDs. Every other table cell, including numerical strings, dates, flags, categories and source hashes, is unchanged. Geometry path metadata is similarly replaced by logical identifiers; numeric coordinates, affine/WKT parameters, weights and source-map hashes are preserved. See `PATH_SANITIZATION.json`, `SOURCE_BINDINGS.json` and `FIGURE_INPUT_MANIFEST.json`.

`reference/` preserves11 reviewed original PNGs and11 prepared PDF assets. The replay was executed on the recorded CPU environment. Nine regenerated PNGs match the reference RGB pixels exactly; FiguresS3 andS4 have matching dimensions and source-bound unchanged plotting values with minor rendering differences, and were visually checked at original resolution. Pixel identity is not promised for arbitrary systems or future library versions. Original reference artwork remains available even when environment-specific antialiasing/layout differs.

The main and supplementary captions are supplied separately. Figure1 is an explanatory diagram drawn with the original explicit specification; Figure2 preserves native georeferencing and the Ne1 coordinate discrepancy. The monthly cores are conditional normalized published cores, not complete daily FFP or field/crop-purity validation. The principal642-case and separate437-date populations remain separate;480 hindcast training cases reuse principal data. Original comparison priorities and all adverse cases are retained.

Project-owned plotting source falls under the release MIT code license. Project-created numerical/figure outputs follow the approved CC-BY4.0 scope, while the three native Chu raster subsets and other provider-derived material retain their source rights, citations and attribution. Source DOI/version and full provider/software/font notices are supplied in the final release rights materials. No library binary, private correspondence, internal author worksheet, credential or unrelated manuscript is included.

This is a prepared release component. No public repository or project DOI is asserted. Actual registration/upload remains subject to the author's final review.
