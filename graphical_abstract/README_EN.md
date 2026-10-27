# Graphical Abstract

This single diagram summarizes the comparison roles in the maize GPP benchmark.
It accompanies the article "Training–Test Representation and Reference-Path Sensitivity in a Satellite-Based Maize Gross Primary Production Benchmark".

- `gagraphic.png`: 660 x 295 pixels, 300 dpi metadata; compact image rendition.
- `gagraphic.pdf`: one-page vector master with embedded font subset.
- `gagraphic_preview.png`: 2640 x 1180 pixels; larger rendition of the same diagram.
- `GRAPHICAL_ABSTRACT_CAPTION_EN.txt`: full definitions, result scope and interpretation limits.

Rows denote training representation and columns denote test input. B uses uniform
cells; C uses normalized monthly retained-core weights. The principal random-forest
callout compares CC versus BB and then contrasts C versus B training at fixed B or
C test input. It does not imply independent GPP truth, a unique causal effect,
complete daily footprint correction or general prediction superiority.

The principal and hindcast cohorts remain distinct. The displayed percentage is
the principal equal-site-season own-fit/own-score VUT/CUT pipeline comparison;
the caption distinguishes the 100 principal and 20 hindcast saved models.

Project-created artwork/caption is CC BY 4.0. See `LICENSE_NOTICE.txt` for the
font exception. Original source data attribution remains in the main record.

The drawing source and requirements are included here. To render with a
licensed Arial-compatible TTF and Poppler installed:

```console
python graphical_abstract/create_graphical_abstract.py --out ../maize-ga-rendered --font /path/to/arial-compatible.ttf --pdftoppm /path/to/pdftoppm
```

Font software is not included. Geometry and labels are fixed; font/library
changes may alter rasterization. Rendering does not execute models or acquire data.
