# Caption moves for Fig. 3, Fig. 4 and Fig. S7

These are display-only replacements. All source values, cohorts, model/reference identities, colors, common axis limits and contrast definitions are retained. No source CSV, manuscript or shared figure file map was edited. The output width is exactly 7.16 inches. Every visible PDF character has a nominal 12-point font, including Unicode superscripts; the floating-point extraction minimum is 11.999999999999986 points (12 points to numerical precision). There is no superscript exception.

## Fig. 3: unchanged 2 × 2 panel mapping

- (a): principal LOYO, Ridge; 642 days, four site-seasons, two regions.
- (b): principal LOYO, random forest (RF); the same 642 days, four site-seasons, two regions.
- (c): separate Mead hindcast, Ridge (secondary); 437 days, three site-seasons, one existing region.
- (d): separate Mead hindcast, RF; the same 437 days, three site-seasons, one existing region.

Verbatim text removed from the artwork and assigned to the caption:

```text
Conditional representation contrasts from saved models
(a) Principal LOYO · Ridge
(b) Principal LOYO · RF
(c) Mead hindcast · Ridge
(secondary)
(d) Mead hindcast · RF
642 days · 4 site-seasons
2 regions
437 days · 3 site-seasons
1 existing region
Negative conditional/pipeline contrasts indicate lower error under C. Interaction is not an error gain.
```

The principal/hindcast subtitle occurred once per relevant panel; the removed caption information applies to both model panels in its row. Keep the existing caption's explanation that negative conditional/pipeline differences favor C, that interaction is not an error gain, and that RMSE nonlinearity and aggregation prevent unique additive causal attribution.

Caption integration: replace the original row/column mapping sentence with “Panels show (a) principal Ridge, (b) principal RF, (c) separate Mead-hindcast Ridge (secondary), and (d) separate Mead-hindcast RF. The principal comparison retains 642 days from four site-seasons in two regions; the hindcast retains 437 days from three site-seasons in one existing region.” The remaining scientific caption content is retained.

Short essential axis-label key, for the caption: “Train C−B | B” and “Train C−B | C” compare training representations at fixed B and C test input, respectively; “Test C−B | B” and “Test C−B | C” compare test inputs within B- and C-trained models, respectively. “Interaction” and “Pipeline” retain the existing formula definitions. Horizontal numerical tick labels suppress the leading zero only; their positions and values remain −0.15, −0.10, −0.05, 0, +0.05 and +0.10. Their fixed limits remain −0.15 to +0.10 g C m⁻² d⁻¹.

## Fig. 4: new letters for the unchanged 2 × 2 panels

The source had four model/aggregation titles without panel letters. Added letters identify the same axes without changing their arrangement:

- (a): Ridge, mean of four site-season RMSEs.
- (b): RF, mean of four site-season RMSEs.
- (c): Ridge, mean of two region-pooled RMSEs.
- (d): RF, mean of two region-pooled RMSEs.

Verbatim text removed from the artwork and assigned to the caption:

```text
Complete year exclusion: input contrasts under two aggregation estimands
Ridge ·
Mean of four site-season RMSEs
RF ·
Mean of four site-season RMSEs
Ridge ·
Mean of two region-pooled RMSEs
RF ·
Mean of two region-pooled RMSEs
Same 642 cases. Negative values favour the first arm.
B − A_shared is a post-review descriptive contrast.
The alternative changes within-region pooling and regional weighting;
C − A_shared remains conditional on point registration.
```

Caption integration: replace the original row/column description with “Panels show (a) Ridge and (b) RF under the mean of four site-season RMSEs, and (c) Ridge and (d) RF under the mean of two region-pooled RMSEs.” Retain the existing 642-case count, post-review status, pooling/weighting distinction, negative-sign definition and absence of population uncertainty bars. Add “C−A_shared remains conditional on point registration. The blue upper and orange lower numerical rows below each axis report the exact four-decimal VUT and CUT values, respectively.”

All 24 four-decimal labels remain. They are positioned in two reference-colored readout rows below each axis rather than at adjacent bar ends; bars and their heights are unchanged. The “shared” line belongs to A in both B−A_shared and C−A_shared. The essential RMSE quantity and units are horizontal above the panels at 12 points rather than rotated along the left edge. The original symmetric y-limit formula and resulting limits are unchanged.

## Fig. S7: unchanged 1 × 2 panel mapping

- (a): principal LOYO, HGB; 642 days, four site-seasons, two regions.
- (b): separate Mead hindcast, HGB; 437 days, three site-seasons, one existing region.

Verbatim text removed from the artwork and assigned to the caption:

```text
Conditional representation contrasts from saved models
(a) Principal LOYO · HGB
(b) Mead hindcast · HGB
642 days · 4 site-seasons
2 regions
437 days · 3 site-seasons
1 existing region
Negative conditional/pipeline contrasts indicate lower error under C. Interaction is not an error gain.
```

Caption integration: add “Panels show (a) principal LOYO, retaining 642 days from four site-seasons in two regions, and (b) the separate Mead hindcast, retaining 437 days from three site-seasons in one existing region.” Retain HGB's secondary role and the existing statements that the populations/references remain separate and the displayed families were not selected using the new crossprediction scores. Use the same short-axis key and unchanged conditional/interaction definitions as Fig. 3.

## Delivery and checks

Final artwork stems in this folder are `Figure3_BC_conditional_contrasts`, `Figure4_primary_LOYO_input_effects` and `FigureS7_HGB_conditional_contrasts`, each supplied as embedded-type PDF, live-text SVG and 600-dpi PNG. Heights are 6.7, 7.0 and 4.2 inches, respectively. The `_final_width_100dpi.png` files are unscaled 716-pixel visual previews representing the final insertion width.

`INDEPENDENT_VALUE_FONT_BBOX_VISUAL_QA.json` independently recovers all 48 + 24 + 24 plotted values from SVG geometry, checks original CSV/export equality, embedded PDF fonts and every visible PDF character at 12 points, exact width and maximum height, 600-dpi dimensions, text and axes bounding boxes, and the final-width visual inspection. All pass. No new fitting or prediction occurred.
