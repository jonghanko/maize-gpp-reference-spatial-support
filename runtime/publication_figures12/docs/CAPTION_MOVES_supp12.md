# Supplementary figure text moved to captions

All descriptive titles and notes listed below have been removed from plotted artwork. Essential model/site/arm identifiers, panel letters, axes, legends and numeric labels remain at least 12 pt at exactly 7.16 inches width. No model fits, predictions, source values, axis limits, color limits or case selections changed.

## Figure S1

Panel mapping:
- single_panel: 20-cell C-minus-B omission heatmap; same row and column order.

Move or retain these statements in the caption:
- C − B sensitivity to site-season omissions from the error summary
- Fixed LOYO out-of-fold outputs; blue means lower C error and orange higher C error.
- Evaluation-summary omissions only: no retraining, new validation or uncertainty intervals.
- Omitting Ro1 leaves only Mead in evaluation; its labels remain in the original training folds.
- The full column contains four site-seasons; each omission leaves three.

## Figure S2

Panel mapping:
- (a): Random forest, primary.
- (b): Ridge, secondary.
- (c): Per-site-season C-minus-B percentages, same 3×4 order.

Move or retain these statements in the caption:
- Earlier-year Mead hindcast: train 2012/2013 and test 2010/2011.
- 480 training days; 437 test days; one existing region.
- Panel (c) season counts: Ne1 2010 n=128, Ne1 2011 n=148, Ne3 2011 n=161.
- VUT/CUT represent paired fit/score reference paths.
- Negative percentages indicate lower RMSE; dependent processed GPP references, fixed monthly 80% DAY cores; retrospective, not a regional or forecast test.

## Figure S3

Panel mapping:
- (a): Principal cohort C minus B.
- (b): Principal cohort B minus shared-scene point.
- (c): Earlier-year Mead hindcast C minus B.

Move or retain these statements in the caption:
- The three original panels are stacked vertically; scientific order and shared limits are unchanged.
- Orange bars are squared-mean-error contributions, blue bars centered-error contributions, and diamonds their sum.
- Both own-reference paths are retained for Ridge, random forest and secondary histogram gradient boosting.

## Figure S4

Panel mapping:
- (a): Weight concentration Nc in monthly maps used by retained-date anchors.
- (b): Realized input differences between C and B on retained dates.

Move or retain these statements in the caption:
- Panels are stacked vertically for readability.
- Month is the actual composite-start month of a retained-date anchor.
- P denotes principal 642 dates; T1 denotes separate 437-date hindcast.
- Blank light-grey months are unused by retained-date anchors.
- The quantity Nc is conditional concentration cell count, not an independent sample size or field-purity metric.

## Figure S5

Panel mapping:
- (a): C minus B input contrast.
- (b): B minus shared-scene point input contrast.
- (c): C minus shared-scene point input contrast.

Move or retain these statements in the caption:
- Dependent reference paths: own-fit/own-score input interactions.
- All five families are retained under the equal-site-season estimand.
- The B-minus-shared-point contrast was added after review.
- Both fitting and scoring reference change together; this does not isolate uStar or establish a causal interaction.

## Figure S6

Panel mapping:
- (a): Ridge C minus B.
- (b): Random forest C minus B.
- (c): Ridge B minus shared-scene point.
- (d): Random forest B minus shared-scene point.
- (e): Ridge C minus shared-scene point.
- (f): Random forest C minus shared-scene point.

Move or retain these statements in the caption:
- All four LOYO site-seasons: directions vary by model and comparator.
- Numeric readouts beneath each panel show VUT in the upper blue row and CUT in the lower orange row; all original 48 numbers and bar heights are retained.
- Site-season sample counts, shared across all panels: Ne1 2012 n=156; Ne1 2013 n=167; Ne3 2013 n=157; Ro1 2015 n=162.
- Shaded Ro1 has identical test inputs across arms; its exact 2015 operation dates remain unknown.
- The same Ro1 train/test membership recurs in LOYO-2015, partial LOSO and joint exclusion.
- No uncertainty intervals are estimated.
