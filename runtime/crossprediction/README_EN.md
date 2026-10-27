# Saved-model B/C crossprediction extension v0.1

This self-contained **public release** package replays the frozen maize B/C training-representation by B/C test-input diagnostic from the 120 original saved models. It contains all136 consumed source files, including the120 saved joblib artifacts, eight unchanged predictor/reference snapshots, folds, protocols, original prediction tables and model manifests. No connection to the original data tree, credentials or network is needed. It does not fit any model, acquire native observations or resolve source/crop/registration holds.

## Execute from an arbitrary directory

Use the recorded CPython3.11.13 Linux environment and `requirements_inference.txt`. Commands do not install dependencies:

```console
python portable_replay.py --mode validate
python portable_replay.py --mode guardtest --out /absolute/new/guardcheck
python portable_replay.py --mode replay --out /absolute/new/replayed_result --start-method spawn
```

The output folder must be new and outside the immutable package. `spawn` is the default and is explicitly exercised by the validation. Available process-start modes depend on the operating system. Cross-platform byte identity is not promised; the supplied original models were tested in the recorded Linux environment. The public guard blocks external research-data/config reads and network connections, independently of any original source-server pathname. A positive-control probe checks that denial, and all120 inference-task records carry a verified worker guard. An independent fresh-worker guard test performs no prediction or fitting.

The adapter checks the manifest before inference. It calls the public-path-adapted original inference/scoring implementation in `original_code/frozen_crossprediction.py`, with separately hashed `RELATIVE_INPUT_BINDINGS.json`. `original_code/CPU_INPUT_BINDINGS.json` supplies logical source IDs in this public export; it is not executed. `SOURCE_TO_PACKAGE_LEDGER.json` records historical source identities using logical IDs. Source models remain byte-identical; metadata-sanitized inputs/configuration have separately verified public hashes. The historical analysis protocol had SHA256 `756e27f0ebf641a73af0bf04c2efa7811a43ab04ff2c9ad91716160cde31680d`. The public derivative replaces private locators and therefore has a different file SHA; the original hash remains an as-of provenance identity, not a claim of current byte identity.

The replay produces112,240 scored records:77,280 from the principal supported roles and34,960 from the separate earlier-year Mead hindcast. These are repeated representations of642 or437 physical dates, not112,240 independent observations. All five fixed families, both fitting and scoring references and all supported original axes remain present. The old diagonals are replayed; the off-diagonals deliberately combine different training/test representations for a conditional computational diagnostic. This does not identify physical footprint improvement, an additive causal decomposition, a forecast or a new preferred operational pipeline.

`reference_results/` retains the original frozen outputs. After replay, the adapter compares every prediction and all six scientific score/contrast tables with their reference frames. Runtime/freezing timestamps, package-relative model paths and added guard fields legitimately differ in execution logs. Verification receipts outside the immutable package document the actual comparison outcomes. A corrupted-copy test should fail the manifest check before inference.

## Scope and access

The archive includes the exact saved models to avoid an upstream refit dependency. Models are used only after their expected SHA256 and lengths have been verified. As with any pickle/joblib artifact, run only a trusted, verified copy. The author approved public distribution under the release file-level rights map; registration/upload is a separate step. The acquired data citations and delivered policy notices in `docs/` must be retained. They do not establish author research funding or coauthorship.

Native EC partitioning, new MODIS/AgERA5 acquisition/QA, daily dynamic FFP, physical crop-purity validation and the broader2,431-case fits are outside this extension. Private author/funding/COI worksheets, correspondence, raw email headers, credentials and unrelated manuscripts/PDFs are excluded. `plot_crossprediction.py` is an optional data-derived plotting script; install the separately recorded plotting requirements if using it, and provide the generated or reference result folder with `--results` plus a new `--out` folder. The tested figure exports are not a universal cross-platform image-hash guarantee.


## Public release adaptation

This is a separately hashed public derivative of the archived internal component. Paths and provenance use logical IDs; original estimator bytes and numerical algorithms remain unchanged. Historical protocol hashes identify the original fitted objects; sanitized protocol files have their own new SHA. The generic offline I/O guard replaces machine-specific paths. Use the top-level release README and file-level rights map for approved MIT/CC-BY4.0 scope. Older as-of statements in archived provenance documents describe their original execution and do not assert public registration.
