# Maize GPP processed-input reproduction package v0.1

This public release component reproduces the two fixed maize GPP experiments and the later error/input-support audits from supplied processed inputs and frozen predictions. It does not reacquire or reprocess native Flux, ONEFlux, MODIS or AgERA5 sources, validate physical daily footprints, resolve site-operation records, or add new scientific cases/experiments. No public repository accession or comprehensive redistribution clearance is implied.

## Run from any location

Use CPython 3.11 with the recorded dependencies. Default replay needs NumPy and pandas only; it does not need PyArrow, a GPU, source-server access, credentials or network access. No installation is performed by these commands.

```console
python run_reproduce.py --mode validate
python run_reproduce.py --mode guardtest --out /absolute/new/guard-test/folder
python run_reproduce.py --mode replay --out /absolute/new/output/folder
python run_reproduce.py --mode refit --out /absolute/new/refit/folder --workers 4
python run_reproduce.py --mode all --out /absolute/new/full/folder --workers 4
```

Output folders must be new. `requirements_replay.txt`, `requirements_refit.txt` and `CURRENT_RUNTIME_SNAPSHOT.json` state the tested runtime; the snapshot is not an unsupported assertion that every historical package version was recorded at original fitting. The default replay produces full fitting-reference × scoring-reference matrices and native pooled/fold/region/site-season metrics plus macro summaries; it compares the available original metric tables. Its `scores/principal/macro_summary.csv` is the reproduced numeric matrix supporting Figure2 interpretations. It also reproduces and compares TablesS21–S27 through explicitly identified CSV/path adapters of the archived recipes. Original plotting code/figures are retained as reference artifacts; default replay does not claim to regenerate every manuscript figure.

The optional refit runs the exact original 200 principal and 20 T1 family/arm/reference tasks, using the immutable original engine, saved protocol parameters, exact positional fold memberships and supplied CSV inputs. These are **reproduction refits**, not new experiments. It compares all 38,640 + 8,740 own-reference prediction records to frozen values and reports model byte hashes only where actually compared. Current numeric comparison tolerance is 1e-10 absolute/relative; it is not a universal cross-platform bitwise guarantee.

## Cohorts and inputs

Principal: 642 daily physical cases from four site-seasons/two regions, four optical arms, two dependent NT references and five fixed families. Five supported fold instances yield 200 tasks; six unsupported planned strata remain visible. Its full cross-reference table has 77,280 repeated scoring records. T1: 480 reused later-year Mead training cases and 437 earlier-year test cases, two optical arms/two references/five families and 20 tasks. Its cross table has 17,480 repeated records. Neither repeated rows nor reused training cases create independent regional replication.

All 12 snapshots and both prediction tables are CSV.gz with 17-significant-digit numeric export; all numeric columns were checked for exact round-trip values during construction. CSV readers parse daily dates without UTC/DST conversion and preserve row order, because the original fold indices are positional. The original eight-predictor order and training-only scaling/model settings remain in the archived engine. `ARTIFACT_MAP.json` contains portable relative execution paths. `SOURCE_EXPORT_RECEIPT.json` links each copied/exported artifact to the original source hash; original code/protocols remain distinct from the new portable adapters/protocols.

The executable outputs retain their historical table IDs. For manuscript v0.8, historical Tables S21–S27 correspond to publication Tables S20–S26; historical T1 S18/S18b correspond to publication S17a/S17b. `PUBLICATION_NUMBERING_CROSSWALK.json` records these aliases without renaming frozen source outputs. Figure IDs are unchanged.

## Explicit native-source boundary

All 84 monthly selected-cell weight records, accepted-composite audits and daily optical inputs are supplied. S24–S27 can therefore be replayed without native TIFFs. Native map-byte validation is **off by default and reported not performed**, not silently passed. If a separately obtained set of the 84 exact original TIFFs is available, set `MAIZE_NATIVE_MAP_ROOT` to their flat containing folder for optional byte-hash validation. The package does not inspect raster bands or reproduce native spatial integration/acquisition. The published monthly retained 80% cores are conditional support, not complete daily flux footprints or crop purity.

Execution installs a generic offline I/O guard in the main process and every refit worker initializer: external research-data/config reads and network access are denied outside the component, explicit output and Python runtime. Worker positive controls verify denial, and `--mode guardtest` separately tests a fresh spawned process without refitting. Private original locations are represented by logical source IDs; they are not execution dependencies. The source tree may still exist on the test machine, so the validation receipt reports blocked/unused access rather than claiming that the tree was physically hidden.

The 2,431-case contextual point benchmark has two aggregate tables in `context_archive/`, copied and checksummed only; they are not regenerated by this minimum package. Rice/German readiness and synthetic agent work are also outside executable reproduction. Author/title/funding/COI worksheets, correspondence, email diagnostics, credentials and unrelated under-review texts are excluded. Data citations/notices in `docs/` must be retained; approved public-license scopes are recorded in the release rights materials.

`MANIFEST.json` binds immutable package files. Validation occurs before replay/refit; a deliberately changed payload in a separate scratch copy must fail validation. Reproduction receipts distinguish validation, frozen-output replay,220 reproduction refits, native-raster checks and the absence of raw acquisition/new scientific eligibility.


## Public release adaptation

This is a separately hashed public derivative of the archived internal component. Paths and provenance use logical IDs; original estimator bytes and numerical algorithms remain unchanged. Historical protocol hashes identify the original fitted objects; sanitized protocol files have their own new SHA. The generic offline I/O guard replaces machine-specific paths. Use the top-level release README and file-level rights map for approved MIT/CC-BY4.0 scope. Older as-of statements in archived provenance documents describe their original execution and do not assert public registration.
