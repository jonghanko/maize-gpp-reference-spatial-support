# Third party notices and data attribution

Prepared27October2026. This release separates original project source code, project-created data/model artifacts, source-provider material and separately installed software. The exact file-level rights map in the release manifest determines the applicable license; the project MIT license does not relicense third-party material.

## Original data sources

The seven dataset records and their persistent identifiers are in DATA_SOURCES_AND_LICENSES.csv and DATA_SOURCES_AND_LICENSES.json. Full dataset citations must remain in the article's citation-tracked References as well as in the release documentation. The four acquired FLUXNET-1F notices identify CC-BY4.0 and are retained verbatim in provider_notices/. Source licenses/attribution apply to redistributed or transformed source material; identify the extraction, QA, aggregation, interpolation, unit conversion or other modifications. No source creator endorses the project by being cited.

Required FLUXNET acknowledgment reproduced from delivered policy version20251226:

FLUXNET data products were produced and harmonized by eddy covariance regional networks and data processing centers, including AmeriFlux, ChinaFlux, European Fluxes Database, ICOS, JapanFlux, KoFlux, OzFlux, SAEON, and TERN. These products also include a modified version of ERA5 hourly data provided by the Copernicus Climate Change Service.

Required AmeriFlux data-service acknowledgment:

Funding for the AmeriFlux data service was provided by the U.S. Department of Energy Office of Science.

The FLUXNET data products used in this work were downloaded from the FLUXNET Shuttle on26October2026. This acquisition date applies to the four stated FLUXNET products, not to satellite, weather or footprint acquisition.

Recommended processing citation: Pastorello, G., Trotta, C., Canfora, E., et al.(2020). The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Scientific Data7,225. https://doi.org/10.1038/s41597-020-0534-3.

MOD09A1 Collection6.1 is distributed by NASA LP DAAC; retain its product DOI and NASA source credit. NASA materials retain their original provider status and any product-specific notices. AgERA5 v2.0 is provided through the ECMWF Climate Data Store on behalf of the Copernicus Climate Change Service. The used footprint dataset is the published Chu et al. supplementary release, version20200622, under CC-BY4.0. Their official policy/metadata links are in the source table. The project's extracted or modified products must be identified as such.

Ro1 and IAC creator strings follow the acquired provider citation metadata, including the supplied Tim/Timothy entries and unusual IAC name fields. This release makes no unsupported correction of provider identities. US-IAC appears only in the separately reported supporting point-input context; it is not an additional site in the principal spatial comparison or Mead hindcast.

## Research funding and source-service support

The research author reported support from Chonnam National University(Grant No.2026-1417-01) and the Basic Science Research Program through the National Research Foundation of Korea(Grant No.NRF-2021R1A2C2004459). These research grants are separate from the required AmeriFlux data-service DOE acknowledgment above. Source-service support does not establish that the author received a DOE research grant or that data providers are coauthors.

## Separately installed software

The project uses the following installed-distribution versions, read from the recorded Linux execution environment on27October2026. No Python interpreter, third-party package binary or environment image is distributed by this source/data release. The entries summarize package metadata; full installed-distribution notices, including bundled component exceptions and font notices, are retained under third_party_licenses/ and inventoried in SOFTWARE_DEPENDENCY_LICENSES.json.

| Distribution | Recorded version | Installed metadata summary |
|---|---|---|
| numpy | 2.4.6 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |
| pandas | 2.3.3 | BSD 3-Clause License |
| pyarrow | 25.0.0 | Apache-2.0 |
| scipy | 1.17.1 | Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers. |
| scikit-learn | 1.9.0 | BSD-3-Clause |
| joblib | 1.5.3 | BSD-3-Clause |
| threadpoolctl | 3.6.0 | BSD-3-Clause |
| matplotlib | 3.11.1 | License agreement for matplotlib versions 1.3.0 and later |
| rasterio | 1.4.4 | BSD |
| pyproj | 3.7.2 | MIT |

SciPy and rasterio metadata include original copyright/license documents rather than a single complete SPDX expression; their full supplied notices govern. NumPy/pandas/Matplotlib and other wheels may carry component-specific notices; a short table label does not replace them. Matplotlib's DejaVu/STIX font notices are included for retained source graphics and embedded/outlining font use. No font is claimed as original project work or relicensed under MIT/CC-BY4.0.

If a future distribution includes third-party wheels, shared libraries, containers or vendored source, preserve that exact distribution's complete licenses and additional notices; this inventory does not clear an unexamined future binary bundle. User-installed dependencies remain governed by their respective licenses.

## Scientific and access scope

Source openness does not resolve crop-operation, tower-registration, GPP-reference dependence or footprint-support limitations. Provider collaboration under the acquired CC-BY4.0 policy is encouraged, not mandatory. An unanswered inquiry is not itself a prohibition on source-data redistribution. The source DOI/license obligations remain independent of authorship and scientific approval.

No private correspondence, raw email headers, credentials, author confirmation worksheets, unrelated unpublished manuscripts or full internal server snapshot is part of the public payload. No public accession, uploaded repository or issued project DOI is asserted by these prepared notices.
