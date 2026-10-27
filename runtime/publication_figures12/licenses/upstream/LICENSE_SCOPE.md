# License scope

Prepared27October2026 for a reviewable public-release candidate. The author approved MIT licensing for original project software and CC-BY4.0 for project-created data, figures and saved-model artifacts. Actual upload/publication remains subject to the author's review of the final payload.

## Original project source code

The project's MIT license applies to original project source scripts, adapters, tests and configuration material identified as project-owned source in the release file-level rights manifest. It does not relicense imported libraries, installed packages, third-party license texts, source datasets, fonts, logos or other third-party assets. Existing source notices must be retained. Citation to a scientific method does not itself grant a license to copied third-party code; any subsequently discovered externally copied source must keep its original license or be excluded until its provenance is resolved.

## Project-created data and models

CC-BY4.0 applies to project-created processed tables, prediction/score/contrast outputs, project documentation/figures and saved trained-model state identified as such in the file-level manifest. Retain author/title/version/project identifier once assigned and indicate modifications. This grant covers only the project's rights; upstream provider attribution and license conditions remain effective. A model file does not license its separately installed Python implementation or embedded third-party rights under a different license.

## Third-party data and fonts

The four acquired AmeriFlux FLUXNET-1F products, AgERA5 and the published footprint dataset retain their CC-BY4.0 terms and source credit. NASA MOD09A1 original materials retain the applicable NASA/LP DAAC provider status and any specific source notices; the project does not impose exclusive copyright in them. Extracted or transformed copies must be labelled and linked to source DOI/version. Provider notices and exact required acknowledgments are in THIRD_PARTY_NOTICES.md and provider_notices/.

Matplotlib/DejaVu/STIX and other software/font materials retain their original notices; see third_party_licenses/. Scientific graphics made with those tools are project figures, while the font software itself is not project-owned code or data. No third-party library binaries or full environment are included in this candidate.

## Manifests and release boundaries

The public candidate must have a fresh file manifest and file-level license map after path/privacy sanitization. Frozen internal packages and their old manifests remain unchanged. Logical source identifiers, original source SHA256 values and dataset identifiers preserve provenance without publishing private internal paths. Public candidate replay/score/figure checks must be completed before upload. A prepared or reserved identifier is not a registered public DOI; insert actual version-specific project repository identifiers only when they exist.

The release reproduces the documented processed-data analyses. It does not claim native eddy-covariance partitioning or every upstream satellite/weather acquisition step, independent field validation, complete reconstruction of the separate2431-case context fits, or portable bit identity under arbitrary future libraries/platforms.
