
# Continuous integration

The mzR project master branch is subjected to CI using travis: 

[![Build Status](https://travis-ci.org/sneumann/mzR.svg?branch=master)](https://travis-ci.org/sneumann/mzR)

The devel branch is tested on multiple architectures by the BioC build farm:

http://bioconductor.org/checkResults/devel/bioc-LATEST/mzR/

# Installation

## Installation on macOS

If you install from source, you first need to install 
the netCDF headers and libraries: `brew install netcdf`, 
and then instal mzR the BioC way via:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("mzR", suppressUpdates = TRUE)
```

## Installation on Linux

We are not shipping the full set of boost headers due to 
size restrictions. This *might* cause compilation failure 
due to missing files. We have tested several OS and compiler 
build environments successfully, please report any compilation failures
at https://github.com/sneumann/mzR/issues
and we'll add the missing files. 

mzR-2.9.1 with boost-1.59.0 has been tested on the following compilers:

On the BioC build farm (http://bioconductor.org/checkResults/devel/bioc-LATEST/mzR/)
* malbec2 (Ubuntu 16.04.1, gcc-5.4.0)
* tokay2 (Windows Server 2012 R2, MinGW-W64-4.9.3)
* Failling oaxaca (Apple clang 3.5svn / 600.0.57)

Also on:
* Debian stretch/sid (gcc-5.4.0)
* Ubuntu 16.04 (clang 3.8.0)
* Ubuntu 16.04 (gcc-5.4.0)
* Ubuntu 14.04 (gcc-4.8.2)
* Ubuntu 12.04 (gcc-4.6.3)

# Contributions

Please read the [contributions guide and code of conduct](./CONTRIBUTIONS.md).

# bio.tools metadata

`inst/biotools.xml` is a bio.tools submission document generated from DESCRIPTION,
`inst/CITATION` and roxygen2 annotations. It describes raw-data and identification
reading, spectrum writing, and metadata-preserving conversion. It does not
register or update an entry on bio.tools automatically.

From the repository root, install the suggested `roxygen2` and `xml2` packages
and run `Rscript tools/update-biotools.R`. This metadata-only command does not
load or compile mzR. The custom roclet is also included in DESCRIPTION's
`Roxygen` configuration, so normal `roxygen2::roxygenise()` runs update the XML
alongside the existing documentation. The development-only implementation lives
in `tools/biotools.R`; it adds no runtime dependencies to mzR.

Annotate public functionality in its roxygen block, for example:

```r
#' Read mass spectrometry raw data
#'
#' Open a file and expose its spectra and metadata.
#' @edam_operation Data parsing (operation_1812)
#' @edam_input_data Mass spectrometry data (data_2536)
#' @edam_input_format mzML (format_3244)
#' @edam_output_data Mass spectrometry data (data_2536)
#' @noRd
```

Each block with `@edam_operation` becomes one bio.tools function. Repeating the
operation tag adds alternative operations to that function. Formats belong to
the preceding input/output data tag in the same block; repeat a data tag to
start another input/output, rather than forming an ambiguous cross product.
Use `@edam_output_format` for file outputs. The title, description and details
become the function note. `@noRd` allows annotations alongside hand-maintained
Rd manuals without replacing them.

Package-wide annotations live in `R/mzR.R`:

| Tag | Value |
| --- | --- |
| `@edam_topic` | EDAM label and identifier, e.g. `Proteomics (topic_0121)` |
| `@biotools_tool_type` | bio.tools type, e.g. `Library` |
| `@biotools_operating_system` | `Linux`, `Windows` or `Mac` |
| `@biotools_language` | Language, e.g. `R` or `C++` |
| `@biotools_link` | `Type \| URL`, e.g. `Software catalogue \| https://bioconductor.org/packages/mzR/` |
| `@biotools_documentation` | `Type \| URL`, e.g. `User manual \| https://example.org/manual` |
| `@biotools_id` | Optional registered bio.tools identifier; omit until assigned |

Use `Package documentation | URL` for the package reference manual. Since
biotoolsSchema has no `Package documentation` type, the roclet emits the
schema-supported `General` type with `Package documentation` as its note.

These tags can be repeated except for the identifier. Topics and package-wide
values are deduplicated. EDAM tags require `Label (branch_NNNN)` with the
appropriate topic, operation, data or format branch. The underscore spellings
also work with roxygen2 7.x through registered prefix parsers.

Name, description, homepage, version, license, author credits (including email,
roles and ORCID), repository and issue-tracker links, and available publication
DOI/PMID/PMCID identifiers are extracted automatically. CITATION entries without
any of those identifiers are not emitted because biotoolsSchema requires one.
Supported platforms and languages are explicit annotations, not guesses based
on the machine running roxygen.

Generation validates the complete XML against the pinned upstream schema in
`tools/biotools/` before replacing the output. Invalid tags, orphan formats,
missing operations or invalid schema values fail explicitly. Run
`Rscript tools/test-biotools.R` for the metadata regression checks. Updating the
XML prepares a submission; publishing it still requires bio.tools registration.
