
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

# Updating Controlled Vocabulary (CV) Files

**WARNING:** The Proteowizard project has been restructured to .NET/CLI. The 
C++ source files (cv.hpp, cv.cpp, CVTranslator.*) are **NO LONGER available**
from their GitHub repository at `https://raw.githubusercontent.com/ProteoWizard/pwiz/master`.
Using wget/curl to download these files will result in empty or corrupted files.

The CV files used by mzR can be updated using the `updateCvParams.sh` script:

```bash
cd src
./updateCvParams.sh
```

This script downloads only the `.obo` files (PSI-MS, UniMod, and Unit ontologies)
from their upstream repositories.

**IMPORTANT:** After updating the OBO files, the C++ source files `cv.hpp` and `cv.cpp`
should be regenerated using the cvgen tool. See
[`src/pwiz/data/common/REGENERATE_CV.md`](./src/pwiz/data/common/REGENERATE_CV.md)
for detailed instructions.

**Do NOT download CVTranslator.hpp, CVTranslator.cpp, or CVTranslatorTest.cpp** -
these are mzR-local files containing specific conflict-resolution logic that must
be preserved in version control.

# CV Version Compatibility

The `.obo` files and generated `cv.hpp`/`cv.cpp` files must be from compatible versions:

```bash
# Check psi-ms.obo version
grep "data-version:" src/pwiz/data/common/psi-ms.obo

# Check cv.hpp version
grep "data-version:" src/pwiz/data/common/cv.hpp
```

If you encounter compilation errors after updating, the cv.hpp/cv.cpp files may need
to be regenerated to match the new .obo file versions.

