
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
source("https://bioconductor.org/biocLite.R")
biocLite("mzR", suppressUpdates = T)
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
