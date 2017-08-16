
# Continuous integration

The mzR project master branch is subjected to CI using travis: 

[![Build Status](https://travis-ci.org/sneumann/mzR.svg?branch=master)](https://travis-ci.org/sneumann/mzR)

The devel branch is tested on multiple architectures by the BioC build farm:

http://bioconductor.org/checkResults/devel/bioc-LATEST/mzR/

# Supported compilers

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

# Compatibility of older mzR versions

Please remember that the Bioconductor project has a rule to fix the versions of packages
to be installed by `biocLite()` to the version of R, e.g. Ubuntu 16.04.2 ships with R version 3.2.3
which is quite old by Bioconductor standards. `biocLite("mzR")` will then try to install
http://bioconductor.org/packages/3.1/bioc/html/mzR.html which is also nearing its second birthday...
Note that this does not only apply to mzR, but all BioC packages. 

If the installation fails with an error similar to `error: no matching function for call to ‘call_once()`
as reported in several issues here, this requires a newer version of mzR,
e.g. from https://cran.r-project.org/ (for Ubuntu see e.g. https://cran.r-project.org/bin/linux/ubuntu/)
which will in turn use a recent BioC release, which will give you a recent
mzR. 




