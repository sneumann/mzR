## Updating PWIZ

- Check out a copy of https://github.com/ProteoWizard/pwiz somewhere
- The pwiz version is in ./pwiz/Version.cpp, e.g. 3_0_20294

- Move old stuff away:
  cp -avx pwiz pwiz-old
  git rm -rf pwiz


## Updating CVs

ProteoWizard retired its C++ tree in pwiz commit 6dd40d8c93 (ProteoWizard/pwiz#4658),
so `cv.cpp` and the `cvgen` generator are no longer on pwiz master.
`../updateCvParams.sh` downloads the latest psi-ms, unimod and unit `.obo` files,
builds `cvgen` from the last C++ pwiz commit (a9ed4ce960) and regenerates
`pwiz/data/common/cv.hpp` and `cv.cpp`. It needs a C++17 compiler and Boost.

If a term used in the C++ sources gets renamed or obsoleted upstream, the build
fails with an unknown `MS_...` identifier; rename it in the sources accordingly.
