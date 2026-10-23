## Updating PWIZ

- Check out a copy of https://github.com/ProteoWizard/pwiz somewhere
- The pwiz version is in ./pwiz/Version.cpp, e.g. 3_0_20294

- Move old stuff away:
  cp -avx pwiz pwiz-old
  git rm -rf pwiz


## Updating CVs

```
#!/bin/bash
set -x
cd src

PWIZGITHUBPREFIX=https://raw.githubusercontent.com/ProteoWizard/pwiz/master

CVFILES="pwiz/data/common/CVTranslator.cpp \
  pwiz/data/common/CVTranslator.hpp \
  pwiz/data/common/CVTranslatorTest.cpp \
  pwiz/data/common/cv.cpp pwiz/data/common/cv.hpp \
  pwiz/data/common/psi-ms.obo"

for F in $CVFILES ; do
  wget -O $F $PWIZGITHUBPREFIX/$F
done
```
