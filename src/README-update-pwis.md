## Updating PWIZ

- Check out a copy of https://github.com/ProteoWizard/pwiz somewhere
- The pwiz version is in ./pwiz/Version.cpp, e.g. 3_0_20294

- Move old stuff away:
  cp -avx pwiz pwiz-old
  git rm -rf pwiz


## Updating CVs

**WARNING** (since 2025): The Proteowizard project has been restructured to .NET/CLI.
The C++ files (cv.hpp, cv.cpp, CVTranslator.*) are **NO LONGER available** from 
`https://raw.githubusercontent.com/ProteoWizard/pwiz/master`. Attempting to download 
them with wget will result in empty files or JSON responses, causing build failures.

### The updateCvParams.sh Script

The `updateCvParams.sh` script (in the repository root) downloads only the `.obo` files:

```bash
#!/bin/bash
# Download only .obo files from their upstream sources
# These are the actual ontology files that need updating

mkdir -p src/pwiz/data/common

# PSI-MS ontology from HUPO-PSI
wget -N -O src/pwiz/data/common/psi-ms.obo \
  https://github.com/HUPO-PSI/psi-ms-CV/raw/master/psi-ms.obo

# UniMod ontology
wget -N -O src/pwiz/data/common/unimod.obo \
  https://www.unimod.org/obo/unimod.obo

# Unit ontology (UO)
wget -N -O src/pwiz/data/common/unit.obo \
  https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/uo.obo
```

### CVTranslator Files Must Be Preserved

**DO NOT** attempt to download CVTranslator.hpp, CVTranslator.cpp, or 
CVTranslatorTest.cpp from external sources. These files are mzR-local
modifications containing conflict-resolution logic for CV term mappings.
They should be maintained in version control and NOT regenerated or downloaded.

### Regenerating cv.hpp and cv.cpp

If you need to regenerate cv.hpp and cv.cpp from updated .obo files, see:
[`pwiz/data/common/REGENERATE_CV.md`](./pwiz/data/common/REGENERATE_CV.md)