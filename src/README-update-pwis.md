## Updating PWIZ

- Check out a copy of https://github.com/ProteoWizard/pwiz somewhere
- The pwiz version is in ./pwiz/Version.cpp, e.g. 3_0_20294
- Move old stuff away: cp -avx pwiz pwiz-old && git rm -rf pwiz

## Updating CVs

**WARNING:** Proteowizard is now .NET/CLI. Do NOT download cv.hpp, cv.cpp, or CVTranslator.*
from GitHub - you'll get empty/corrupted files.

### The updateCvParams.sh Script

Run from the repository root:

```bash
./updateCvParams.sh
```

This downloads only the `.obo` files (psi-ms.obo, unimod.obo, unit.obo) from their upstream sources.

**Do NOT download CVTranslator.hpp, CVTranslator.cpp, or CVTranslatorTest.cpp** - these are mzR-local files with conflict-resolution logic.