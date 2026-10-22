# Regenerating CV Files

## Files

- **OBO sources:** `psi-ms.obo`, `unimod.obo`, `unit.obo` (downloaded via updateCvParams.sh)
- **Generated:** `cv.hpp`, `cv.cpp` (regenerated from .obo files using cvgen)
- **mzR-local:** `CVTranslator.hpp`, `CVTranslator.cpp`, `CVTranslatorTest.cpp` (never download/regenerate)

## CRITICAL: boost::utility::singleton.hpp

**cvgen requires `boost::utility/singleton.hpp` which was removed from Boost 1.70+.**

### Workaround

A bundled version exists at `src/boost/utility/singleton.hpp`. Ensure your boost include path
prioritizes this file.

## Update Process

1. **Download .obo files:**
   ```bash
   cd src
   ./updateCvParams.sh
   ```

2. **Check version compatibility:**
   ```bash
   grep "data-version:" src/pwiz/data/common/psi-ms.obo
   grep "data-version:" src/pwiz/data/common/cv.hpp
   ```
   If they match, skip to step 5.

3. **(Advanced) Regenerate cv.hpp/cv.cpp:**
   ```bash
   cd src/pwiz/data/common
   mkdir -p build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j4 cvgen
   ./cvgen ../psi-ms.obo ../unimod.obo ../unit.obo
   ```

4. **Update CVTranslator.cpp if needed:**
   After regeneration, update `shouldIgnore()` and `shouldReplace()` functions
   to handle any new conflicts.

5. **Verify:**
   ```bash
   R CMD INSTALL .
   Rscript -e "library(mzR); packageVersion('mzR')"
   ```

## CI/CD Note

CI systems should assume cv.hpp, cv.cpp, and CVTranslator.* are already present and correct.
Only update .obo files if regeneration infrastructure is available.