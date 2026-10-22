#!/bin/bash
# ============================================================================
# Update CV (Controlled Vocabulary) files for mzR
# ============================================================================
#
# IMPORTANT: The Proteowizard project has been restructured to .NET/CLI.
# The C++ files (cv.hpp, cv.cpp, CVTranslator.*) are NO LONGER available
# from https://github.com/ProteoWizard/pwiz. Attempting to download them
# will result in empty or corrupted files.
#
# This script ONLY downloads the .obo ontology files from their upstream
# sources. The generated C++ files (cv.hpp, cv.cpp) and the mzR-specific
# CVTranslator files are maintained locally in the repository.
# ============================================================================

set -e  # Exit on error

# Ensure we're in the right directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p src/pwiz/data/common

echo "=== Downloading PSI-MS ontology (latest) ==="
wget -N -O src/pwiz/data/common/psi-ms.obo \
  https://github.com/HUPO-PSI/psi-ms-CV/raw/master/psi-ms.obo

echo "=== Downloading UniMod ontology (latest) ==="
wget -N -O src/pwiz/data/common/unimod.obo \
  https://www.unimod.org/obo/unimod.obo

echo "=== Downloading Unit Ontology (UO) (latest) ==="
wget -N -O src/pwiz/data/common/unit.obo \
  https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/uo.obo

echo ""
echo "=== CV files updated successfully ==="
echo ""
echo "NOTES:"
echo "- cv.hpp and cv.cpp are generated files that must be regenerated"
echo "  using cvgen (see src/pwiz/data/common/REGENERATE_CV.md)"
echo "- CVTranslator.hpp, CVTranslator.cpp, CVTranslatorTest.cpp are"
echo "  mzR-local files that should NOT be downloaded or regenerated."
echo "- If you need to regenerate cv.hpp/cv.cpp, see the documentation"