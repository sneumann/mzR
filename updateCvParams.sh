#!/bin/bash
## Regenerate src/pwiz/data/common/cv.{hpp,cpp} from the latest CV releases.
##
## ProteoWizard retired its C++ tree (ProteoWizard/pwiz#4658, commit 6dd40d8c93),
## so cv.cpp and cvgen are no longer available on pwiz master. We pin the C++
## generator (cvgen) to the last commit before that change and run it on the
## latest upstream .obo files.
##
## Requires: a C++17 compiler and Boost (filesystem, xpressive), e.g. via Homebrew.
set -ex

PWIZ_COMMIT=a9ed4ce96066fccc0ffbcb2ce3e97c3d5eac65dc
PWIZ_PREFIX=https://raw.githubusercontent.com/ProteoWizard/pwiz/$PWIZ_COMMIT/pwiz/data/common
BOOST_PREFIX=${BOOST_PREFIX:-$(brew --prefix 2>/dev/null || echo /usr/local)}

ROOT=$(cd "$(dirname "$0")" && pwd)
CVDIR=$ROOT/src/pwiz/data/common
WORK=$(mktemp -d "${TMPDIR:-/tmp}/cvgen.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

## CV releases (order matters: psi-ms, unimod, unit)
wget -O $CVDIR/psi-ms.obo https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
wget -O $CVDIR/unimod.obo https://www.unimod.org/obo/unimod.obo
wget -O $CVDIR/unit.obo   https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/uo.obo

## Generator sources from the last C++ pwiz commit.
## cvgen.cpp is kept as is: mzR patches it to skip duplicate synonym enum names.
for F in cvgen_common.hpp obo.cpp obo.hpp cv.inl CVTranslator.cpp CVTranslator.hpp ; do
  wget -O $CVDIR/$F $PWIZ_PREFIX/$F
done

## Build cvgen and regenerate cv.hpp / cv.cpp (written next to the executable)
${CXX:-c++} -std=c++17 -O1 -w \
  -I$BOOST_PREFIX/include -I$ROOT/src -I$ROOT/src/boost_aux \
  -o $WORK/cvgen $CVDIR/cvgen.cpp $CVDIR/obo.cpp \
  -L$BOOST_PREFIX/lib -lboost_filesystem
$WORK/cvgen $CVDIR/psi-ms.obo $CVDIR/unimod.obo $CVDIR/unit.obo
cp $WORK/cv.hpp $WORK/cv.cpp $CVDIR/
