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
