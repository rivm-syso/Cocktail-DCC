#!/usr/bin/env bash
set -e

rm -rf build
cmake -S . -B build
cmake --build build

export COCKTAIL_DCC_ICRP_SJ_DIR='build/icrp_sj-src'

#! if output needs to be somewhere, that could be added in the nuclide_decay.f90 file at line 89.
./build/src/nuclide_decay sparse

#! needs the sj-zip-2-ani-49-2.zip unzipped
./build/src/test_cocktail_dcc

#In libpinpoint ln 68 the location to the matrices is to be set
#in libdcc ln 195-199 is where the values for DCC values is given
#in libdcc ln 116-125 is where the different locations for ground dcc can be given.