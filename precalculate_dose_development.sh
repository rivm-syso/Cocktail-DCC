#!/usr/bin/env bash
set -e

./build.sh

./build/src/NuclideDecay03c.exe sparse
./build/src/testcocktailDCC02
...
