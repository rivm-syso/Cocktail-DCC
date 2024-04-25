#!/usr/bin/env bash
set -e

build_type=release

rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE="${build_type}"
cmake --build build --config "${build_type}"
