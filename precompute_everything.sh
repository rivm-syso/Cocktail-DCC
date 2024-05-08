#!/usr/bin/env bash
set -e

build_dir="${PWD}/build"
build_type=RelWithDebInfo

rm -rf "${build_dir}"
rm -rf matrices
cmake -S . -B "${build_dir}" -DCMAKE_BUILD_TYPE="${build_type}"
cmake --build "${build_dir}" --config "${build_type}"
cmake --build "${build_dir}" --config "${build_type}" --target precompute_everything

#In libpinpoint ln 68 the location to the matrices is to be set
#in libdcc ln 195-199 is where the values for DCC values is given
#in libdcc ln 116-125 is where the different locations for ground dcc can be given.
