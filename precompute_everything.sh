#!/usr/bin/env bash
set -e

build_dir="${PWD}/build"
build_type=RelWithDebInfo

rm -rf "${build_dir}"
rm -rf matrices
cmake -S . -B "${build_dir}" -DCMAKE_BUILD_TYPE="${build_type}"
cmake --build "${build_dir}" --config "${build_type}"
cmake --build "${build_dir}" --config "${build_type}" --target precompute_everything
