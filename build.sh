#!/usr/bin/env bash
set -e

rm -rf build
cmake -S . -B build
cmake --build build
