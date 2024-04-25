cmake_minimum_required(VERSION 3.20.2 FATAL_ERROR)
project(cocktaildcc LANGUAGES Fortran)

set(CMAKE_Fortran_EXTENSIONS OFF)

add_subdirectory(external)

include(cmake/CompileOptions.cmake)

add_subdirectory(src)
