add_library(matrixexponential dmexp.f90)
target_include_directories(matrixexponential PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
