add_library(matrixexponential dmexp.f90)
target_include_directories(matrixexponential PUBLIC ${CMAKE_CURRENT_BINARY_DIR})


include(FetchContent)
set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR})

FetchContent_Declare(
  icrp_sj
  URL      https://www.icrp.org/docs/sj-zip-2-ani-49-2.zip
  URL_HASH SHA512=ac0a330fa903c53e1c151825bdc01bef0c93a04fd998b0e0ce3acee769dc4adbd5a35be2aabe80c14002862b1e4947da95ff9727ba7ac6f67fa965e16bca8309
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(icrp_sj)
