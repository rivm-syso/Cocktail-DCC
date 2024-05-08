add_library(matrixexponential dmexp.f90)
target_include_directories(matrixexponential PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

include(FetchContent)

FetchContent_Declare(
  icrp_sj
  URL https://www.icrp.org/docs/sj-zip-2-ani-49-2.zip
  URL_HASH SHA512=ac0a330fa903c53e1c151825bdc01bef0c93a04fd998b0e0ce3acee769dc4adbd5a35be2aabe80c14002862b1e4947da95ff9727ba7ac6f67fa965e16bca8309
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(icrp_sj)

FetchContent_Declare(
  endf
  URL https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_decay.zip
  URL_HASH SHA512=bc13bd87f3a91ebd4bb0897ae4b2c4dcac7ed0d3dd63b8b2e11e49aad34e4badb6f43ad98dcb7e4f690ee6c69afdafa9880a6673e1d76d7d7547585bd7452b75
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  PATCH_COMMAND < "${CMAKE_CURRENT_LIST_DIR}/endf_viii_rename_cs.patch" patch
  UPDATE_DISCONNECTED 1
)
FetchContent_MakeAvailable(endf)
