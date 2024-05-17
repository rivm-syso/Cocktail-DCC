add_library(matrixexponential dmexp.f90)
target_include_directories(matrixexponential PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

include(FetchContent)

# Work around a known CMake problem, where `cmake -E tar xzf foo.zip` fails
# for some zip files on some file systems (like on AlmaLinux 9).
#   https://stackoverflow.com/questions/61488479/cmake-error-problem-with-archive-write-finish-entry-cant-restore-time
#
# Therefore use `DOWNLOAD_NO_EXTRACT`, and extract and apply patches with separate commands:
#
FetchContent_Declare(
  _icrp_sj_unextracted
  URL https://www.icrp.org/docs/sj-zip-2-ani-49-2.zip
  URL_HASH SHA512=ac0a330fa903c53e1c151825bdc01bef0c93a04fd998b0e0ce3acee769dc4adbd5a35be2aabe80c14002862b1e4947da95ff9727ba7ac6f67fa965e16bca8309
  DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_MakeAvailable(_icrp_sj_unextracted)
set(icrp_sj_SOURCE_DIR ${CMAKE_BINARY_DIR}/icrp_sj)
file(MAKE_DIRECTORY ${icrp_sj_SOURCE_DIR})
execute_process(
  COMMAND unzip -qq -d ${icrp_sj_SOURCE_DIR} -- ${_icrp_sj_unextracted_SOURCE_DIR}/sj-zip-2-ani-49-2.zip
  COMMAND_ERROR_IS_FATAL ANY)
set(icrp_sj_SOURCE_DIR "${icrp_sj_SOURCE_DIR}/Supplemental Files -v4" PARENT_SCOPE)


set(endf_SOURCE_DIR ${CMAKE_BINARY_DIR}/endf)
FetchContent_Declare(
  _endf_unextracted
  URL https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_decay.zip
  URL_HASH SHA512=bc13bd87f3a91ebd4bb0897ae4b2c4dcac7ed0d3dd63b8b2e11e49aad34e4badb6f43ad98dcb7e4f690ee6c69afdafa9880a6673e1d76d7d7547585bd7452b75
  DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_MakeAvailable(_endf_unextracted)
file(MAKE_DIRECTORY ${endf_SOURCE_DIR})
execute_process(
  COMMAND unzip -qq -d ${endf_SOURCE_DIR} -- ${_endf_unextracted_SOURCE_DIR}/ENDF-B-VIII.0_decay.zip
  COMMAND_ERROR_IS_FATAL ANY)
set(endf_SOURCE_DIR "${endf_SOURCE_DIR}/ENDF-B-VIII.0_decay" PARENT_SCOPE)
execute_process(
  COMMAND sh -c "< ${CMAKE_CURRENT_LIST_DIR}/endf_viii_rename_cs.patch patch"
  WORKING_DIRECTORY ${endf_SOURCE_DIR}/ENDF-B-VIII.0_decay
  COMMAND_ERROR_IS_FATAL ANY)
