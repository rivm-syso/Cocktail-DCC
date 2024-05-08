cmake_minimum_required(VERSION 3.20.2 FATAL_ERROR)
project(cocktaildcc LANGUAGES Fortran)

set(CMAKE_Fortran_EXTENSIONS OFF)

add_subdirectory(external)
add_subdirectory(src)

add_custom_target(
    precompute_everything
    COMMAND
        ${CMAKE_COMMAND} -E env COCKTAIL_DCC_ICRP_SJ_DIR=${icrp_sj_SOURCE_DIR}
                                COCKTAIL_DCC_ENDF_DIR=${endf_SOURCE_DIR}
                                COCKTAIL_DCC_SOURCES_DIR=${CMAKE_SOURCE_DIR}
        $<TARGET_FILE:nuclide_decay> sparse
    COMMAND
        ${CMAKE_COMMAND} -E env COCKTAIL_DCC_ICRP_SJ_DIR=${icrp_sj_SOURCE_DIR}
                                COCKTAIL_DCC_ENDF_DIR=${endf_SOURCE_DIR}
                                COCKTAIL_DCC_SOURCES_DIR=${CMAKE_SOURCE_DIR}
        $<TARGET_FILE:test_cocktail_dcc>
    DEPENDS nuclide_decay test_cocktail_dcc
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
