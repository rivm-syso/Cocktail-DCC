cmake_minimum_required(VERSION 3.20.2 FATAL_ERROR)
project(cocktaildcc LANGUAGES Fortran)

set(CMAKE_Fortran_EXTENSIONS OFF)

add_subdirectory(external)

include(cmake/CompileOptions.cmake)

add_subdirectory(src)

set(output_dir "$(CMAKE_SOURCE_DIR)/matrices")
add_custom_target(
    precompute_everything
    COMMAND
        $(CMAKE_COMMAND) -E make_directory ${output_dir}
    COMMAND
        ${CMAKE_COMMAND} -E env COCKTAIL_DCC_ICRP_SJ_DIR=${icrp_sj_SOURCE_DIR}
                                COCKTAIL_DCC_ENDF_DIR=${endf_SOURCE_DIR}
                                COCKTAIL_DCC_SOURCES_DIR=${CMAKE_SOURCE_DIR}
                                COCKTAIL_DCC_OUTPUT_DIR=${output_dir}
        $<TARGET_FILE:nuclide_decay> sparse
    COMMAND
        ${CMAKE_COMMAND} -E env COCKTAIL_DCC_ICRP_SJ_DIR=${icrp_sj_SOURCE_DIR}
                                COCKTAIL_DCC_ENDF_DIR=${endf_SOURCE_DIR}
                                COCKTAIL_DCC_SOURCES_DIR=${CMAKE_SOURCE_DIR}
                                COCKTAIL_DCC_OUTPUT_DIR=${output_dir}
        $<TARGET_FILE:test_cocktail_dcc>
    COMMAND
        # It is silly that we need to move these files manually. The Fortran code should not have placed
        # in an input directory in the first place, but that would require some non-trivial refactoring.
        mv ${CMAKE_SOURCE_DIR}/resources/*_Pinpoint*.txt ${output_dir}
    DEPENDS nuclide_decay test_cocktail_dcc
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
