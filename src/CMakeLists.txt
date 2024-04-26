# The library, shared by the excectables:
file(GLOB decay_source_files CONFIGURE_DEPENDS lib/*.f90)
add_library(decay ${decay_source_files})
set(PUFF_PROGRAMS
   nuclide_decay
   source_term_dose
   test_cocktail_dcc
   mature_nuclides
   )
target_link_libraries(decay matrixexponential)

# The excecutables:
foreach(prog ${PUFF_PROGRAMS})
   cmake_path(GET prog FILENAME filename)
   add_executable("${prog}" "${prog}.f90")
   target_link_libraries("${prog}" decay)
endforeach()
