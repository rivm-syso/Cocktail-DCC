# The library, shared by the excectables:
file(GLOB decay_source_files CONFIGURE_DEPENDS lib/*.f90)
add_library(decay ${decay_source_files})
target_link_libraries(decay matrixexponential)

# The excecutables:
foreach(_prog    nuclide_decay source_term_dose test_cocktail_dcc mature_nuclides)
   add_executable("${_prog}" "${_prog}.f90")
   target_link_libraries("${_prog}" decay)
endforeach()
