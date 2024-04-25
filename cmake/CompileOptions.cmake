# TODO:  turn on all sensible warnings and runtime checks,
#        and fix the errors one by one
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "Ref. the CMake docs." FORCE)
endif()
message("* Current build type is : ${CMAKE_BUILD_TYPE}")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
   set(basic_flags
      -g
      -fbacktrace
      # -std=f2018
      -pedantic
      # -Werror
      -Wall
      -Wextra
      -Werror=shadow
      -Werror=intrinsic-shadow
      -Wno-uninitialized
      -Warray-bounds
      -Wunreachable-code
      -Wconversion
      -Wstrict-overflow=5
      -Wunused-dummy-argument

      -Wno-maybe-uninitialized  # because it complains about auto-allocatation in Fortran 2003+
      -Wno-aggregate-return
      -Wno-compare-reals  # it can sometimes be valid to compare to 0 to check whether a step can be skipped
      -Waliasing
      -Wampersand
      -Wc-binding-type
      -Wcharacter-truncation
      -Wconversion
      -Wdo-subscript
      -Wimplicit-interface
      -Wimplicit-procedure
      -Wintrinsic-shadow
      -Wintrinsics-std
      -Wline-truncation
      -Wno-tabs
      -Wreal-q-constant
      -Wsurprising
      -Wunderflow
      -Wunused-parameter
      -Wfrontend-loop-interchange
      -Wtarget-lifetime
   )
   set(debug_flags
      -O0
      -ftrapv
      -fsanitize-address-use-after-scope
      -ffpe-trap=invalid,zero,overflow
      -fcheck=all
   )
   set(release_flags
      -O2
      -march=native
      -faggressive-function-elimination
      -Wno-error=function-elimination
      -fno-protect-parens
      -Wno-strict-overflow
      -Wno-array-bounds
      -Wno-uninitialized
   )
else ()
   message(FATAL_ERROR "Unknown CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
endif ()

foreach(opt ${basic_flags})
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${opt}")
endforeach()
foreach(opt ${debug_flags})
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${opt}")
endforeach()
foreach(opt ${release_flags})
   set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${opt}")
   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} ${opt}")
endforeach()
